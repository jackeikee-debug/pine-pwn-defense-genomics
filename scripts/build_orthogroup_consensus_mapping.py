#!/usr/bin/env python3
"""Classify transcript mappings against a multi-species orthogroup reference."""

from __future__ import annotations

import argparse
import csv
import logging
from collections import defaultdict
from pathlib import Path


FIELDS = [
    "transcript_id",
    "orthogroup_id",
    "best_orthogroup_id",
    "second_best_orthogroup_id",
    "best_target_gene_id",
    "mapping_status",
    "supporting_species_count",
    "supporting_species",
    "best_bitscore",
    "second_best_bitscore",
    "relative_bitscore_margin",
    "pident",
    "alignment_length_aa",
    "query_translation_coverage",
    "subject_coverage",
    "evalue",
    "evidence_type",
]
EVIDENCE_TYPE = "multispecies_orthogroup_consensus_homology"


def _read_manifest(path: Path) -> dict[str, dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    manifest = {row["gene_id"]: row for row in rows}
    if len(manifest) != len(rows):
        raise ValueError("Duplicate gene_id in reference manifest")
    return manifest


def _format(value: float) -> str:
    return f"{value:.12g}"


def _parse_hit(values: list[str], manifest: dict[str, dict[str, str]]) -> dict[str, object]:
    if len(values) != 10:
        raise ValueError(f"Expected 10 DIAMOND columns, found {len(values)}")
    query, subject, pident, length, qlen, _qstart, _qend, slen, evalue, bitscore = values
    if subject not in manifest:
        raise ValueError(f"DIAMOND subject absent from manifest: {subject}")
    length_n = float(length)
    qlen_n = float(qlen)
    slen_n = float(slen)
    return {
        "query": query,
        "subject": subject,
        "orthogroup": manifest[subject]["orthogroup_id"],
        "species": manifest[subject]["species_id"],
        "pident": float(pident),
        "length": length_n,
        "query_cov": 3 * length_n / qlen_n if qlen_n else 0.0,
        "subject_cov": length_n / slen_n if slen_n else 0.0,
        "evalue": float(evalue),
        "bitscore": float(bitscore),
    }


def _passes(hit: dict[str, object]) -> bool:
    return bool(
        hit["pident"] >= 35
        and hit["query_cov"] >= 0.30
        and hit["evalue"] <= 1e-10
        and hit["bitscore"] >= 80
    )


def _rank(hit: dict[str, object]) -> tuple[float, float, float, float, str]:
    return (-hit["bitscore"], hit["evalue"], -hit["pident"], -hit["length"], hit["subject"])


def build_mapping(blast_path: Path, manifest_path: Path, output_path: Path) -> list[dict[str, str]]:
    """Assign transcripts only when a best OG has multi-species support and separation."""
    manifest = _read_manifest(manifest_path)
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    with Path(blast_path).open(encoding="utf-8", newline="") as handle:
        for values in csv.reader(handle, delimiter="\t"):
            hit = _parse_hit(values, manifest)
            grouped[hit["query"]].append(hit)

    rows: list[dict[str, str]] = []
    for transcript in sorted(grouped):
        all_hits = grouped[transcript]
        passing = [hit for hit in all_hits if _passes(hit)]
        considered = passing or all_hits
        best_overall = sorted(considered, key=_rank)[0]
        by_og: dict[str, list[dict[str, object]]] = defaultdict(list)
        for hit in passing:
            by_og[hit["orthogroup"]].append(hit)
        og_best = sorted(
            ((orthogroup, sorted(hits, key=_rank)[0]) for orthogroup, hits in by_og.items()),
            key=lambda item: _rank(item[1]),
        )
        best_og = og_best[0][0] if og_best else ""
        best_hit = og_best[0][1] if og_best else best_overall
        second_og = og_best[1][0] if len(og_best) > 1 else ""
        second_score = float(og_best[1][1]["bitscore"]) if len(og_best) > 1 else 0.0
        best_score = float(best_hit["bitscore"])
        margin = (best_score - second_score) / best_score if best_score and second_og else 1.0
        species = sorted({str(hit["species"]) for hit in by_og.get(best_og, [])})
        if not passing:
            status = "threshold_failed"
        elif margin < 0.10:
            status = "ambiguous_orthogroup"
        elif len(species) < 2:
            status = "weak_species_support"
        else:
            status = "accepted_consensus"
        rows.append(
            {
                "transcript_id": transcript,
                "orthogroup_id": best_og if status == "accepted_consensus" else "",
                "best_orthogroup_id": best_og,
                "second_best_orthogroup_id": second_og,
                "best_target_gene_id": str(best_hit["subject"]),
                "mapping_status": status,
                "supporting_species_count": str(len(species)),
                "supporting_species": ";".join(species),
                "best_bitscore": _format(best_score),
                "second_best_bitscore": _format(second_score) if second_og else "",
                "relative_bitscore_margin": _format(margin),
                "pident": _format(float(best_hit["pident"])),
                "alignment_length_aa": _format(float(best_hit["length"])),
                "query_translation_coverage": _format(float(best_hit["query_cov"])),
                "subject_coverage": _format(float(best_hit["subject_cov"])),
                "evalue": _format(float(best_hit["evalue"])),
                "evidence_type": EVIDENCE_TYPE,
            }
        )

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with Path(output_path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Wrote %d consensus mapping audit rows", len(rows))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blast", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    build_mapping(args.blast, args.manifest, args.output)


if __name__ == "__main__":
    main()
