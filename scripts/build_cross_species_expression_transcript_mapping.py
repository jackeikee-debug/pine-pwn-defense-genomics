#!/usr/bin/env python3
"""Filter targeted DIAMOND hits into auditable transcript-to-OG mappings."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


FIELDS = [
    "transcript_id", "target_gene_id", "orthogroup_id", "target_species_id",
    "mapping_status", "pident", "alignment_length_aa", "query_translation_coverage",
    "subject_coverage", "evalue", "bitscore", "evidence_type",
]


def _read_manifest(path: Path, species: str) -> dict[str, str]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return {
            row["gene_id"]: row["orthogroup_id"]
            for row in csv.DictReader(handle, delimiter="\t")
            if row["species_id"] == species
        }


def _hit(values: list[str], target_to_og: dict[str, str]) -> dict[str, object]:
    if len(values) != 10:
        raise ValueError(f"Expected 10 DIAMOND columns, found {len(values)}")
    query, subject, pident, length, qlen, _qstart, _qend, slen, evalue, bitscore = values
    if subject not in target_to_og:
        raise ValueError(f"DIAMOND subject absent from target manifest: {subject}")
    length_n, qlen_n, slen_n = float(length), float(qlen), float(slen)
    query_cov = 3 * length_n / qlen_n if qlen_n else 0.0
    return {
        "query": query, "subject": subject, "orthogroup": target_to_og[subject],
        "pident": float(pident), "length": length_n, "query_cov": query_cov,
        "subject_cov": length_n / slen_n if slen_n else 0.0,
        "evalue": float(evalue), "bitscore": float(bitscore),
    }


def _passes(hit: dict[str, object]) -> bool:
    return bool(hit["pident"] >= 35 and hit["query_cov"] >= 0.30 and hit["evalue"] <= 1e-10 and hit["bitscore"] >= 80)


def _rank(hit: dict[str, object]) -> tuple[float, float, float, float, str]:
    return (-hit["bitscore"], hit["evalue"], -hit["pident"], -hit["length"], hit["subject"])


def _support(hit: dict[str, object]) -> tuple[float, float, float, float]:
    return (hit["bitscore"], hit["evalue"], hit["pident"], hit["length"])


def _format(value: float) -> str:
    return f"{value:.12g}"


def build_mapping(
    blast_path: Path,
    manifest_path: Path,
    species: str,
    evidence_type: str,
    output_path: Path,
) -> list[dict[str, str]]:
    target_to_og = _read_manifest(manifest_path, species)
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    with Path(blast_path).open(encoding="utf-8", newline="") as handle:
        for values in csv.reader(handle, delimiter="\t"):
            hit = _hit(values, target_to_og)
            grouped[hit["query"]].append(hit)

    rows: list[dict[str, str]] = []
    for query in sorted(grouped):
        hits = grouped[query]
        passing = [hit for hit in hits if _passes(hit)]
        considered = passing or hits
        best = sorted(considered, key=_rank)[0]
        top_ogs = {hit["orthogroup"] for hit in considered if _support(hit) == _support(best)}
        status = "threshold_failed" if not passing else "accepted"
        orthogroup = best["orthogroup"] if passing else ""
        target = best["subject"]
        if passing and len(top_ogs) > 1:
            status, orthogroup, target = "ambiguous_cross_orthogroup", "", ""
        rows.append({
            "transcript_id": query,
            "target_gene_id": target,
            "orthogroup_id": orthogroup,
            "target_species_id": species,
            "mapping_status": status,
            "pident": _format(best["pident"]),
            "alignment_length_aa": _format(best["length"]),
            "query_translation_coverage": _format(best["query_cov"]),
            "subject_coverage": _format(best["subject_cov"]),
            "evalue": _format(best["evalue"]),
            "bitscore": _format(best["bitscore"]),
            "evidence_type": evidence_type,
        })

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with Path(output_path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blast", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--species", choices=["pden", "ptae"], required=True)
    parser.add_argument("--evidence-type", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = build_mapping(args.blast, args.manifest, args.species, args.evidence_type, args.output)
    counts = {status: sum(row["mapping_status"] == status for row in rows) for status in ("accepted", "ambiguous_cross_orthogroup", "threshold_failed")}
    print(f"Wrote {len(rows)} transcript mappings: {counts}")


if __name__ == "__main__":
    main()
