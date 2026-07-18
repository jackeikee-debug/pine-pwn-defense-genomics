#!/usr/bin/env python3
"""Summarize route 1 P. densiflora Trinity blastx hits against current Pd proteins."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


FIELDS = [
    "transcript_id",
    "orthogroup_ids",
    "effector_target_network_seed_symbols",
    "sequence_count",
    "extracted_sequence_ids",
    "mapping_status",
    "top_pd_hit",
    "top_hit_bitscore",
    "top_hit_pident",
    "top_hit_evalue",
    "top_hit_query_aa_coverage",
    "best_effector_target_network_candidate_pd_gene",
    "best_effector_target_network_candidate_bitscore",
    "best_effector_target_network_candidate_pident",
    "best_effector_target_network_candidate_evalue",
    "effector_target_network_candidate_pd_genes",
    "claim_ceiling",
]

CLAIM_CEILING = (
    "blastx sequence-level bridge from Trinity transcript to current Pd protein; "
    "not differential-expression evidence and not experimental validation"
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_effector_target_network_pden_trinity_to_pd_mapping(
    blastx_path: Path,
    sequence_manifest_path: Path,
    network_seed_path: Path,
    output_path: Path,
    markdown_path: Path,
    min_bitscore: float,
    min_pident: float,
) -> list[dict[str, str]]:
    candidate_genes = effector_target_network_candidate_genes(read_tsv(network_seed_path))
    hits = read_blastx_hits(blastx_path, min_bitscore=min_bitscore, min_pident=min_pident)
    rows = []
    for item in read_tsv(sequence_manifest_path):
        rows.append(build_mapping_row(item, hits.get(item["transcript_id"], []), candidate_genes))
    write_tsv(output_path, rows)
    write_markdown(markdown_path, rows, min_bitscore=min_bitscore, min_pident=min_pident)
    return rows


def effector_target_network_candidate_genes(seed_rows: list[dict[str, str]]) -> dict[str, set[str]]:
    grouped: dict[str, set[str]] = defaultdict(set)
    for row in seed_rows:
        orthogroup_id = row.get("orthogroup_id", "")
        for gene_id in split_semicolon(row.get("pine_gene_ids", "")):
            if gene_id.startswith("pden|"):
                grouped[orthogroup_id].add(gene_id.split("|", 1)[1])
    return grouped


def read_blastx_hits(path: Path, min_bitscore: float, min_pident: float) -> dict[str, list[dict[str, str]]]:
    hits: dict[str, list[dict[str, str]]] = defaultdict(list)
    if not path.exists():
        return hits
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            qseqid, sseqid, pident, length, evalue, bitscore, qlen, slen, qstart, qend, sstart, send = fields[:12]
            if float(bitscore) < min_bitscore or float(pident) < min_pident:
                continue
            base_id = trinity_base_id(qseqid)
            hit = {
                "qseqid": qseqid,
                "sseqid": sseqid,
                "pident": format_float(float(pident)),
                "length": length,
                "evalue": evalue,
                "bitscore": format_float(float(bitscore)),
                "qlen": qlen,
                "slen": slen,
                "query_aa_coverage": query_aa_coverage(float(length), float(qlen)),
            }
            hits[base_id].append(hit)
    for base_id in hits:
        hits[base_id].sort(key=lambda row: float(row["bitscore"]), reverse=True)
    return hits


def build_mapping_row(
    manifest_row: dict[str, str],
    hits: list[dict[str, str]],
    candidate_genes: dict[str, set[str]],
) -> dict[str, str]:
    orthogroup_ids = split_semicolon(manifest_row.get("orthogroup_ids", ""))
    effector_target_network_candidates = sorted({gene for og in orthogroup_ids for gene in candidate_genes.get(og, set())})
    candidate_hits = [hit for hit in hits if hit["sseqid"] in effector_target_network_candidates]
    top_hit = hits[0] if hits else {}
    best_candidate = candidate_hits[0] if candidate_hits else {}
    return {
        "transcript_id": manifest_row.get("transcript_id", ""),
        "orthogroup_ids": manifest_row.get("orthogroup_ids", ""),
        "effector_target_network_seed_symbols": manifest_row.get("effector_target_network_seed_symbols", ""),
        "sequence_count": manifest_row.get("sequence_count", ""),
        "extracted_sequence_ids": manifest_row.get("extracted_sequence_ids", ""),
        "mapping_status": classify_mapping_status(hits, candidate_hits),
        "top_pd_hit": top_hit.get("sseqid", ""),
        "top_hit_bitscore": top_hit.get("bitscore", ""),
        "top_hit_pident": top_hit.get("pident", ""),
        "top_hit_evalue": top_hit.get("evalue", ""),
        "top_hit_query_aa_coverage": top_hit.get("query_aa_coverage", ""),
        "best_effector_target_network_candidate_pd_gene": best_candidate.get("sseqid", ""),
        "best_effector_target_network_candidate_bitscore": best_candidate.get("bitscore", ""),
        "best_effector_target_network_candidate_pident": best_candidate.get("pident", ""),
        "best_effector_target_network_candidate_evalue": best_candidate.get("evalue", ""),
        "effector_target_network_candidate_pd_genes": ";".join(effector_target_network_candidates),
        "claim_ceiling": CLAIM_CEILING,
    }


def classify_mapping_status(hits: list[dict[str, str]], candidate_hits: list[dict[str, str]]) -> str:
    if candidate_hits:
        return "current_pd_effector_target_network_candidate_hit"
    if hits:
        return "non_effector_target_network_pd_hit_only"
    return "no_blastx_hit"


def query_aa_coverage(alignment_length: float, query_nt_length: float) -> str:
    if query_nt_length <= 0:
        return ""
    return format_float(alignment_length / (query_nt_length / 3.0))


def trinity_base_id(transcript_id: str) -> str:
    match = re.match(r"(.+)_i\d+$", transcript_id)
    if match:
        return match.group(1)
    return transcript_id


def split_semicolon(value: str) -> set[str]:
    return {item.strip() for item in value.split(";") if item.strip()}


def format_float(value: float) -> str:
    return f"{value:.6g}"


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(path: Path, rows: list[dict[str, str]], min_bitscore: float, min_pident: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    counts = Counter(row["mapping_status"] for row in rows)
    candidate_rows = [row for row in rows if row["mapping_status"] == "current_pd_effector_target_network_candidate_hit"]
    lines = [
        "# Route 1 P. densiflora Trinity-to-Pd Mapping",
        "",
        "This report summarizes DIAMOND blastx mappings from route 1 Trinity bridge transcripts to the current P. densiflora HA protein set.",
        "",
        f"- Transcript targets summarized: {len(rows)}",
        f"- Current Pd effector_target_network candidate hits: {counts.get('current_pd_effector_target_network_candidate_hit', 0)}",
        f"- Non-effector_target_network Pd hits only: {counts.get('non_effector_target_network_pd_hit_only', 0)}",
        f"- No blastx hit: {counts.get('no_blastx_hit', 0)}",
        f"- Filtering: bitscore >= {min_bitscore:g}, pident >= {min_pident:g}",
        "",
        "Evidence boundary: sequence-level mapping is not differential-expression support.",
        "",
    ]
    if candidate_rows:
        lines.extend(["## Candidate Gene Hits", ""])
    for row in candidate_rows:
        lines.append(
            f"- {row['orthogroup_ids']} / {row['transcript_id']} -> "
            f"{row['best_effector_target_network_candidate_pd_gene']} "
            f"(bitscore {row['best_effector_target_network_candidate_bitscore']}, pident {row['best_effector_target_network_candidate_pident']})"
        )
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blastx", type=Path, required=True)
    parser.add_argument("--sequence-manifest", type=Path, required=True)
    parser.add_argument("--network-seed", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    parser.add_argument("--min-bitscore", type=float, default=80.0)
    parser.add_argument("--min-pident", type=float, default=35.0)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = build_effector_target_network_pden_trinity_to_pd_mapping(
        blastx_path=args.blastx,
        sequence_manifest_path=args.sequence_manifest,
        network_seed_path=args.network_seed,
        output_path=args.output,
        markdown_path=args.markdown,
        min_bitscore=args.min_bitscore,
        min_pident=args.min_pident,
    )
    counts = Counter(row["mapping_status"] for row in rows)
    print(
        "Route 1 P. densiflora Trinity-to-Pd rows written: "
        f"{len(rows)}; candidate hits: {counts.get('current_pd_effector_target_network_candidate_hit', 0)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
