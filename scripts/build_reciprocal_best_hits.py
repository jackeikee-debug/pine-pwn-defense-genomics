#!/usr/bin/env python3
"""Build a conservative reciprocal-best-hit crosswalk from pine genes to TAIR loci."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path


DIAMOND_FIELDS = ["qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore"]
OUTPUT_FIELDS = [
    "pmas_gene_id",
    "uniprot_accession",
    "tair_locus",
    "identity",
    "query_coverage",
    "subject_coverage",
    "evalue",
    "bitscore",
    "relative_bitscore_margin",
    "mapping_status",
]


def uniprot_accession(identifier: str) -> str:
    parts = identifier.split("|")
    return parts[1] if len(parts) >= 3 and parts[0] in {"sp", "tr"} else identifier


def pmas_gene_id(identifier: str) -> str:
    return identifier.split("|", 1)[1] if identifier.startswith("pmas|") else identifier


def safe_float(value: str, default: float = math.nan) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def hit_metrics(hit: dict[str, str]) -> dict[str, float]:
    aligned = safe_float(hit.get("length", ""), 0.0)
    qlen = safe_float(hit.get("qlen", ""), 0.0)
    slen = safe_float(hit.get("slen", ""), 0.0)
    return {
        "identity": safe_float(hit.get("pident", ""), 0.0),
        "query_coverage": 100.0 * aligned / qlen if qlen > 0 else 0.0,
        "subject_coverage": 100.0 * aligned / slen if slen > 0 else 0.0,
        "evalue": safe_float(hit.get("evalue", ""), math.inf),
        "bitscore": safe_float(hit.get("bitscore", ""), 0.0),
    }


def passes_thresholds(hit: dict[str, str]) -> bool:
    metrics = hit_metrics(hit)
    return (
        metrics["identity"] >= 30.0
        and metrics["query_coverage"] >= 50.0
        and metrics["subject_coverage"] >= 50.0
        and metrics["evalue"] <= 1e-10
    )


def unique_best_hits(
    hits: list[dict[str, str]],
    relative_margin: float = 0.05,
) -> dict[str, dict[str, object]]:
    by_query: dict[str, list[dict[str, str]]] = defaultdict(list)
    for hit in hits:
        by_query[hit["qseqid"]].append(hit)

    decisions: dict[str, dict[str, object]] = {}
    for query, query_hits in by_query.items():
        passing = [hit for hit in query_hits if passes_thresholds(hit)]
        if not passing:
            decisions[query] = {"status": "threshold_failed", "hit": max(query_hits, key=lambda row: safe_float(row["bitscore"], 0.0)), "margin": 0.0}
            continue
        passing.sort(key=lambda row: (-safe_float(row["bitscore"], 0.0), row["sseqid"]))
        top = passing[0]
        top_score = safe_float(top["bitscore"], 0.0)
        margin = 1.0
        if len(passing) > 1 and top_score > 0:
            margin = (top_score - safe_float(passing[1]["bitscore"], 0.0)) / top_score
        if len(passing) > 1 and margin < relative_margin:
            decisions[query] = {"status": "ambiguous_best_hit", "hit": top, "margin": margin}
        else:
            decisions[query] = {"status": "unique_best", "hit": top, "margin": margin}
    return decisions


def select_reciprocal_best_hits(
    forward_hits: list[dict[str, str]],
    reverse_hits: list[dict[str, str]],
    accession_to_tair: dict[str, str],
    all_pmas_ids: set[str] | None = None,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    forward = unique_best_hits(forward_hits)
    reverse = unique_best_hits(reverse_hits)
    audit: list[dict[str, str]] = []

    for forward_query in sorted(forward):
        decision = forward[forward_query]
        hit = decision["hit"]
        metrics = hit_metrics(hit)
        accession = uniprot_accession(hit["sseqid"])
        tair = accession_to_tair.get(accession, "").upper()
        status = str(decision["status"])
        if status == "unique_best":
            reverse_decision = reverse.get(hit["sseqid"])
            if reverse_decision is None:
                reverse_decision = reverse.get(accession)
            if reverse_decision is None or reverse_decision["status"] != "unique_best":
                status = "no_unique_reverse_best"
            elif pmas_gene_id(reverse_decision["hit"]["sseqid"]) != pmas_gene_id(forward_query):
                status = "nonreciprocal_best"
            elif not tair:
                status = "missing_tair_locus"
            else:
                status = "reciprocal_best"
        audit.append(
            {
                "pmas_gene_id": pmas_gene_id(forward_query),
                "uniprot_accession": accession,
                "tair_locus": tair,
                "identity": f"{metrics['identity']:.3f}",
                "query_coverage": f"{metrics['query_coverage']:.3f}",
                "subject_coverage": f"{metrics['subject_coverage']:.3f}",
                "evalue": f"{metrics['evalue']:.6g}",
                "bitscore": f"{metrics['bitscore']:.3f}",
                "relative_bitscore_margin": f"{float(decision['margin']):.6f}",
                "mapping_status": status,
            }
        )

    for identifier in sorted((all_pmas_ids or set()) - set(forward)):
        audit.append(
            {
                "pmas_gene_id": pmas_gene_id(identifier),
                "uniprot_accession": "",
                "tair_locus": "",
                "identity": "",
                "query_coverage": "",
                "subject_coverage": "",
                "evalue": "",
                "bitscore": "",
                "relative_bitscore_margin": "",
                "mapping_status": "no_diamond_hit",
            }
        )

    accepted_tair_counts = Counter(
        row["tair_locus"] for row in audit if row["mapping_status"] == "reciprocal_best"
    )
    for row in audit:
        if row["mapping_status"] == "reciprocal_best" and accepted_tair_counts[row["tair_locus"]] > 1:
            row["mapping_status"] = "ambiguous_tair_locus"
    accepted = [dict(row) for row in audit if row["mapping_status"] == "reciprocal_best"]
    return accepted, sorted(audit, key=lambda row: row["pmas_gene_id"])


def read_diamond(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t", fieldnames=DIAMOND_FIELDS))


def read_uniprot_tair_map(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accession = row.get("Entry", "").strip()
            ordered_locus = row.get("Gene Names (ordered locus)", "").strip().split()
            if accession and ordered_locus:
                mapping[accession] = ordered_locus[0].upper()
    return mapping


def read_fasta_ids(path: Path) -> set[str]:
    identifiers: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                identifiers.add(line[1:].split()[0])
    return identifiers


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--forward-hits", type=Path, required=True)
    parser.add_argument("--reverse-hits", type=Path, required=True)
    parser.add_argument("--uniprot-metadata", type=Path, required=True)
    parser.add_argument("--pmas-fasta", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit-output", type=Path, required=True)
    args = parser.parse_args()

    accepted, audit = select_reciprocal_best_hits(
        read_diamond(args.forward_hits),
        read_diamond(args.reverse_hits),
        read_uniprot_tair_map(args.uniprot_metadata),
        read_fasta_ids(args.pmas_fasta),
    )
    write_tsv(args.output, accepted)
    write_tsv(args.audit_output, audit)
    status_counts = Counter(row["mapping_status"] for row in audit)
    print(
        f"Accepted {len(accepted)} reciprocal TAIR mappings from {len(audit)} pine queries; "
        + "; ".join(f"{key}={value}" for key, value in sorted(status_counts.items())),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
