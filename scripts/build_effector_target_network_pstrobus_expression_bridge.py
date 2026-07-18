#!/usr/bin/env python3
"""Bridge P. strobus time-course expression to Route 1 orthogroups by homology."""

from __future__ import annotations

import argparse
import csv
import logging
import math
from pathlib import Path


DE_COLUMNS = ["gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
OUTPUT_FIELDS = [
    "pstrobus_transcript_id",
    "pstrobus_gene_id",
    "ptae_protein_id",
    "orthogroup_id",
    "effector_target_network_role",
    "host_module",
    "host_network_role",
    "pident",
    "alignment_length_aa",
    "query_translation_coverage",
    "subject_coverage_descriptive",
    "evalue",
    "bitscore",
    "two_week_log2_fold_change",
    "two_week_padj",
    "two_week_fdr_status",
    "four_week_log2_fold_change",
    "four_week_padj",
    "four_week_fdr_status",
    "bridge_evidence_type",
    "claim_ceiling",
]
CLAIM_CEILING = (
    "Sequence similarity to a P. taeda Route 1 orthogroup is a cross-species homology proxy; "
    "time-course differential expression supports infection response only and is not direct effector targeting, PPI, or causal tolerance evidence"
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_deseq2(path: Path) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for values in csv.reader(handle, delimiter="\t"):
            if len(values) != len(DE_COLUMNS):
                raise ValueError(f"Expected {len(DE_COLUMNS)} columns in {path}, found {len(values)}")
            rows[values[0]] = dict(zip(DE_COLUMNS, values))
    return rows


def number(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def fdr_status(row: dict[str, str] | None) -> str:
    padj = number(row.get("padj", "")) if row else math.nan
    if not math.isfinite(padj):
        return "not_tested_or_filtered"
    return "fdr_lt_0.05" if padj < 0.05 else "fdr_ge_0.05"


def format_value(value: float) -> str:
    return f"{value:.10g}"


def best_hits(blast_path: Path) -> dict[str, dict[str, str]]:
    candidates: dict[str, list[dict[str, str]]] = {}
    with blast_path.open("r", encoding="utf-8", newline="") as handle:
        for values in csv.reader(handle, delimiter="\t"):
            if len(values) != 10:
                raise ValueError(f"Expected 10 DIAMOND columns in {blast_path}, found {len(values)}")
            query, subject, pident, length, qlen, _qstart, _qend, slen, evalue, bitscore = values
            pident_number = number(pident)
            length_number = number(length)
            qlen_number = number(qlen)
            subject_length = number(slen)
            evalue_number = number(evalue)
            bitscore_number = number(bitscore)
            query_coverage = 3 * length_number / qlen_number if qlen_number else math.nan
            if not (
                pident_number >= 35
                and query_coverage >= 0.30
                and evalue_number <= 1e-10
                and bitscore_number >= 80
            ):
                continue
            candidates.setdefault(query, []).append(
                {
                    "query": query,
                    "subject": subject,
                    "pident": format_value(pident_number),
                    "length": format_value(length_number),
                    "query_coverage": format_value(query_coverage),
                    "subject_coverage": format_value(length_number / subject_length) if subject_length else "",
                    "evalue": evalue,
                    "bitscore": format_value(bitscore_number),
                }
            )
    return {
        query: sorted(rows, key=lambda row: (-number(row["bitscore"]), -number(row["pident"]), number(row["evalue"]), row["subject"]))[0]
        for query, rows in candidates.items()
    }


def build_expression_bridge(
    blast_path: Path,
    bridge_path: Path,
    tx2gene_path: Path,
    two_week_path: Path,
    four_week_path: Path,
    output_path: Path,
) -> list[dict[str, str]]:
    """Write best supported P. strobus-to-Route 1 homology bridges with DE results."""
    bridge_by_protein = {row["ptae_protein_id"]: row for row in read_tsv(bridge_path)}
    transcript_to_gene = {row["transcript_id"]: row["gene_id"] for row in read_tsv(tx2gene_path)}
    two_week = read_deseq2(two_week_path)
    four_week = read_deseq2(four_week_path)
    rows: list[dict[str, str]] = []
    for transcript_id, hit in sorted(best_hits(blast_path).items()):
        bridge = bridge_by_protein.get(hit["subject"])
        gene_id = transcript_to_gene.get(transcript_id)
        if bridge is None or gene_id is None:
            continue
        two = two_week.get(gene_id)
        four = four_week.get(gene_id)
        rows.append(
            {
                "pstrobus_transcript_id": transcript_id,
                "pstrobus_gene_id": gene_id,
                "ptae_protein_id": hit["subject"],
                "orthogroup_id": bridge["orthogroup_id"],
                "effector_target_network_role": bridge["effector_target_network_role"],
                "host_module": bridge["host_module"],
                "host_network_role": bridge["host_network_role"],
                "pident": hit["pident"],
                "alignment_length_aa": hit["length"],
                "query_translation_coverage": hit["query_coverage"],
                "subject_coverage_descriptive": hit["subject_coverage"],
                "evalue": hit["evalue"],
                "bitscore": hit["bitscore"],
                "two_week_log2_fold_change": two.get("log2FoldChange", "") if two else "",
                "two_week_padj": two.get("padj", "") if two else "",
                "two_week_fdr_status": fdr_status(two),
                "four_week_log2_fold_change": four.get("log2FoldChange", "") if four else "",
                "four_week_padj": four.get("padj", "") if four else "",
                "four_week_fdr_status": fdr_status(four),
                "bridge_evidence_type": "P_strobus_transcript_to_P_taeda_effector_target_network_orthogroup_homology_proxy",
                "claim_ceiling": CLAIM_CEILING,
            }
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Wrote %d threshold-passing P. strobus Route 1 homology bridges", len(rows))
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blast", type=Path, required=True)
    parser.add_argument("--bridge", type=Path, required=True)
    parser.add_argument("--tx2gene", type=Path, required=True)
    parser.add_argument("--two-week", type=Path, required=True)
    parser.add_argument("--four-week", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = parse_args()
    build_expression_bridge(args.blast, args.bridge, args.tx2gene, args.two_week, args.four_week, args.output)


if __name__ == "__main__":
    main()
