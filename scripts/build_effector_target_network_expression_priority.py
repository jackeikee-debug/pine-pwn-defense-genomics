#!/usr/bin/env python3
"""Integrate Route 1 candidates with three tximport-DESeq2 contrasts."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


DE_COLUMNS = ["transcript_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
OUTPUT_FIELDS = [
    "transcript_id",
    "orthogroup_ids",
    "effector_target_network_seed_symbols",
    "mapping_status",
    "best_effector_target_network_candidate_pd_gene",
    "primary_log2_fold_change",
    "primary_padj",
    "bxyl_vs_water_log2_fold_change",
    "bxyl_vs_water_padj",
    "bthai_vs_water_log2_fold_change",
    "bthai_vs_water_padj",
    "expression_evidence_tier",
    "claim_ceiling",
]
CLAIM_CEILING = (
    "Differential expression is infection-response support only; it is not direct effector-target, "
    "causal resistance, or cross-region expression evidence"
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


def is_significant(row: dict[str, str] | None, threshold: float = 0.05) -> bool:
    if row is None:
        return False
    try:
        return float(row["padj"]) < threshold
    except (KeyError, TypeError, ValueError):
        return False


def value(row: dict[str, str] | None, field: str) -> str:
    return row.get(field, "") if row else ""


def evidence_tier(primary: dict[str, str] | None, bxyl_water: dict[str, str] | None, bthai_water: dict[str, str] | None) -> str:
    if is_significant(primary) and is_significant(bxyl_water) and not is_significant(bthai_water):
        return "pathogen_associated_expression_candidate"
    if is_significant(primary):
        return "pathogen_vs_nonpathogen_expression_candidate"
    if is_significant(bxyl_water):
        return "bxyl_vs_water_expression_candidate"
    return "no_fdr_expression_support"


def build_priority_table(
    candidates_path: Path,
    primary_path: Path,
    bxyl_water_path: Path,
    bthai_water_path: Path,
    output_path: Path,
    report_path: Path,
) -> list[dict[str, str]]:
    primary = read_deseq2(primary_path)
    bxyl_water = read_deseq2(bxyl_water_path)
    bthai_water = read_deseq2(bthai_water_path)
    rows: list[dict[str, str]] = []
    for candidate in read_tsv(candidates_path):
        transcript_id = candidate["transcript_id"]
        primary_row = primary.get(transcript_id)
        bxyl_water_row = bxyl_water.get(transcript_id)
        bthai_water_row = bthai_water.get(transcript_id)
        rows.append(
            {
                "transcript_id": transcript_id,
                "orthogroup_ids": candidate.get("orthogroup_ids", ""),
                "effector_target_network_seed_symbols": candidate.get("effector_target_network_seed_symbols", ""),
                "mapping_status": candidate.get("mapping_status", ""),
                "best_effector_target_network_candidate_pd_gene": candidate.get("best_effector_target_network_candidate_pd_gene", ""),
                "primary_log2_fold_change": value(primary_row, "log2FoldChange"),
                "primary_padj": value(primary_row, "padj"),
                "bxyl_vs_water_log2_fold_change": value(bxyl_water_row, "log2FoldChange"),
                "bxyl_vs_water_padj": value(bxyl_water_row, "padj"),
                "bthai_vs_water_log2_fold_change": value(bthai_water_row, "log2FoldChange"),
                "bthai_vs_water_padj": value(bthai_water_row, "padj"),
                "expression_evidence_tier": evidence_tier(primary_row, bxyl_water_row, bthai_water_row),
                "claim_ceiling": CLAIM_CEILING,
            }
        )
    write_tsv(output_path, rows)
    write_report(report_path, rows, primary_path, bxyl_water_path, bthai_water_path)
    return rows


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_report(path: Path, rows: list[dict[str, str]], primary: Path, bxyl_water: Path, bthai_water: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tiers = Counter(row["expression_evidence_tier"] for row in rows)
    lines = [
        "# Route 1 Expression Priority",
        "",
        f"- Primary contrast: `{primary}` (B. xylophilus vs B. thailandae)",
        f"- Secondary contrast: `{bxyl_water}` (B. xylophilus vs water)",
        f"- Secondary contrast: `{bthai_water}` (B. thailandae vs water)",
        f"- Candidate Trinity genes: {len(rows)}",
        f"- Pathogen-associated expression candidates: {tiers['pathogen_associated_expression_candidate']}",
        f"- Pathogen-vs-nonpathogen candidates: {tiers['pathogen_vs_nonpathogen_expression_candidate']}",
        f"- Bxyl-vs-water-only candidates: {tiers['bxyl_vs_water_expression_candidate']}",
        f"- Candidates without FDR support: {tiers['no_fdr_expression_support']}",
        "",
        "Expression tiers are prioritization evidence, not direct effector-target validation or cross-region expression evidence.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--bxyl-water", type=Path, required=True)
    parser.add_argument("--bthai-water", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    build_priority_table(args.candidates, args.primary, args.bxyl_water, args.bthai_water, args.output, args.report)


if __name__ == "__main__":
    main()
