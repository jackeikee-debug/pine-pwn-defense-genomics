#!/usr/bin/env python3
"""Overlay P. massoniana Route 1 members with the 1-dpi DESeq2 interaction test."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


INTERACTION_FIELDS = [
    "interaction_base_mean",
    "interaction_log2_fold_change",
    "interaction_lfc_se",
    "interaction_pvalue",
    "interaction_padj",
    "interaction_fdr_status",
    "response_difference_direction",
    "interaction_coefficient",
    "interaction_claim_ceiling",
]
CLAIM = (
    "A genotype-by-inoculum expression interaction at 1 dpi is response-divergence evidence only; "
    "it does not establish direct effector targeting or causal resistance/tolerance."
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def numeric(value: str | None) -> float | None:
    if value in (None, "", "NA", "NaN"):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def interaction_status(row: dict[str, str] | None) -> str:
    if row is None:
        return "not_tested_or_filtered"
    padj = numeric(row.get("interaction_padj"))
    fold_change = numeric(row.get("interaction_log2_fold_change"))
    if padj is None or fold_change is None:
        return "not_tested_or_filtered"
    if padj < 0.05 and abs(fold_change) >= 1:
        return "fdr_lt_0.05_abs_lfc_ge_1"
    if padj < 0.05:
        return "fdr_lt_0.05_abs_lfc_lt_1"
    return "fdr_ge_0.05"


def response_direction(row: dict[str, str] | None) -> str:
    if row is None:
        return "not_tested_or_filtered"
    fold_change = numeric(row.get("interaction_log2_fold_change"))
    if fold_change is None:
        return "not_tested_or_filtered"
    if fold_change > 0:
        return "stronger_or_less_repressed_in_susceptible"
    if fold_change < 0:
        return "weaker_or_more_repressed_in_susceptible"
    return "no_estimated_response_difference"


def build_overlay(members_path: Path, interaction_path: Path, output_path: Path) -> list[dict[str, str]]:
    members = read_tsv(members_path)
    interactions = {row["feature_id"]: row for row in read_tsv(interaction_path)}
    if not members:
        raise ValueError("Route 1 member table is empty")

    rows: list[dict[str, str]] = []
    for member in members:
        interaction = interactions.get(member["pmas_feature_id"])
        row = dict(member)
        row.update(
            {
                "interaction_base_mean": interaction.get("baseMean", "") if interaction else "",
                "interaction_log2_fold_change": interaction.get("interaction_log2_fold_change", "") if interaction else "",
                "interaction_lfc_se": interaction.get("interaction_lfc_se", "") if interaction else "",
                "interaction_pvalue": interaction.get("interaction_pvalue", "") if interaction else "",
                "interaction_padj": interaction.get("interaction_padj", "") if interaction else "",
                "interaction_fdr_status": interaction_status(interaction),
                "response_difference_direction": response_direction(interaction),
                "interaction_coefficient": interaction.get("interaction_coefficient", "") if interaction else "",
                "interaction_claim_ceiling": CLAIM,
            }
        )
        rows.append(row)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(members[0]) + INTERACTION_FIELDS
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--members", type=Path, required=True)
    parser.add_argument("--interaction", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = build_overlay(args.members, args.interaction, args.output)
    print(f"P. massoniana Route 1 interaction overlay: {len(rows)} exact members")


if __name__ == "__main__":
    main()
