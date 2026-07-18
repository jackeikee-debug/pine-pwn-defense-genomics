#!/usr/bin/env python3
"""Rank regional candidate orthogroups by combined evidence layers."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


FIELDS = [
    "evidence_rank",
    "orthogroup_id",
    "evidence_tier",
    "evidence_score",
    "regional_copy_directions",
    "module_ids",
    "effector_classes",
    "matched_keywords",
    "best_abs_log2_ratio",
    "best_log2_ratio",
    "review_tier",
    "review_flags",
    "cafe_family_pvalue",
    "cafe_family_significant_0_05",
    "significant_branch_count",
    "max_abs_branch_change",
    "significant_expansion_branches",
    "significant_contraction_branches",
    "all_significant_branches",
    "east_asia_species_counts",
    "north_america_species_counts",
    "outgroup_species_counts",
    "alignment_path",
    "tree_path",
    "interpretation_note",
]


TIER_ORDER = {
    "strong_candidate": 0,
    "moderate_candidate": 1,
    "caution_candidate": 2,
    "weak_or_unresolved": 3,
}


CAUTION_REVIEW_TIERS = {
    "caution_single_species_driver",
    "caution_large_single_species_family",
    "caution_large_family",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def classify_candidate(
    review_tier: str,
    cafe_significant: str,
    significant_branch_count: int,
    max_abs_branch_change: int,
) -> str:
    has_model_support = cafe_significant == "yes" and significant_branch_count > 0 and max_abs_branch_change > 0
    if not has_model_support:
        return "weak_or_unresolved"
    if review_tier == "priority_low_risk":
        return "strong_candidate"
    if review_tier in CAUTION_REVIEW_TIERS:
        return "caution_candidate"
    return "moderate_candidate"


def evidence_score(
    tier: str,
    cafe_significant: str,
    significant_branch_count: int,
    max_abs_branch_change: int,
    best_abs_log2_ratio: float,
) -> int:
    score = 0
    if cafe_significant == "yes":
        score += 3
    if significant_branch_count > 0:
        score += 2
    if max_abs_branch_change >= 3:
        score += 1
    if best_abs_log2_ratio >= 1.5:
        score += 1
    if tier == "strong_candidate":
        score += 1
    if tier == "caution_candidate":
        score -= 1
    return score


def rank_candidate_evidence(
    regional_top_path: Path,
    validation_summary_path: Path,
    cafe_intersection_path: Path,
    output_path: Path,
    strong_output_path: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    regional_by_og = {row["orthogroup_id"]: row for row in read_tsv(regional_top_path)}
    validation_rows = read_tsv(validation_summary_path)
    cafe_by_og = {row["orthogroup_id"]: row for row in read_tsv(cafe_intersection_path)}

    rows: list[dict[str, str]] = []
    for validation in validation_rows:
        orthogroup_id = validation["orthogroup_id"]
        regional = regional_by_og.get(orthogroup_id, {})
        cafe = cafe_by_og.get(orthogroup_id, {})
        branch_count = to_int(cafe.get("significant_branch_count", "0"))
        max_change = to_int(cafe.get("max_abs_branch_change", "0"))
        cafe_significant = cafe.get("cafe_family_significant_0_05", "not_tested")
        best_ratio = to_float(validation.get("best_abs_log2_ratio", regional.get("best_abs_log2_ratio", "0")))
        tier = classify_candidate(
            review_tier=validation.get("review_tier", ""),
            cafe_significant=cafe_significant,
            significant_branch_count=branch_count,
            max_abs_branch_change=max_change,
        )
        score = evidence_score(
            tier=tier,
            cafe_significant=cafe_significant,
            significant_branch_count=branch_count,
            max_abs_branch_change=max_change,
            best_abs_log2_ratio=best_ratio,
        )
        rows.append(
            {
                "evidence_rank": "0",
                "orthogroup_id": orthogroup_id,
                "evidence_tier": tier,
                "evidence_score": str(score),
                "regional_copy_directions": validation.get("regional_copy_directions", regional.get("regional_copy_directions", "")),
                "module_ids": validation.get("module_ids", regional.get("module_ids", "")),
                "effector_classes": regional.get("effector_classes", ""),
                "matched_keywords": validation.get("matched_keywords", regional.get("matched_keywords", "")),
                "best_abs_log2_ratio": validation.get("best_abs_log2_ratio", regional.get("best_abs_log2_ratio", "")),
                "best_log2_ratio": validation.get("best_log2_ratio", regional.get("best_log2_ratio", "")),
                "review_tier": validation.get("review_tier", ""),
                "review_flags": validation.get("review_flags", ""),
                "cafe_family_pvalue": cafe.get("cafe_family_pvalue", ""),
                "cafe_family_significant_0_05": cafe_significant,
                "significant_branch_count": str(branch_count),
                "max_abs_branch_change": str(max_change),
                "significant_expansion_branches": cafe.get("significant_expansion_branches", ""),
                "significant_contraction_branches": cafe.get("significant_contraction_branches", ""),
                "all_significant_branches": cafe.get("all_significant_branches", ""),
                "east_asia_species_counts": regional.get("east_asia_species_counts", ""),
                "north_america_species_counts": regional.get("north_america_species_counts", ""),
                "outgroup_species_counts": regional.get("outgroup_species_counts", ""),
                "alignment_path": validation.get("alignment_path", ""),
                "tree_path": validation.get("tree_path", ""),
                "interpretation_note": build_interpretation_note(tier),
            }
        )

    rows.sort(key=rank_sort_key)
    for rank, row in enumerate(rows, start=1):
        row["evidence_rank"] = str(rank)
    strong_rows = [row for row in rows if row["evidence_tier"] == "strong_candidate"]
    write_tsv(output_path, FIELDS, rows)
    write_tsv(strong_output_path, FIELDS, strong_rows)
    return rows, strong_rows


def rank_sort_key(row: dict[str, str]) -> tuple[int, int, float, int, str]:
    return (
        TIER_ORDER.get(row["evidence_tier"], 99),
        -to_int(row["evidence_score"]),
        -to_float(row["best_abs_log2_ratio"]),
        to_int(row.get("evidence_rank", "0")),
        row["orthogroup_id"],
    )


def build_interpretation_note(tier: str) -> str:
    if tier == "strong_candidate":
        return "Regional copy-number contrast has CAFE5 family and branch support with low review-risk flags."
    if tier == "moderate_candidate":
        return "Regional copy-number contrast has CAFE5 support but requires context-aware manual review."
    if tier == "caution_candidate":
        return "CAFE5 support overlaps a family flagged for single-species or large-family effects; inspect manually before interpretation."
    return "Regional candidate lacks CAFE5 model support or was not tested in the current pilot run."


def to_int(value: str) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def to_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regional-top", type=Path, required=True)
    parser.add_argument("--validation-summary", type=Path, required=True)
    parser.add_argument("--cafe-intersections", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--strong-output", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows, strong_rows = rank_candidate_evidence(
        regional_top_path=args.regional_top,
        validation_summary_path=args.validation_summary,
        cafe_intersection_path=args.cafe_intersections,
        output_path=args.output,
        strong_output_path=args.strong_output,
    )
    print(f"Candidate evidence ranking rows written: {len(rows)}")
    print(f"Strong candidate rows written: {len(strong_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
