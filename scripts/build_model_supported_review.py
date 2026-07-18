#!/usr/bin/env python3
"""Create a focused review packet for model-supported candidate orthogroups."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


SUPPORTED_TIERS = {"moderate_candidate", "caution_candidate"}

FIELDS = [
    "review_rank",
    "orthogroup_id",
    "evidence_tier",
    "evidence_score",
    "regional_copy_directions",
    "module_ids",
    "effector_classes",
    "matched_keywords",
    "best_abs_log2_ratio",
    "review_tier",
    "review_flags",
    "cafe_family_pvalue",
    "significant_branch_count",
    "max_abs_branch_change",
    "significant_expansion_branches",
    "significant_contraction_branches",
    "east_asia_species_counts",
    "north_america_species_counts",
    "outgroup_species_counts",
    "alignment_path",
    "tree_path",
    "recommended_action",
    "interpretation_note",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_model_supported_review(
    ranking_path: Path,
    output_path: Path,
    markdown_path: Path,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for row in read_tsv(ranking_path):
        if row.get("evidence_tier") not in SUPPORTED_TIERS:
            continue
        rows.append(
            {
                "review_rank": "0",
                "orthogroup_id": row.get("orthogroup_id", ""),
                "evidence_tier": row.get("evidence_tier", ""),
                "evidence_score": row.get("evidence_score", ""),
                "regional_copy_directions": row.get("regional_copy_directions", ""),
                "module_ids": row.get("module_ids", ""),
                "effector_classes": row.get("effector_classes", ""),
                "matched_keywords": row.get("matched_keywords", ""),
                "best_abs_log2_ratio": row.get("best_abs_log2_ratio", ""),
                "review_tier": row.get("review_tier", ""),
                "review_flags": row.get("review_flags", ""),
                "cafe_family_pvalue": row.get("cafe_family_pvalue", ""),
                "significant_branch_count": row.get("significant_branch_count", ""),
                "max_abs_branch_change": row.get("max_abs_branch_change", ""),
                "significant_expansion_branches": row.get("significant_expansion_branches", ""),
                "significant_contraction_branches": row.get("significant_contraction_branches", ""),
                "east_asia_species_counts": row.get("east_asia_species_counts", ""),
                "north_america_species_counts": row.get("north_america_species_counts", ""),
                "outgroup_species_counts": row.get("outgroup_species_counts", ""),
                "alignment_path": row.get("alignment_path", ""),
                "tree_path": row.get("tree_path", ""),
                "recommended_action": recommended_action(row),
                "interpretation_note": row.get("interpretation_note", ""),
            }
        )

    rows.sort(key=review_sort_key)
    for index, row in enumerate(rows, start=1):
        row["review_rank"] = str(index)
    write_tsv(output_path, FIELDS, rows)
    write_markdown(markdown_path, rows)
    return rows


def recommended_action(row: dict[str, str]) -> str:
    tier = row.get("evidence_tier", "")
    flags = row.get("review_flags", "")
    if tier == "caution_candidate" or "single_species_driver" in flags or "large_family" in flags:
        return "inspect_single_species_or_large_family_driver"
    if "outgroup_present" in flags:
        return "inspect_outgroup_context_and_branch_support"
    return "inspect_branch_support_and_family_tree"


def review_sort_key(row: dict[str, str]) -> tuple[int, int, float, str]:
    tier_rank = 0 if row.get("evidence_tier") == "moderate_candidate" else 1
    return (
        tier_rank,
        -to_int(row.get("evidence_score", "")),
        -to_float(row.get("best_abs_log2_ratio", "")),
        row.get("orthogroup_id", ""),
    )


def write_markdown(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Model-Supported Candidate Review",
        "",
        "This review packet contains the regional candidate orthogroups with CAFE5 pilot support.",
        "These are prioritization targets, not final mechanistic claims.",
        "",
        f"Total candidates for manual review: {len(rows)}",
        "",
    ]
    for row in rows:
        lines.extend(
            [
                f"## {row['review_rank']}. {row['orthogroup_id']} ({row['evidence_tier']})",
                "",
                f"- Direction: {row['regional_copy_directions']}",
                f"- Modules: {row['module_ids']}",
                f"- Effector classes: {row['effector_classes']}",
                f"- Matched keywords: {row['matched_keywords']}",
                f"- Evidence score: {row['evidence_score']}",
                f"- CAFE5 family p-value: {row['cafe_family_pvalue']}",
                f"- Significant branches: {row['significant_branch_count']}",
                f"- Expansion branches: {row['significant_expansion_branches'] or 'none'}",
                f"- Contraction branches: {row['significant_contraction_branches'] or 'none'}",
                f"- Review flags: {row['review_flags']}",
                f"- Recommended action: {row['recommended_action']}",
                f"- Alignment: `{row['alignment_path']}`",
                f"- Tree: `{row['tree_path']}`",
                "",
            ]
        )
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ranking", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = build_model_supported_review(
        ranking_path=args.ranking,
        output_path=args.output,
        markdown_path=args.markdown,
    )
    print(f"Model-supported candidate review rows written: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
