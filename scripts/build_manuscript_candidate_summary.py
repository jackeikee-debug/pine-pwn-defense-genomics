#!/usr/bin/env python3
"""Build a manuscript-facing summary for reviewed candidate families."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


FIELDS = [
    "summary_rank",
    "orthogroup_id",
    "manuscript_section",
    "claim_ceiling",
    "regional_direction",
    "module_ids",
    "effector_classes",
    "matched_keywords",
    "evidence_tier",
    "manuscript_use",
    "preliminary_review_decision",
    "driver_risk",
    "tree_status",
    "cafe_family_pvalue",
    "significant_branch_count",
    "east_asia_species_counts",
    "north_america_species_counts",
    "outgroup_species_counts",
    "safe_wording",
    "avoid_wording",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_manuscript_candidate_summary(
    review_path: Path,
    decisions_path: Path,
    output_path: Path,
    markdown_path: Path,
) -> list[dict[str, str]]:
    decisions_by_og = {row["orthogroup_id"]: row for row in read_tsv(decisions_path)}
    rows: list[dict[str, str]] = []
    for review in read_tsv(review_path):
        decision = decisions_by_og.get(review["orthogroup_id"], {})
        manuscript_section, claim_ceiling = classify_manuscript_use(decision)
        rows.append(
            {
                "summary_rank": review.get("review_rank", ""),
                "orthogroup_id": review.get("orthogroup_id", ""),
                "manuscript_section": manuscript_section,
                "claim_ceiling": claim_ceiling,
                "regional_direction": review.get("regional_copy_directions", ""),
                "module_ids": review.get("module_ids", ""),
                "effector_classes": review.get("effector_classes", ""),
                "matched_keywords": review.get("matched_keywords", ""),
                "evidence_tier": review.get("evidence_tier", ""),
                "manuscript_use": decision.get("manuscript_use", ""),
                "preliminary_review_decision": decision.get("preliminary_review_decision", ""),
                "driver_risk": decision.get("driver_risk", ""),
                "tree_status": decision.get("focal_region_cluster_status", ""),
                "cafe_family_pvalue": review.get("cafe_family_pvalue", ""),
                "significant_branch_count": review.get("significant_branch_count", ""),
                "east_asia_species_counts": review.get("east_asia_species_counts", ""),
                "north_america_species_counts": review.get("north_america_species_counts", ""),
                "outgroup_species_counts": review.get("outgroup_species_counts", ""),
                "safe_wording": safe_wording(review, decision),
                "avoid_wording": avoid_wording(decision),
            }
        )
    write_tsv(output_path, FIELDS, rows)
    write_markdown(markdown_path, rows)
    return rows


def classify_manuscript_use(decision: dict[str, str]) -> tuple[str, str]:
    use = decision.get("manuscript_use", "")
    if use == "candidate_for_manual_tree_review":
        return (
            "main_text_priority_review",
            "model_supported_candidate_for_manual_review_not_regional_expansion",
        )
    if use == "candidate_for_functional_context":
        return (
            "main_text_functional_context",
            "functional_context_candidate_not_topology_supported_regional_family",
        )
    if use == "candidate_for_supplementary_review":
        return (
            "supplement_unresolved_candidate",
            "model_supported_unresolved_family_for_supplement_only",
        )
    return (
        "supplement_background",
        "background_or_supplement_only_no_main_claim",
    )


def safe_wording(review: dict[str, str], decision: dict[str, str]) -> str:
    return (
        f"{review.get('orthogroup_id', '')} is an effector-guided, CAFE5-pilot-supported "
        f"{review.get('matched_keywords', '')} candidate in {review.get('module_ids', '')}, "
        f"but current tree review status is {decision.get('focal_region_cluster_status', '')}."
    )


def avoid_wording(decision: dict[str, str]) -> str:
    if decision.get("manuscript_use") == "candidate_for_manual_tree_review":
        return "Do not call this a confirmed regional expansion or resistance mechanism."
    if decision.get("manuscript_use") == "candidate_for_functional_context":
        return "Do not imply topology-supported East Asia vs North America divergence."
    return "Do not use as a primary mechanistic claim."


def write_markdown(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Manuscript Candidate Summary",
        "",
        "These candidates must not be described as validated resistance genes or direct host-nematode coevolution evidence.",
        "",
    ]
    for section in [
        "main_text_priority_review",
        "main_text_functional_context",
        "supplement_unresolved_candidate",
        "supplement_background",
    ]:
        section_rows = [row for row in rows if row["manuscript_section"] == section]
        if not section_rows:
            continue
        lines.extend([f"## {section}", ""])
        for row in section_rows:
            lines.extend(
                [
                    f"- {row['orthogroup_id']} ({row['matched_keywords']}): {row['safe_wording']}",
                    f"  Claim ceiling: {row['claim_ceiling']}.",
                    f"  Avoid: {row['avoid_wording']}",
                ]
            )
        lines.append("")
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--review", type=Path, required=True)
    parser.add_argument("--decisions", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = build_manuscript_candidate_summary(
        review_path=args.review,
        decisions_path=args.decisions,
        output_path=args.output,
        markdown_path=args.markdown,
    )
    print(f"Manuscript candidate summary rows written: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
