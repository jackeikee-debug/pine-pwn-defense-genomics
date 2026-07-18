#!/usr/bin/env python3
"""Summarize conservative manual-review decisions for candidate families."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


FIELDS = [
    "decision_rank",
    "orthogroup_id",
    "evidence_tier",
    "regional_direction",
    "module_ids",
    "matched_keywords",
    "driver_risk",
    "focal_region_cluster_status",
    "focal_region_mrca_purity",
    "tree_review_priority",
    "preliminary_review_decision",
    "manuscript_use",
    "decision_rationale",
    "alignment_path",
    "tree_path",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_candidate_review_decisions(
    checklist_path: Path,
    metrics_path: Path,
    output_path: Path,
    markdown_path: Path,
) -> list[dict[str, str]]:
    metrics_by_og = {row["orthogroup_id"]: row for row in read_tsv(metrics_path)}
    rows: list[dict[str, str]] = []
    for checklist in read_tsv(checklist_path):
        metric = metrics_by_og.get(checklist["orthogroup_id"], {})
        decision, manuscript_use, rationale = classify_decision(checklist, metric)
        rows.append(
            {
                "decision_rank": checklist.get("check_rank", ""),
                "orthogroup_id": checklist.get("orthogroup_id", ""),
                "evidence_tier": checklist.get("evidence_tier", ""),
                "regional_direction": checklist.get("regional_direction", ""),
                "module_ids": checklist.get("module_ids", ""),
                "matched_keywords": checklist.get("matched_keywords", ""),
                "driver_risk": checklist.get("driver_risk", ""),
                "focal_region_cluster_status": metric.get("focal_region_cluster_status", ""),
                "focal_region_mrca_purity": metric.get("focal_region_mrca_purity", ""),
                "tree_review_priority": metric.get("manual_tree_review_priority", ""),
                "preliminary_review_decision": decision,
                "manuscript_use": manuscript_use,
                "decision_rationale": rationale,
                "alignment_path": checklist.get("alignment_path", ""),
                "tree_path": checklist.get("tree_path", ""),
            }
        )
    write_tsv(output_path, FIELDS, rows)
    write_markdown(markdown_path, rows)
    return rows


def classify_decision(checklist: dict[str, str], metric: dict[str, str]) -> tuple[str, str, str]:
    tier = checklist.get("evidence_tier", "")
    driver_risk = checklist.get("driver_risk", "")
    tree_status = metric.get("focal_region_cluster_status", "")

    if tier == "caution_candidate" or driver_risk == "high":
        return (
            "deprioritize_species_or_large_family_driver",
            "background_or_supplement_only",
            "CAFE5 support overlaps high driver risk or caution-tier family flags.",
        )
    if tree_status == "region_biased_treewide_distribution":
        return (
            "priority_manual_review_no_regional_claim",
            "candidate_for_manual_tree_review",
            "Model support and regional copy bias are present, but the focal MRCA spans the tree.",
        )
    if tier == "moderate_candidate" and driver_risk == "low":
        return (
            "retain_functional_candidate_tree_mixed",
            "candidate_for_functional_context",
            "Driver risk is low, but tree topology is mixed rather than a regional clade.",
        )
    return (
        "retain_as_unresolved_manual_review",
        "candidate_for_supplementary_review",
        "Evidence is model-supported but topology and/or driver risk remain unresolved.",
    )


def write_markdown(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Candidate Review Decisions",
        "",
        "These decisions are conservative manual-review labels, not final biological conclusions.",
        "",
    ]
    for row in rows:
        lines.extend(
            [
                f"## {row['decision_rank']}. {row['orthogroup_id']}",
                "",
                f"- Decision: {row['preliminary_review_decision']}",
                f"- Manuscript use: {row['manuscript_use']}",
                f"- Evidence tier: {row['evidence_tier']}",
                f"- Direction: {row['regional_direction']}",
                f"- Module/keyword: {row['module_ids']} / {row['matched_keywords']}",
                f"- Driver risk: {row['driver_risk']}",
                f"- Tree status: {row['focal_region_cluster_status']} (purity {row['focal_region_mrca_purity']})",
                f"- Rationale: {row['decision_rationale']}",
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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checklist", type=Path, required=True)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = build_candidate_review_decisions(
        checklist_path=args.checklist,
        metrics_path=args.metrics,
        output_path=args.output,
        markdown_path=args.markdown,
    )
    print(f"Candidate review decision rows written: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
