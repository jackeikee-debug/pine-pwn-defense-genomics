#!/usr/bin/env python3
"""Build a manual review checklist for model-supported candidate families."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


EAST_ASIA_SPECIES = {"pden", "pmas", "ptab"}
NORTH_AMERICA_SPECIES = {"plam", "plon", "ptae"}

FIELDS = [
    "check_rank",
    "orthogroup_id",
    "evidence_tier",
    "regional_direction",
    "module_ids",
    "effector_classes",
    "matched_keywords",
    "nematode_role_in_evidence",
    "focal_region",
    "focal_species_counts",
    "focal_total_count",
    "focal_max_species",
    "focal_max_count",
    "focal_max_share",
    "driver_risk",
    "cafe_family_pvalue",
    "significant_branch_count",
    "significant_expansion_branches",
    "significant_contraction_branches",
    "review_flags",
    "review_question",
    "promote_before_manual_tree_review",
    "alignment_path",
    "tree_path",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_candidate_manual_checklist(review_path: Path, output_path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for source in read_tsv(review_path):
        direction = source.get("regional_copy_directions", "")
        focal_region, focal_counts = focal_region_counts(source, direction)
        total = sum(focal_counts.values())
        max_species, max_count = max_count_item(focal_counts)
        max_share = max_count / total if total else 0.0
        driver_risk = classify_driver_risk(max_share, source.get("review_flags", ""))
        rows.append(
            {
                "check_rank": source.get("review_rank", ""),
                "orthogroup_id": source.get("orthogroup_id", ""),
                "evidence_tier": source.get("evidence_tier", ""),
                "regional_direction": direction,
                "module_ids": source.get("module_ids", ""),
                "effector_classes": source.get("effector_classes", ""),
                "matched_keywords": source.get("matched_keywords", ""),
                "nematode_role_in_evidence": "indirect_effector_guided_anchor_not_cafe_covariate",
                "focal_region": focal_region,
                "focal_species_counts": format_counts(focal_counts),
                "focal_total_count": str(total),
                "focal_max_species": max_species,
                "focal_max_count": str(max_count),
                "focal_max_share": f"{max_share:.3f}",
                "driver_risk": driver_risk,
                "cafe_family_pvalue": source.get("cafe_family_pvalue", ""),
                "significant_branch_count": source.get("significant_branch_count", ""),
                "significant_expansion_branches": source.get("significant_expansion_branches", ""),
                "significant_contraction_branches": source.get("significant_contraction_branches", ""),
                "review_flags": source.get("review_flags", ""),
                "review_question": build_review_question(source, driver_risk, max_species, max_share),
                "promote_before_manual_tree_review": "no",
                "alignment_path": source.get("alignment_path", ""),
                "tree_path": source.get("tree_path", ""),
            }
        )
    write_tsv(output_path, FIELDS, rows)
    return rows


def focal_region_counts(source: dict[str, str], direction: str) -> tuple[str, dict[str, int]]:
    if "North_America_enriched" in direction:
        return "North_America", parse_counts(source.get("north_america_species_counts", ""))
    if "East_Asia_enriched" in direction:
        return "East_Asia", parse_counts(source.get("east_asia_species_counts", ""))
    return "mixed_or_balanced", {}


def parse_counts(value: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for item in value.split(";"):
        if not item or ":" not in item:
            continue
        key, raw_count = item.split(":", 1)
        try:
            counts[key] = int(float(raw_count))
        except ValueError:
            counts[key] = 0
    return counts


def format_counts(counts: dict[str, int]) -> str:
    return ";".join(f"{species}:{counts[species]}" for species in sorted(counts))


def max_count_item(counts: dict[str, int]) -> tuple[str, int]:
    if not counts:
        return "", 0
    return max(counts.items(), key=lambda item: (item[1], item[0]))


def classify_driver_risk(max_share: float, review_flags: str) -> str:
    if "single_species_driver" in review_flags or "large_family" in review_flags or max_share >= 0.70:
        return "high"
    if max_share >= 0.50:
        return "moderate"
    return "low"


def build_review_question(source: dict[str, str], driver_risk: str, max_species: str, max_share: float) -> str:
    if driver_risk == "high":
        return (
            f"Does the family tree support a regional signal, or is the contrast mainly driven by "
            f"single species {max_species} ({max_share:.1%} of focal-region copies)?"
        )
    if "outgroup_present" in source.get("review_flags", ""):
        return "Do outgroup sequences indicate an older conifer-family pattern rather than a Pinus regional change?"
    return "Do alignment and tree topology support treating this as a coherent regional candidate family?"


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--review", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = build_candidate_manual_checklist(args.review, args.output)
    print(f"Candidate manual checklist rows written: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
