#!/usr/bin/env python3
"""Parse CAFE5 Base-model outputs and intersect them with regional candidates."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


MODEL_FIELDS = [
    "model",
    "final_negative_log_likelihood",
    "lambda",
    "maximum_possible_lambda",
    "attempted_values",
    "rejected_percent",
]

CLADE_FIELDS = ["taxon_id", "taxon_label", "node_id", "increase", "decrease"]
FAMILY_FIELDS = ["family_id", "family_pvalue", "significant_0_05"]
BRANCH_FIELDS = ["family_id", "branch_id", "taxon_label", "node_id", "change", "direction", "branch_pvalue"]
TOP_INTERSECTION_FIELDS = [
    "orthogroup_id",
    "review_tier",
    "regional_copy_directions",
    "module_ids",
    "matched_keywords",
    "best_abs_log2_ratio",
    "cafe_family_pvalue",
    "cafe_family_significant_0_05",
    "significant_branch_count",
    "significant_expansion_branches",
    "significant_contraction_branches",
    "max_abs_branch_change",
    "all_significant_branches",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_cafe5_outputs(cafe_dir: Path) -> tuple[dict[str, str], list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    model = parse_model_summary(cafe_dir / "Base_results.txt")
    clades = parse_clade_results(cafe_dir / "Base_clade_results.txt")
    families = parse_family_results(cafe_dir / "Base_family_results.txt")
    branches = parse_branch_tables(cafe_dir / "Base_change.tab", cafe_dir / "Base_branch_probabilities.tab")
    return model, clades, families, branches


def parse_model_summary(path: Path) -> dict[str, str]:
    text = path.read_text(encoding="utf-8")
    model_match = re.search(r"Model\s+(\S+)\s+Final Likelihood \(-lnL\):\s+(\S+)", text)
    lambda_match = re.search(r"Lambda:\s+(\S+)", text)
    max_lambda_match = re.search(r"Maximum possible lambda for this topology:\s+(\S+)", text)
    attempted_match = re.search(r"(\d+) values were attempted \(([^%]+)% rejected\)", text)
    return {
        "model": model_match.group(1) if model_match else "",
        "final_negative_log_likelihood": model_match.group(2) if model_match else "",
        "lambda": lambda_match.group(1) if lambda_match else "",
        "maximum_possible_lambda": max_lambda_match.group(1) if max_lambda_match else "",
        "attempted_values": attempted_match.group(1) if attempted_match else "",
        "rejected_percent": attempted_match.group(2) if attempted_match else "",
    }


def parse_clade_results(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for row in read_hash_header_tsv(path):
        label, node = split_branch_id(row["Taxon_ID"])
        rows.append(
            {
                "taxon_id": row["Taxon_ID"],
                "taxon_label": label,
                "node_id": node,
                "increase": row["Increase"],
                "decrease": row["Decrease"],
            }
        )
    return rows


def parse_family_results(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for row in read_hash_header_tsv(path):
        rows.append(
            {
                "family_id": row["FamilyID"],
                "family_pvalue": row["pvalue"],
                "significant_0_05": "yes" if row["Significant at 0.05"].strip().lower() == "y" else "no",
            }
        )
    return rows


def parse_branch_tables(change_path: Path, probability_path: Path) -> list[dict[str, str]]:
    probability_by_family = {row["FamilyID"]: row for row in read_tsv(probability_path)}
    rows: list[dict[str, str]] = []
    for change_row in read_tsv(change_path):
        family_id = change_row["FamilyID"]
        probability_row = probability_by_family.get(family_id, {})
        for branch_id, change in change_row.items():
            if branch_id == "FamilyID":
                continue
            label, node = split_branch_id(branch_id)
            rows.append(
                {
                    "family_id": family_id,
                    "branch_id": branch_id,
                    "taxon_label": label,
                    "node_id": node,
                    "change": change,
                    "direction": classify_change(change),
                    "branch_pvalue": probability_row.get(branch_id, ""),
                }
            )
    return rows


def read_hash_header_tsv(path: Path) -> list[dict[str, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if lines and lines[0].startswith("#"):
        lines[0] = lines[0].lstrip("#")
    reader = csv.DictReader(lines, delimiter="\t")
    return list(reader)


def split_branch_id(branch_id: str) -> tuple[str, str]:
    match = re.match(r"^(?:(?P<label>[^<]+))?<(?P<node>\d+)>$", branch_id)
    if not match:
        return branch_id, ""
    return match.group("label") or "", match.group("node")


def classify_change(value: str) -> str:
    number = to_int(value)
    if number > 0:
        return "expansion"
    if number < 0:
        return "contraction"
    return "no_change"


def intersect_top_candidates(
    top_candidates_path: Path,
    family_results: list[dict[str, str]],
    branch_results: list[dict[str, str]],
    branch_pvalue_threshold: float = 0.05,
) -> list[dict[str, str]]:
    family_by_id = {row["family_id"]: row for row in family_results}
    branches_by_family: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in branch_results:
        if is_significant_branch(row.get("branch_pvalue", ""), branch_pvalue_threshold):
            branches_by_family[row["family_id"]].append(row)

    rows: list[dict[str, str]] = []
    for top in read_tsv(top_candidates_path):
        orthogroup_id = top["orthogroup_id"]
        family = family_by_id.get(orthogroup_id, {})
        significant_branches = sorted(
            branches_by_family.get(orthogroup_id, []),
            key=lambda row: (-abs(to_int(row["change"])), row["branch_id"]),
        )
        rows.append(
            {
                "orthogroup_id": orthogroup_id,
                "review_tier": top.get("review_tier", ""),
                "regional_copy_directions": top.get("regional_copy_directions", ""),
                "module_ids": top.get("module_ids", ""),
                "matched_keywords": top.get("matched_keywords", ""),
                "best_abs_log2_ratio": top.get("best_abs_log2_ratio", ""),
                "cafe_family_pvalue": family.get("family_pvalue", ""),
                "cafe_family_significant_0_05": family.get("significant_0_05", "not_tested"),
                "significant_branch_count": str(len(significant_branches)),
                "significant_expansion_branches": format_branch_list(significant_branches, "expansion"),
                "significant_contraction_branches": format_branch_list(significant_branches, "contraction"),
                "max_abs_branch_change": str(max((abs(to_int(row["change"])) for row in significant_branches), default=0)),
                "all_significant_branches": format_branch_list(significant_branches),
            }
        )
    return rows


def is_significant_branch(value: str, threshold: float) -> bool:
    if value == "N/A":
        return False
    try:
        return float(value) <= threshold
    except (TypeError, ValueError):
        return False


def format_branch_list(rows: list[dict[str, str]], direction: str | None = None) -> str:
    filtered = [row for row in rows if direction is None or row["direction"] == direction]
    return ";".join(f"{row['branch_id']}:{row['change']}" for row in filtered)


def write_one_row(path: Path, fields: list[str], row: dict[str, str]) -> None:
    write_tsv(path, fields, [row])


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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cafe-dir", type=Path, required=True)
    parser.add_argument("--top-candidates", type=Path, required=True)
    parser.add_argument("--model-summary-output", type=Path, required=True)
    parser.add_argument("--clade-output", type=Path, required=True)
    parser.add_argument("--family-output", type=Path, required=True)
    parser.add_argument("--branch-output", type=Path, required=True)
    parser.add_argument("--top-intersection-output", type=Path, required=True)
    parser.add_argument("--branch-pvalue-threshold", type=float, default=0.05)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    model, clades, families, branches = parse_cafe5_outputs(args.cafe_dir)
    top_rows = intersect_top_candidates(
        args.top_candidates,
        families,
        branches,
        branch_pvalue_threshold=args.branch_pvalue_threshold,
    )
    write_one_row(args.model_summary_output, MODEL_FIELDS, model)
    write_tsv(args.clade_output, CLADE_FIELDS, clades)
    write_tsv(args.family_output, FAMILY_FIELDS, families)
    write_tsv(args.branch_output, BRANCH_FIELDS, branches)
    write_tsv(args.top_intersection_output, TOP_INTERSECTION_FIELDS, top_rows)
    print(f"CAFE5 families parsed: {len(families)}")
    print(f"CAFE5 branch rows parsed: {len(branches)}")
    print(f"Top candidate CAFE intersections written: {len(top_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
