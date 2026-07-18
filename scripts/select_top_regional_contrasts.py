#!/usr/bin/env python3
"""Select top regional contrast candidates for downstream family and sequence tests."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


ASSOCIATION_PREFIX_FIELDS = [
    "regional_contrast_rank",
    "regional_group_rank",
    "selection_group",
    "abs_log2_ratio",
    "selection_reason",
]

ORTHOGROUP_FIELDS = [
    "orthogroup_rank",
    "orthogroup_id",
    "regional_copy_directions",
    "module_ids",
    "effector_classes",
    "best_abs_log2_ratio",
    "best_log2_ratio",
    "best_selection_group",
    "association_count",
    "matched_keywords",
    "east_asia_species_counts",
    "north_america_species_counts",
    "outgroup_species_counts",
    "max_host_priority_score",
    "best_mean_bitscore",
    "candidate_interpretation",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def select_top_regional_contrasts(
    contrast_path: Path,
    association_output: Path,
    orthogroup_output: Path,
    max_per_direction_module_effector: int = 5,
    include_balanced: bool = False,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    rows = [
        row
        for row in read_tsv(contrast_path)
        if include_balanced or row.get("regional_copy_direction") != "balanced"
    ]

    selected = select_association_rows(rows, max_per_direction_module_effector)
    orthogroups = summarize_unique_orthogroups(selected)

    association_fields = ASSOCIATION_PREFIX_FIELDS + list(rows[0].keys()) if rows else ASSOCIATION_PREFIX_FIELDS
    write_tsv(association_output, association_fields, selected)
    write_tsv(orthogroup_output, ORTHOGROUP_FIELDS, orthogroups)
    return selected, orthogroups


def select_association_rows(
    rows: list[dict[str, str]],
    max_per_direction_module_effector: int,
) -> list[dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        group = selection_group(row)
        grouped[group].append(row)

    selected: list[dict[str, str]] = []
    for group in sorted(grouped):
        group_rows = sorted(grouped[group], key=association_sort_key)
        cap = max_per_direction_module_effector if max_per_direction_module_effector > 0 else len(group_rows)
        for group_rank, row in enumerate(group_rows[:cap], start=1):
            enriched = {
                **row,
                "regional_contrast_rank": "0",
                "regional_group_rank": str(group_rank),
                "selection_group": group,
                "abs_log2_ratio": format_float(abs(to_float(row.get("log2_east_asia_vs_north_america_mean_copy_ratio", "")))),
                "selection_reason": build_selection_reason(row),
            }
            selected.append(enriched)

    selected.sort(key=association_sort_key)
    for rank, row in enumerate(selected, start=1):
        row["regional_contrast_rank"] = str(rank)
    return selected


def summarize_unique_orthogroups(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["orthogroup_id"]].append(row)

    summaries: list[dict[str, str]] = []
    for orthogroup_id, group_rows in grouped.items():
        best = sorted(group_rows, key=association_sort_key)[0]
        summaries.append(
            {
                "orthogroup_rank": "0",
                "orthogroup_id": orthogroup_id,
                "regional_copy_directions": join_values(row["regional_copy_direction"] for row in group_rows),
                "module_ids": join_values(row["module_id"] for row in group_rows),
                "effector_classes": join_values(row["effector_class"] for row in group_rows),
                "best_abs_log2_ratio": format_float(
                    abs(to_float(best.get("log2_east_asia_vs_north_america_mean_copy_ratio", "")))
                ),
                "best_log2_ratio": best.get("log2_east_asia_vs_north_america_mean_copy_ratio", ""),
                "best_selection_group": best.get("selection_group", selection_group(best)),
                "association_count": str(len(group_rows)),
                "matched_keywords": join_split_values(row.get("matched_keywords", "") for row in group_rows),
                "east_asia_species_counts": best.get("east_asia_species_counts", ""),
                "north_america_species_counts": best.get("north_america_species_counts", ""),
                "outgroup_species_counts": best.get("outgroup_species_counts", ""),
                "max_host_priority_score": format_float(max(to_float(row.get("host_priority_score", "")) for row in group_rows)),
                "best_mean_bitscore": best.get("mean_bitscore", ""),
                "candidate_interpretation": best.get("interpretation", ""),
            }
        )

    summaries.sort(
        key=lambda row: (
            -to_float(row["best_abs_log2_ratio"]),
            -to_float(row["max_host_priority_score"]),
            -to_float(row["best_mean_bitscore"]),
            row["orthogroup_id"],
        )
    )
    for rank, row in enumerate(summaries, start=1):
        row["orthogroup_rank"] = str(rank)
    return summaries


def selection_group(row: dict[str, str]) -> str:
    return f"{row['regional_copy_direction']}|{row['module_id']}|{row['effector_class']}"


def association_sort_key(row: dict[str, str]) -> tuple[float, float, float, int, str]:
    return (
        -abs(to_float(row.get("log2_east_asia_vs_north_america_mean_copy_ratio", ""))),
        -to_float(row.get("host_priority_score", "")),
        -to_float(row.get("mean_bitscore", "")),
        to_int(row.get("host_shortlist_rank", "")),
        row.get("orthogroup_id", ""),
    )


def build_selection_reason(row: dict[str, str]) -> str:
    return (
        f"direction={row.get('regional_copy_direction', '')};"
        f"abs_log2_ratio={abs(to_float(row.get('log2_east_asia_vs_north_america_mean_copy_ratio', ''))):.2f};"
        f"group={selection_group(row)}"
    )


def split_values(value: str) -> list[str]:
    if not value.strip():
        return []
    for separator in [";", ","]:
        if separator in value:
            return [part.strip() for part in value.split(separator) if part.strip()]
    return [value.strip()]


def join_split_values(values: object) -> str:
    parts: set[str] = set()
    for value in values:
        parts.update(split_values(str(value)))
    return ";".join(sorted(parts))


def join_values(values: object) -> str:
    return ";".join(sorted({str(value) for value in values if str(value)}))


def to_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def to_int(value: str) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 10**9


def format_float(value: float) -> str:
    return f"{value:.2f}"


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contrast", type=Path, required=True)
    parser.add_argument("--association-output", type=Path, required=True)
    parser.add_argument("--orthogroup-output", type=Path, required=True)
    parser.add_argument("--max-per-direction-module-effector", type=int, default=5)
    parser.add_argument("--include-balanced", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    association_rows, orthogroup_rows = select_top_regional_contrasts(
        contrast_path=args.contrast,
        association_output=args.association_output,
        orthogroup_output=args.orthogroup_output,
        max_per_direction_module_effector=args.max_per_direction_module_effector,
        include_balanced=args.include_balanced,
    )
    print(f"Top regional contrast associations written: {len(association_rows)}")
    print(f"Top regional contrast orthogroups written: {len(orthogroup_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
