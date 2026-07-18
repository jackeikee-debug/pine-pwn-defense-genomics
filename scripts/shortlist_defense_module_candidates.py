#!/usr/bin/env python3
"""Create a strict directional shortlist from defense-module candidate priorities."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


DIRECTIONAL_BIASES = {
    "susceptible_enriched",
    "tolerant_enriched",
    "susceptible_only",
    "tolerant_only",
}
SHORTLIST_FIELDS = [
    "shortlist_rank",
    "shortlist_group",
    "shortlist_group_rank",
    "shortlist_reason",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def shortlist_candidates(
    priority_path: Path,
    output_path: Path,
    max_per_module_direction: int,
    min_mean_bitscore: float,
    min_core_retention: float,
    min_core_jaccard: float,
    allowed_tiers: list[str],
) -> list[dict[str, str]]:
    rows = [
        row
        for row in read_tsv(priority_path)
        if passes_filters(
            row=row,
            min_mean_bitscore=min_mean_bitscore,
            min_core_retention=min_core_retention,
            min_core_jaccard=min_core_jaccard,
            allowed_tiers=allowed_tiers,
        )
    ]
    rows.sort(
        key=lambda row: (
            row["module_id"],
            row["copy_bias_direction"],
            -int(row["priority_score"]),
            -float(row["mean_bitscore"]),
            int(row["priority_rank"]),
        )
    )

    group_counts: defaultdict[str, int] = defaultdict(int)
    selected: list[dict[str, str]] = []
    for row in rows:
        group = f"{row['module_id']}|{row['copy_bias_direction']}"
        if max_per_module_direction > 0 and group_counts[group] >= max_per_module_direction:
            continue
        group_counts[group] += 1
        selected.append(
            {
                **row,
                "shortlist_rank": "0",
                "shortlist_group": group,
                "shortlist_group_rank": str(group_counts[group]),
                "shortlist_reason": build_reason(
                    row,
                    min_mean_bitscore=min_mean_bitscore,
                    min_core_retention=min_core_retention,
                    min_core_jaccard=min_core_jaccard,
                ),
            }
        )

    selected.sort(key=lambda row: int(row["priority_rank"]))
    for index, row in enumerate(selected, start=1):
        row["shortlist_rank"] = str(index)

    fields = SHORTLIST_FIELDS + [field for field in rows[0].keys()] if rows else SHORTLIST_FIELDS
    write_tsv(output_path, fields, selected)
    return selected


def passes_filters(
    row: dict[str, str],
    min_mean_bitscore: float,
    min_core_retention: float,
    min_core_jaccard: float,
    allowed_tiers: list[str],
) -> bool:
    return (
        row["priority_tier"] in allowed_tiers
        and row["copy_bias_direction"] in DIRECTIONAL_BIASES
        and float(row["mean_bitscore"]) >= min_mean_bitscore
        and float(row["pilot_core_retention"]) >= min_core_retention
        and float(row["core_jaccard"]) >= min_core_jaccard
    )


def build_reason(
    row: dict[str, str],
    min_mean_bitscore: float,
    min_core_retention: float,
    min_core_jaccard: float,
) -> str:
    return (
        f"directional={row['copy_bias_direction']};"
        f"tier={row['priority_tier']};"
        f"mean_bitscore>={min_mean_bitscore:g};"
        f"pilot_core_retention>={min_core_retention:g};"
        f"core_jaccard>={min_core_jaccard:g}"
    )


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--priority", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-per-module-direction", type=int, default=20)
    parser.add_argument("--min-mean-bitscore", type=float, default=150.0)
    parser.add_argument("--min-core-retention", type=float, default=0.95)
    parser.add_argument("--min-core-jaccard", type=float, default=0.80)
    parser.add_argument("--allowed-tiers", nargs="+", default=["high"])
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = shortlist_candidates(
        priority_path=args.priority,
        output_path=args.output,
        max_per_module_direction=args.max_per_module_direction,
        min_mean_bitscore=args.min_mean_bitscore,
        min_core_retention=args.min_core_retention,
        min_core_jaccard=args.min_core_jaccard,
        allowed_tiers=args.allowed_tiers,
    )
    print(f"Defense module candidate shortlist rows written: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
