#!/usr/bin/env python3
"""Aggregate P. massoniana interaction statistics at orthogroup level."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


FIELDS = [
    "orthogroup_id",
    "member_feature_ids",
    "tested_member_count",
    "finite_pvalue_member_count",
    "significant_member_count",
    "significant_feature_ids",
    "min_member_padj",
    "max_abs_interaction_lfc",
    "direction_status",
    "simes_pvalue",
    "simes_padj",
    "module_ids",
    "mechanism_axes",
    "pmas_family_size",
    "mapping_completeness",
]


def _number(value: str) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _format(value: float | None) -> str:
    return "" if value is None else f"{value:.12g}"


def _split(value: str) -> set[str]:
    return {part.strip() for part in (value or "").split(";") if part.strip()}


def simes_pvalue(values: list[float]) -> float | None:
    """Return the Simes combined p-value for finite p-values."""
    finite = sorted(value for value in values if math.isfinite(value))
    if not finite:
        return None
    count = len(finite)
    return min(1.0, min(count * value / rank for rank, value in enumerate(finite, 1)))


def _bh_adjust(values: list[float | None]) -> list[float | None]:
    finite = [(index, value) for index, value in enumerate(values) if value is not None]
    adjusted: list[float | None] = [None] * len(values)
    running = 1.0
    total = len(finite)
    for rank, (index, value) in reversed(list(enumerate(sorted(finite, key=lambda item: item[1]), 1))):
        running = min(running, value * total / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def build_orthogroup_summary(
    features_path: Path, families_path: Path, output_path: Path
) -> list[dict[str, str]]:
    with Path(families_path).open(encoding="utf-8", newline="") as handle:
        family_sizes = {
            row["orthogroup_id"]: row.get("pmas_count", "")
            for row in csv.DictReader(handle, delimiter="\t")
        }

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    with Path(features_path).open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            orthogroup = row.get("orthogroup_id", "").strip()
            if orthogroup:
                grouped[orthogroup].append(row)

    rows: list[dict[str, str]] = []
    simes_values: list[float | None] = []
    for orthogroup in sorted(grouped):
        members = sorted(grouped[orthogroup], key=lambda row: row["feature_id"])
        pvalues = [value for row in members if (value := _number(row.get("interaction_pvalue", ""))) is not None]
        padj = [value for row in members if (value := _number(row.get("interaction_padj", ""))) is not None]
        effects = [value for row in members if (value := _number(row.get("interaction_log2_fold_change", ""))) is not None]
        significant = [row for row in members if row.get("foreground_850", "").lower() == "yes"]
        significant_effects = [
            value
            for row in significant
            if (value := _number(row.get("interaction_log2_fold_change", ""))) is not None
        ]
        signs = {1 if value > 0 else -1 if value < 0 else 0 for value in significant_effects}
        signs.discard(0)
        direction = "not_estimable"
        if signs == {1}:
            direction = "positive"
        elif signs == {-1}:
            direction = "negative"
        elif len(signs) > 1:
            direction = "mixed"

        modules = sorted(set().union(*(_split(row.get("module_ids", "")) for row in members)))
        axes = sorted(set().union(*(_split(row.get("mechanism_axes", "")) for row in members)))
        family_size = _number(family_sizes.get(orthogroup, ""))
        simes = simes_pvalue(pvalues)
        simes_values.append(simes)
        rows.append(
            {
                "orthogroup_id": orthogroup,
                "member_feature_ids": ";".join(row["feature_id"] for row in members),
                "tested_member_count": str(len(members)),
                "finite_pvalue_member_count": str(len(pvalues)),
                "significant_member_count": str(len(significant)),
                "significant_feature_ids": ";".join(row["feature_id"] for row in significant),
                "min_member_padj": _format(min(padj) if padj else None),
                "max_abs_interaction_lfc": _format(max(map(abs, effects)) if effects else None),
                "direction_status": direction,
                "simes_pvalue": _format(simes),
                "simes_padj": "",
                "module_ids": ";".join(modules),
                "mechanism_axes": ";".join(axes),
                "pmas_family_size": _format(family_size),
                "mapping_completeness": _format(len(members) / family_size if family_size else None),
            }
        )

    for row, adjusted in zip(rows, _bh_adjust(simes_values)):
        row["simes_padj"] = _format(adjusted)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with Path(output_path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--families", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = build_orthogroup_summary(args.features, args.families, args.output)
    print(f"Wrote {len(rows)} orthogroup summaries to {args.output}")


if __name__ == "__main__":
    main()
