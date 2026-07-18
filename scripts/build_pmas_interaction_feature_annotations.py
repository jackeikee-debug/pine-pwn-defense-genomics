#!/usr/bin/env python3
"""Annotate P. massoniana interaction features for mechanism competition."""

from __future__ import annotations

import argparse
import csv
import logging
import math
from collections import defaultdict
from pathlib import Path


LOGGER = logging.getLogger(__name__)
ADDED_FIELDS = [
    "orthogroup_id",
    "orthogroup_mapping_status",
    "module_ids",
    "mechanism_axes",
    "annotation_sources",
    "finite_wald",
    "finite_padj",
    "foreground_850",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def finite_number(value: str | None) -> float | None:
    if value in (None, "", "NA", "NaN"):
        return None
    try:
        number = float(value)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def read_mechanism_modules(path: Path) -> dict[str, set[str]]:
    module_to_axes: dict[str, set[str]] = defaultdict(set)
    for row in read_tsv(path):
        axis = row["mechanism_axis"].strip()
        for module in filter(None, (value.strip() for value in row.get("host_modules", "").split(";"))):
            module_to_axes[module].add(axis)
    return dict(module_to_axes)


def index_orthogroups(path: Path) -> dict[str, str]:
    feature_to_orthogroup: dict[str, str] = {}
    for row in read_tsv(path):
        if row.get("species_id") != "pmas":
            continue
        feature_id = row["gene_id"].removeprefix("pmas|")
        orthogroup = row["orthogroup_id"]
        previous = feature_to_orthogroup.get(feature_id)
        if previous is not None and previous != orthogroup:
            raise ValueError(
                f"Feature {feature_id} has conflicting exact orthogroup assignments: {previous}, {orthogroup}"
            )
        feature_to_orthogroup[feature_id] = orthogroup
    return feature_to_orthogroup


def index_defense_hits(path: Path) -> dict[str, dict[str, set[str]]]:
    by_orthogroup: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: {"modules": set(), "sources": set()}
    )
    for row in read_tsv(path):
        orthogroup = row.get("orthogroup_id", "").strip()
        if not orthogroup:
            continue
        module = row.get("module_id", "").strip()
        source = row.get("annotation_source", "").strip()
        if module:
            by_orthogroup[orthogroup]["modules"].add(module)
        if source:
            by_orthogroup[orthogroup]["sources"].add(source)
    return dict(by_orthogroup)


def build_feature_annotations(
    interaction_path: Path,
    orthogroup_genes_path: Path,
    defense_hits_path: Path,
    mechanism_axes_path: Path,
    output_path: Path,
) -> list[dict[str, str]]:
    interactions = read_tsv(interaction_path)
    feature_to_orthogroup = index_orthogroups(orthogroup_genes_path)
    defense_by_orthogroup = index_defense_hits(defense_hits_path)
    module_to_axes = read_mechanism_modules(mechanism_axes_path)

    seen: set[str] = set()
    output_rows: list[dict[str, str]] = []
    for interaction in interactions:
        feature_id = interaction["feature_id"]
        if feature_id in seen:
            raise ValueError(f"Duplicate interaction feature_id: {feature_id}")
        seen.add(feature_id)

        orthogroup = feature_to_orthogroup.get(feature_id, "")
        defense = defense_by_orthogroup.get(orthogroup, {"modules": set(), "sources": set()})
        modules = sorted(defense["modules"])
        axes = sorted({axis for module in modules for axis in module_to_axes.get(module, set())})
        statistic = finite_number(interaction.get("interaction_stat"))
        adjusted_p = finite_number(interaction.get("interaction_padj"))
        fold_change = finite_number(interaction.get("interaction_log2_fold_change"))
        foreground = adjusted_p is not None and adjusted_p < 0.05 and fold_change is not None and abs(fold_change) >= 1

        row = dict(interaction)
        row.update(
            {
                "orthogroup_id": orthogroup,
                "orthogroup_mapping_status": "exact_stable_orthogroup" if orthogroup else "unmapped",
                "module_ids": ";".join(modules),
                "mechanism_axes": ";".join(axes),
                "annotation_sources": ";".join(sorted(defense["sources"])),
                "finite_wald": "yes" if statistic is not None else "no",
                "finite_padj": "yes" if adjusted_p is not None else "no",
                "foreground_850": "yes" if foreground else "no",
            }
        )
        output_rows.append(row)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(interactions[0]) + ADDED_FIELDS if interactions else ADDED_FIELDS
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)
    return output_rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--interaction", type=Path, required=True)
    parser.add_argument("--orthogroup-genes", type=Path, required=True)
    parser.add_argument("--defense-hits", type=Path, required=True)
    parser.add_argument("--mechanism-axes", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = build_feature_annotations(
        args.interaction,
        args.orthogroup_genes,
        args.defense_hits,
        args.mechanism_axes,
        args.output,
    )
    LOGGER.info("Wrote %d feature annotations to %s", len(rows), args.output)
    return 0


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    raise SystemExit(main())
