#!/usr/bin/env python3
"""Prepare predeclared occupancy-filtered CAFE5 diagnostic inputs."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path

from Bio import Phylo


SUMMARY_FIELDS = [
    "panel_id", "species", "root_group_a", "root_group_b", "min_species_present",
    "max_total_size", "max_species_count", "total_orthogroups", "written_families",
    "filtered_missing_root_side", "filtered_by_min_species_present",
    "filtered_by_max_total_size", "filtered_by_max_species_count",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def count_members(value: str) -> int:
    return sum(1 for member in value.split(",") if member.strip())


def validate_panel(species: list[str], root_group_a: list[str], root_group_b: list[str]) -> None:
    if not species or len(species) != len(set(species)):
        raise ValueError("Panel species must be non-empty and unique")
    if not root_group_a or not root_group_b:
        raise ValueError("Both root groups must be non-empty")
    overlap = set(root_group_a) & set(root_group_b)
    if overlap:
        raise ValueError(f"Root groups overlap: {', '.join(sorted(overlap))}")
    outside = (set(root_group_a) | set(root_group_b)) - set(species)
    if outside:
        raise ValueError(f"Root-group species outside panel: {', '.join(sorted(outside))}")


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def prepare_diagnostic_inputs(
    orthogroups_path: Path,
    species: list[str],
    root_group_a: list[str],
    root_group_b: list[str],
    min_species_present: int,
    max_total_size: int,
    max_species_count: int,
    counts_output: Path,
    summary_output: Path,
    panel_id: str,
) -> tuple[list[dict[str, str]], dict[str, str]]:
    validate_panel(species, root_group_a, root_group_b)
    if min_species_present < 2 or min_species_present > len(species):
        raise ValueError("min_species_present must be between 2 and panel size")
    if max_total_size <= 0 or max_species_count <= 0:
        raise ValueError("Family-size ceilings must be positive")

    retained: list[dict[str, str]] = []
    excluded = {
        "filtered_missing_root_side": 0,
        "filtered_by_min_species_present": 0,
        "filtered_by_max_total_size": 0,
        "filtered_by_max_species_count": 0,
    }
    total = 0
    with orthogroups_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [name for name in ["Orthogroup", *species] if name not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"Missing OrthoFinder columns: {', '.join(missing)}")
        for row in reader:
            total += 1
            counts = {name: count_members(row.get(name, "")) for name in species}
            if not any(counts[name] > 0 for name in root_group_a) or not any(counts[name] > 0 for name in root_group_b):
                excluded["filtered_missing_root_side"] += 1
                continue
            present = sum(value > 0 for value in counts.values())
            if present < min_species_present:
                excluded["filtered_by_min_species_present"] += 1
                continue
            total_size = sum(counts.values())
            if total_size > max_total_size:
                excluded["filtered_by_max_total_size"] += 1
                continue
            if max(counts.values()) > max_species_count:
                excluded["filtered_by_max_species_count"] += 1
                continue
            retained.append({
                "Desc": "(null)", "Family ID": row["Orthogroup"],
                **{name: str(counts[name]) for name in species},
            })

    write_tsv(counts_output, ["Desc", "Family ID", *species], retained)
    summary = {
        "panel_id": panel_id, "species": ";".join(species),
        "root_group_a": ";".join(root_group_a), "root_group_b": ";".join(root_group_b),
        "min_species_present": str(min_species_present), "max_total_size": str(max_total_size),
        "max_species_count": str(max_species_count), "total_orthogroups": str(total),
        "written_families": str(len(retained)),
        **{key: str(value) for key, value in excluded.items()},
    }
    write_tsv(summary_output, SUMMARY_FIELDS, [summary])
    logging.info("Panel %s retained %d of %d families", panel_id, len(retained), total)
    return retained, summary


def prune_time_tree(source: Path, retained_species: list[str], output: Path) -> None:
    tree = Phylo.read(source, "newick")
    source_species = {terminal.name for terminal in tree.get_terminals()}
    missing = sorted(set(retained_species) - source_species)
    if missing:
        raise ValueError(f"Species missing from time tree: {', '.join(missing)}")
    for terminal in list(tree.get_terminals()):
        if terminal.name not in set(retained_species):
            tree.prune(terminal)
    observed = {terminal.name for terminal in tree.get_terminals()}
    if observed != set(retained_species):
        raise ValueError("Pruned time tree has unexpected terminal membership")
    distances = [tree.distance(terminal) for terminal in tree.get_terminals()]
    delta = max(distances) - min(distances)
    if delta >= 1e-6:
        raise ValueError(f"Pruned time tree is not ultrametric: root-to-tip delta={delta}")
    output.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, output, "newick")
    text = output.read_text(encoding="utf-8").strip()
    output.write_text(text + "\n", encoding="utf-8", newline="\n")


def split_list(value: str) -> list[str]:
    return [item.strip() for item in value.split(";") if item.strip()]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orthogroups", type=Path, required=True)
    parser.add_argument("--tree", type=Path, required=True)
    parser.add_argument("--panels", type=Path, required=True)
    parser.add_argument("--panel-id", required=True)
    parser.add_argument("--counts-output", type=Path, required=True)
    parser.add_argument("--tree-output", type=Path, required=True)
    parser.add_argument("--summary-output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    matches = [row for row in read_tsv(args.panels) if row["panel_id"] == args.panel_id]
    if len(matches) != 1:
        raise ValueError(f"Expected one panel row for {args.panel_id}, found {len(matches)}")
    panel = matches[0]
    species = split_list(panel["species"])
    prepare_diagnostic_inputs(
        args.orthogroups, species, split_list(panel["root_group_a"]),
        split_list(panel["root_group_b"]), int(panel["min_species_present"]),
        int(panel["max_total_size"]), int(panel["max_species_count"]),
        args.counts_output, args.summary_output, args.panel_id,
    )
    prune_time_tree(args.tree, species, args.tree_output)


if __name__ == "__main__":
    main()
