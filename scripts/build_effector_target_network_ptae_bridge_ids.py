#!/usr/bin/env python3
"""Expand Route 1 seed orthogroups to their P. taeda bridge proteins."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


FIELDS = [
    "orthogroup_id",
    "ptae_protein_id",
    "effector_target_network_role",
    "host_module",
    "host_network_role",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_bridge_ids(seeds_path: Path, orthogroups_path: Path, output_path: Path) -> list[dict[str, str]]:
    """Write one row per P. taeda protein in a Route 1 seed orthogroup."""
    seeds = {row["orthogroup_id"]: row for row in read_tsv(seeds_path)}
    orthogroups = {row["Orthogroup"]: row for row in read_tsv(orthogroups_path)}
    rows: list[dict[str, str]] = []
    for orthogroup_id in sorted(seeds):
        seed = seeds[orthogroup_id]
        members = orthogroups.get(orthogroup_id, {}).get("ptae", "")
        for protein_id in sorted(member.strip() for member in members.split(",") if member.strip()):
            rows.append(
                {
                    "orthogroup_id": orthogroup_id,
                    "ptae_protein_id": protein_id,
                    "effector_target_network_role": seed.get("effector_target_network_role", ""),
                    "host_module": seed.get("host_module", ""),
                    "host_network_role": seed.get("host_network_role", ""),
                }
            )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Wrote %d P. taeda bridge proteins for %d Route 1 seed orthogroups", len(rows), len(seeds))
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=Path, required=True)
    parser.add_argument("--orthogroups", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = parse_args()
    build_bridge_ids(args.seeds, args.orthogroups, args.output)


if __name__ == "__main__":
    main()
