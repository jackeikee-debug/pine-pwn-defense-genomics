#!/usr/bin/env python3
"""Prepare CAFE/CAFE5 pilot gene-family count inputs from OrthoFinder output."""

from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path


SUMMARY_FIELDS = [
    "total_orthogroups",
    "written_families",
    "filtered_by_max_family_size",
    "filtered_by_min_species_present",
    "max_family_size",
    "min_species_present",
    "species",
    "tree_source",
    "tree_note",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def prepare_cafe_inputs(
    orthogroups_path: Path,
    manifest_path: Path,
    tree_path: Path,
    output_path: Path,
    summary_path: Path,
    tree_output_path: Path,
    max_family_size: int = 100,
    min_species_present: int = 2,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    species = [row["species_id"] for row in read_tsv(manifest_path)]
    count_rows: list[dict[str, str]] = []
    filtered_by_max = 0
    filtered_by_min_present = 0
    total = 0

    with orthogroups_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            total += 1
            counts = {species_id: count_members(row.get(species_id, "")) for species_id in species}
            family_size = sum(counts.values())
            species_present = sum(1 for count in counts.values() if count > 0)
            if max_family_size > 0 and family_size > max_family_size:
                filtered_by_max += 1
                continue
            if species_present < min_species_present:
                filtered_by_min_present += 1
                continue
            count_rows.append(
                {
                    "Desc": "(null)",
                    "Family ID": row["Orthogroup"],
                    **{species_id: str(counts[species_id]) for species_id in species},
                }
            )

    write_tsv(output_path, ["Desc", "Family ID", *species], count_rows)
    copy_tree(tree_path, tree_output_path)

    summary_rows = [
        {
            "total_orthogroups": str(total),
            "written_families": str(len(count_rows)),
            "filtered_by_max_family_size": str(filtered_by_max),
            "filtered_by_min_species_present": str(filtered_by_min_present),
            "max_family_size": str(max_family_size),
            "min_species_present": str(min_species_present),
            "species": ";".join(species),
            "tree_source": str(tree_path),
            "tree_note": "OrthoFinder rooted species tree for pilot only; replace with time-calibrated ultrametric tree for formal CAFE inference.",
        }
    ]
    write_tsv(summary_path, SUMMARY_FIELDS, summary_rows)
    return count_rows, summary_rows


def count_members(value: str) -> int:
    if not value.strip():
        return 0
    return len([part for part in value.split(",") if part.strip()])


def copy_tree(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    text = source.read_text(encoding="utf-8").strip()
    destination.write_text(f"{text}\n", encoding="utf-8")


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orthogroups", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--tree", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--tree-output", type=Path, required=True)
    parser.add_argument("--max-family-size", type=int, default=100)
    parser.add_argument("--min-species-present", type=int, default=2)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows, summary_rows = prepare_cafe_inputs(
        orthogroups_path=args.orthogroups,
        manifest_path=args.manifest,
        tree_path=args.tree,
        output_path=args.output,
        summary_path=args.summary,
        tree_output_path=args.tree_output,
        max_family_size=args.max_family_size,
        min_species_present=args.min_species_present,
    )
    print(f"CAFE pilot families written: {len(rows)}")
    print(f"CAFE pilot summary rows written: {len(summary_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
