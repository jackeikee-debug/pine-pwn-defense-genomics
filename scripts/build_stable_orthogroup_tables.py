#!/usr/bin/env python3
"""Build reusable tables for stable core orthogroups."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


BASE_FIELDS = [
    "orthogroup_id",
    "best_expansion_orthogroup",
    "species_present",
    "total_core_genes",
]
METRIC_FIELDS = ["pilot_core_retention", "core_jaccard"]
GENE_FIELDS = ["orthogroup_id", "species_id", "gene_id"]


def split_genes(value: str) -> list[str]:
    if not value.strip():
        return []
    return [gene.strip() for gene in value.split(",") if gene.strip()]


def read_stable_rows(path: Path) -> dict[str, dict[str, str]]:
    stable_rows: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("status") == "stable":
                stable_rows[row["pilot_orthogroup"]] = row
    return stable_rows


def build_tables(
    orthogroups_path: Path,
    stability_path: Path,
    species_ids: list[str],
    orthogroup_output: Path,
    gene_output: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    stable_rows = read_stable_rows(stability_path)
    count_fields = [f"{species_id}_count" for species_id in species_ids]
    orthogroup_fields = BASE_FIELDS + count_fields + METRIC_FIELDS

    orthogroup_rows: list[dict[str, str]] = []
    gene_rows: list[dict[str, str]] = []

    with orthogroups_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            orthogroup_id = row["Orthogroup"]
            stability = stable_rows.get(orthogroup_id)
            if not stability:
                continue

            genes_by_species = {
                species_id: split_genes(row.get(species_id, ""))
                for species_id in species_ids
            }
            species_present = sum(1 for genes in genes_by_species.values() if genes)
            total_core_genes = sum(len(genes) for genes in genes_by_species.values())

            orthogroup_row = {
                "orthogroup_id": orthogroup_id,
                "best_expansion_orthogroup": stability["best_expansion_orthogroup"],
                "species_present": str(species_present),
                "total_core_genes": str(total_core_genes),
                "pilot_core_retention": stability["pilot_core_retention"],
                "core_jaccard": stability["core_jaccard"],
            }
            for species_id, genes in genes_by_species.items():
                orthogroup_row[f"{species_id}_count"] = str(len(genes))
                for gene_id in genes:
                    gene_rows.append(
                        {
                            "orthogroup_id": orthogroup_id,
                            "species_id": species_id,
                            "gene_id": gene_id,
                        }
                    )
            orthogroup_rows.append(orthogroup_row)

    write_tsv(orthogroup_output, orthogroup_fields, orthogroup_rows)
    write_tsv(gene_output, GENE_FIELDS, gene_rows)
    return orthogroup_rows, gene_rows


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orthogroups", type=Path, required=True)
    parser.add_argument("--stability", type=Path, required=True)
    parser.add_argument("--species", nargs="+", required=True)
    parser.add_argument("--orthogroup-output", type=Path, required=True)
    parser.add_argument("--gene-output", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    orthogroup_rows, gene_rows = build_tables(
        orthogroups_path=args.orthogroups,
        stability_path=args.stability,
        species_ids=args.species,
        orthogroup_output=args.orthogroup_output,
        gene_output=args.gene_output,
    )
    print(f"Stable orthogroups written: {len(orthogroup_rows)}")
    print(f"Stable orthogroup genes written: {len(gene_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
