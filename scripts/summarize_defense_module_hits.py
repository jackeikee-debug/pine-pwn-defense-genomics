#!/usr/bin/env python3
"""Summarize defense-module hits by species and module."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


SPECIES_FIELDS = [
    "module_id",
    "species_id",
    "pwd_response",
    "response_score",
    "orthogroup_count",
    "gene_count",
    "matched_keywords",
    "annotation_sources",
    "evidence_types",
]
MODULE_FIELDS = [
    "module_id",
    "species_count",
    "orthogroup_count",
    "gene_count",
    "species_ids",
    "matched_keywords",
    "annotation_sources",
    "evidence_types",
]


def split_values(value: str) -> list[str]:
    values: list[str] = []
    for separator in [";", ","]:
        if separator in value:
            return [part.strip() for part in value.split(separator) if part.strip()]
    if value.strip():
        values.append(value.strip())
    return values


def read_phenotypes(path: Path) -> dict[str, dict[str, str]]:
    phenotypes: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            phenotypes[row["species_id"]] = row
    return phenotypes


def summarize_hits(
    matrix_path: Path,
    phenotype_path: Path,
    species_output: Path,
    module_output: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    phenotypes = read_phenotypes(phenotype_path)
    by_species: dict[tuple[str, str], dict[str, object]] = defaultdict(new_bucket)
    by_module: dict[str, dict[str, object]] = defaultdict(new_bucket)

    with matrix_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            module_id = row["module_id"]
            species_id = row["species_id"]
            gene_count = int(row.get("gene_count", "0") or 0)

            species_bucket = by_species[(module_id, species_id)]
            update_bucket(species_bucket, row, gene_count, species_id)

            module_bucket = by_module[module_id]
            update_bucket(module_bucket, row, gene_count, species_id)

    species_rows: list[dict[str, str]] = []
    for module_id, species_id in sorted(by_species):
        bucket = by_species[(module_id, species_id)]
        phenotype = phenotypes.get(species_id, {})
        species_rows.append(
            {
                "module_id": module_id,
                "species_id": species_id,
                "pwd_response": phenotype.get("pwd_response", "NA"),
                "response_score": phenotype.get("response_score", "NA"),
                "orthogroup_count": str(len(bucket["orthogroups"])),
                "gene_count": str(bucket["gene_count"]),
                "matched_keywords": join_set(bucket["matched_keywords"]),
                "annotation_sources": join_set(bucket["annotation_sources"]),
                "evidence_types": join_set(bucket["evidence_types"]),
            }
        )

    module_rows: list[dict[str, str]] = []
    for module_id in sorted(by_module):
        bucket = by_module[module_id]
        module_rows.append(
            {
                "module_id": module_id,
                "species_count": str(len(bucket["species_ids"])),
                "orthogroup_count": str(len(bucket["orthogroups"])),
                "gene_count": str(bucket["gene_count"]),
                "species_ids": join_set(bucket["species_ids"]),
                "matched_keywords": join_set(bucket["matched_keywords"]),
                "annotation_sources": join_set(bucket["annotation_sources"]),
                "evidence_types": join_set(bucket["evidence_types"]),
            }
        )

    write_tsv(species_output, SPECIES_FIELDS, species_rows)
    write_tsv(module_output, MODULE_FIELDS, module_rows)
    return species_rows, module_rows


def new_bucket() -> dict[str, object]:
    return {
        "orthogroups": set(),
        "species_ids": set(),
        "gene_count": 0,
        "matched_keywords": set(),
        "annotation_sources": set(),
        "evidence_types": set(),
    }


def update_bucket(bucket: dict[str, object], row: dict[str, str], gene_count: int, species_id: str) -> None:
    bucket["orthogroups"].add(row["orthogroup_id"])
    bucket["species_ids"].add(species_id)
    bucket["gene_count"] += gene_count
    bucket["matched_keywords"].update(split_values(row.get("matched_keywords", "")))
    bucket["annotation_sources"].update(split_values(row.get("annotation_sources", "")))
    bucket["evidence_types"].update(split_values(row.get("evidence_type", "")))


def join_set(value: object) -> str:
    return ";".join(sorted(value))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--phenotype", type=Path, required=True)
    parser.add_argument("--species-output", type=Path, required=True)
    parser.add_argument("--module-output", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    species_rows, module_rows = summarize_hits(
        matrix_path=args.matrix,
        phenotype_path=args.phenotype,
        species_output=args.species_output,
        module_output=args.module_output,
    )
    print(f"Defense module species summary rows written: {len(species_rows)}")
    print(f"Defense module summary rows written: {len(module_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
