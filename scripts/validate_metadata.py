#!/usr/bin/env python3
"""Validate project metadata tables for the pine PWN comparative genomics repo."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

try:
    import yaml
except ImportError:  # pragma: no cover - exercised only in minimal environments
    yaml = None


SPECIES_REQUIRED_COLUMNS = [
    "species_id",
    "species_name",
    "common_name",
    "region",
    "group",
    "data_type",
    "protein_fasta",
    "gff3",
    "genome_fasta",
    "source",
    "accession",
    "version",
    "notes",
]

PHENOTYPE_REQUIRED_COLUMNS = [
    "species_id",
    "species_name",
    "pwd_response",
    "response_score",
    "evidence_type",
    "evidence_strength",
    "key_reference",
    "notes",
]

EXPANSION_CANDIDATE_REQUIRED_COLUMNS = [
    "species_id",
    "species_name",
    "common_name",
    "region",
    "phylogenetic_group",
    "role",
    "recommended_use",
    "protein_fasta",
    "source",
    "accession",
    "version",
    "download_status",
    "notes",
]

VALID_RESPONSE_SCORES = {"0", "1", "2", "NA", ""}
MISSING_PATH_VALUES = {"", "NA", "na", "N/A", "n/a"}
DEFAULT_EXPANSION_RECOMMENDED_USES = {"outgroup_or_sensitivity_analysis"}


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def missing_columns(header: list[str], required: list[str]) -> list[str]:
    return [column for column in required if column not in header]


def duplicate_values(rows: list[dict[str, str]], column: str) -> list[str]:
    seen: set[str] = set()
    duplicated: set[str] = set()
    for row in rows:
        value = row.get(column, "").strip()
        if not value:
            continue
        if value in seen:
            duplicated.add(value)
        seen.add(value)
    return sorted(duplicated)


def validate_species_metadata(path: str | Path, strict_paths: bool = False, base_dir: str | Path = ".") -> list[str]:
    path = Path(path)
    base_dir = Path(base_dir)
    errors: list[str] = []
    header, rows = read_tsv(path)

    for column in missing_columns(header, SPECIES_REQUIRED_COLUMNS):
        errors.append(f"Missing required column in {path}: {column}")

    for species_id in duplicate_values(rows, "species_id"):
        errors.append(f"Duplicate species_id in {path}: {species_id}")

    if strict_paths and "protein_fasta" in header:
        for index, row in enumerate(rows, start=2):
            protein_fasta = row.get("protein_fasta", "").strip()
            if protein_fasta in MISSING_PATH_VALUES:
                continue
            fasta_path = Path(protein_fasta)
            if not fasta_path.is_absolute():
                fasta_path = base_dir / fasta_path
            if not fasta_path.exists():
                errors.append(f"Missing protein_fasta file in {path} line {index}: {protein_fasta}")

    return errors


def validate_phenotype_matrix(path: str | Path) -> list[str]:
    path = Path(path)
    errors: list[str] = []
    header, rows = read_tsv(path)

    for column in missing_columns(header, PHENOTYPE_REQUIRED_COLUMNS):
        errors.append(f"Missing required column in {path}: {column}")

    for species_id in duplicate_values(rows, "species_id"):
        errors.append(f"Duplicate species_id in {path}: {species_id}")

    for index, row in enumerate(rows, start=2):
        score = row.get("response_score", "").strip()
        if score not in VALID_RESPONSE_SCORES:
            errors.append(f"Invalid response_score in {path} line {index}: {score}")

    return errors


def validate_expansion_candidates(
    path: str | Path,
    strict_paths: bool = False,
    base_dir: str | Path = ".",
) -> list[str]:
    path = Path(path)
    base_dir = Path(base_dir)
    errors: list[str] = []
    header, rows = read_tsv(path)

    for column in missing_columns(header, EXPANSION_CANDIDATE_REQUIRED_COLUMNS):
        errors.append(f"Missing required column in {path}: {column}")

    for species_id in duplicate_values(rows, "species_id"):
        errors.append(f"Duplicate species_id in {path}: {species_id}")

    if strict_paths and "protein_fasta" in header:
        for index, row in enumerate(rows, start=2):
            if row.get("download_status", "").strip() != "downloaded":
                continue
            protein_fasta = row.get("protein_fasta", "").strip()
            if protein_fasta in MISSING_PATH_VALUES:
                errors.append(f"Downloaded expansion candidate missing protein_fasta in {path} line {index}")
                continue
            fasta_path = Path(protein_fasta)
            if not fasta_path.is_absolute():
                fasta_path = base_dir / fasta_path
            if not fasta_path.exists():
                errors.append(f"Missing expansion protein_fasta file in {path} line {index}: {protein_fasta}")

    return errors


def selected_expansion_species(path: str | Path) -> list[dict[str, str]]:
    _, rows = read_tsv(Path(path))
    return [
        row
        for row in rows
        if row.get("download_status", "").strip() == "downloaded"
        and row.get("recommended_use", "").strip() in DEFAULT_EXPANSION_RECOMMENDED_USES
    ]


def load_config(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        if yaml is not None:
            return yaml.safe_load(handle)
        return load_minimal_metadata_config(handle.readlines())


def load_minimal_metadata_config(lines: list[str]) -> dict:
    """Read only the simple metadata mapping when PyYAML is unavailable."""
    metadata: dict[str, str] = {}
    in_metadata = False

    for raw_line in lines:
        line = raw_line.rstrip()
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if not line.startswith(" ") and stripped.endswith(":"):
            in_metadata = stripped == "metadata:"
            continue
        if in_metadata and line.startswith("  ") and ":" in stripped:
            key, value = stripped.split(":", 1)
            metadata[key.strip()] = value.strip().strip("'\"")

    return {"metadata": metadata}


def validate_from_config(config_path: str | Path, strict_paths: bool = False) -> list[str]:
    config_path = Path(config_path)
    config = load_config(config_path)
    metadata = config.get("metadata", {})
    errors: list[str] = []

    species_path = Path(metadata.get("species", "config/species_metadata.tsv"))
    phenotype_path = Path(metadata.get("phenotype", "config/phenotype_matrix.tsv"))
    expansion_candidates_value = metadata.get("expansion_candidates", "").strip()

    errors.extend(validate_species_metadata(species_path, strict_paths=strict_paths, base_dir=Path.cwd()))
    errors.extend(validate_phenotype_matrix(phenotype_path))
    if expansion_candidates_value:
        expansion_candidates_path = Path(expansion_candidates_value)
        errors.extend(
            validate_expansion_candidates(
                expansion_candidates_path,
                strict_paths=strict_paths,
                base_dir=Path.cwd(),
            )
        )
    return errors


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/config.yaml", help="Path to project config YAML.")
    parser.add_argument("--strict-paths", action="store_true", help="Require non-placeholder FASTA paths to exist.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    errors = validate_from_config(args.config, strict_paths=args.strict_paths)

    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1

    print("Metadata validation passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
