#!/usr/bin/env python3
"""Build a multi-species protein reference for Tier A/B candidate orthogroups."""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
from pathlib import Path


FIELDS = [
    "orthogroup_id",
    "species_id",
    "gene_id",
    "fasta_id",
    "sequence_length",
    "source_fasta",
]


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _read_fasta(path: Path) -> dict[str, str]:
    opener = gzip.open if path.suffix == ".gz" else open
    records: dict[str, str] = {}
    current = ""
    sequence: list[str] = []
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    if current in records:
                        raise ValueError(f"Duplicate FASTA ID: {current}")
                    records[current] = "".join(sequence)
                current = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line.strip())
    if current:
        if current in records:
            raise ValueError(f"Duplicate FASTA ID: {current}")
        records[current] = "".join(sequence)
    return records


def _resolve_fasta(value: str, metadata_path: Path) -> Path:
    path = Path(value)
    if path.is_absolute() or path.exists():
        return path
    project_relative = metadata_path.parent.parent / path
    return project_relative


def build_reference(
    shortlist_path: Path,
    genes_path: Path,
    metadata_path: Path,
    fasta_output: Path,
    manifest_output: Path,
) -> list[dict[str, str]]:
    """Extract all stable-proteome members of Tier A/B candidate orthogroups."""
    candidates = {
        row["orthogroup_id"]
        for row in _read_tsv(shortlist_path)
        if row.get("evidence_tier") in {"Tier A", "Tier B"}
        and row.get("candidate_level", "orthogroup") == "orthogroup"
    }
    proteomes = {
        row["species_id"]: _resolve_fasta(row["protein_fasta"], Path(metadata_path))
        for row in _read_tsv(metadata_path)
        if row.get("protein_fasta", "") not in {"", "NA"}
        and row.get("group", "") != "outgroup"
    }
    selected = sorted(
        (
            row
            for row in _read_tsv(genes_path)
            if row["orthogroup_id"] in candidates and row["species_id"] in proteomes
        ),
        key=lambda row: (row["orthogroup_id"], row["species_id"], row["gene_id"]),
    )
    sequence_cache = {species: _read_fasta(path) for species, path in proteomes.items()}
    rows: list[dict[str, str]] = []
    records: list[tuple[str, str]] = []
    seen_genes: set[str] = set()
    for row in selected:
        gene_id = row["gene_id"]
        if gene_id in seen_genes:
            raise ValueError(f"Duplicate stable gene assignment: {gene_id}")
        seen_genes.add(gene_id)
        species = row["species_id"]
        fasta_id = gene_id.split("|", 1)[-1]
        sequence = sequence_cache[species].get(fasta_id)
        if sequence is None:
            raise ValueError(f"Stable gene {gene_id} missing from {proteomes[species]}")
        records.append((gene_id, sequence))
        rows.append(
            {
                "orthogroup_id": row["orthogroup_id"],
                "species_id": species,
                "gene_id": gene_id,
                "fasta_id": fasta_id,
                "sequence_length": str(len(sequence)),
                "source_fasta": str(proteomes[species]),
            }
        )

    Path(fasta_output).parent.mkdir(parents=True, exist_ok=True)
    with Path(fasta_output).open("w", encoding="utf-8", newline="\n") as handle:
        for gene_id, sequence in records:
            handle.write(f">{gene_id}\n{sequence}\n")
    Path(manifest_output).parent.mkdir(parents=True, exist_ok=True)
    with Path(manifest_output).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Wrote %d proteins from %d candidate orthogroups", len(rows), len({row["orthogroup_id"] for row in rows}))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shortlist", type=Path, required=True)
    parser.add_argument("--genes", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--fasta-output", type=Path, required=True)
    parser.add_argument("--manifest-output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    build_reference(args.shortlist, args.genes, args.metadata, args.fasta_output, args.manifest_output)


if __name__ == "__main__":
    main()
