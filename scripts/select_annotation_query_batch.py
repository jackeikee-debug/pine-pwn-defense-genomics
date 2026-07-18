#!/usr/bin/env python3
"""Select a small stable-protein batch for external annotation."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


MANIFEST_FIELDS = ["orthogroup_id", "species_id", "gene_id", "query_fasta"]


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    header: str | None = None
    sequence_parts: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(sequence_parts)
                header = line[1:].split()[0]
                sequence_parts = []
            else:
                sequence_parts.append(line)
    if header is not None:
        records[header] = "".join(sequence_parts)
    return records


def wrap_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[index : index + width] for index in range(0, len(sequence), width))


def select_query_batch(
    stable_genes_path: Path,
    proteome_dir: Path,
    output_fasta: Path,
    manifest_output: Path,
    max_total: int,
    max_per_species: int,
) -> list[dict[str, str]]:
    proteomes: dict[str, dict[str, str]] = {}
    species_counts: defaultdict[str, int] = defaultdict(int)
    selected_rows: list[dict[str, str]] = []
    selected_sequences: list[tuple[str, str]] = []

    with stable_genes_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if max_total > 0 and len(selected_rows) >= max_total:
                break
            species_id = row["species_id"]
            if max_per_species > 0 and species_counts[species_id] >= max_per_species:
                continue
            if species_id not in proteomes:
                proteomes[species_id] = read_fasta(proteome_dir / f"{species_id}.faa")
            sequence = proteomes[species_id].get(row["gene_id"])
            if not sequence:
                continue

            species_counts[species_id] += 1
            selected_sequences.append((row["gene_id"], sequence))
            selected_rows.append(
                {
                    "orthogroup_id": row["orthogroup_id"],
                    "species_id": species_id,
                    "gene_id": row["gene_id"],
                    "query_fasta": str(output_fasta),
                }
            )

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with output_fasta.open("w", encoding="utf-8") as handle:
        for gene_id, sequence in selected_sequences:
            handle.write(f">{gene_id}\n{wrap_sequence(sequence)}\n")
    write_tsv(manifest_output, MANIFEST_FIELDS, selected_rows)
    return selected_rows


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stable-genes", type=Path, required=True)
    parser.add_argument("--proteome-dir", type=Path, required=True)
    parser.add_argument("--output-fasta", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--max-total", type=int, default=600)
    parser.add_argument("--max-per-species", type=int, default=100)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = select_query_batch(
        stable_genes_path=args.stable_genes,
        proteome_dir=args.proteome_dir,
        output_fasta=args.output_fasta,
        manifest_output=args.manifest,
        max_total=args.max_total,
        max_per_species=args.max_per_species,
    )
    print(f"Selected annotation query proteins: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
