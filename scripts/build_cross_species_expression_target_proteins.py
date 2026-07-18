#!/usr/bin/env python3
"""Extract pine protein references for current Tier B orthogroups."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


FIELDS = ["orthogroup_id", "species_id", "gene_id", "fasta_id", "sequence_length", "source_fasta"]


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _read_fasta(path: Path) -> dict[str, str]:
    opener = gzip.open if Path(path).suffix == ".gz" else open
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


def build_target_references(
    shortlist_path: Path,
    genes_path: Path,
    proteomes: dict[str, Path],
    manifest_path: Path,
    pden_output: Path,
    ptae_output: Path,
) -> list[dict[str, str]]:
    followup_orthogroups = {
        row["orthogroup_id"]
        for row in _read_tsv(shortlist_path)
        if row["evidence_tier"] in {"Tier A", "Tier B"}
    }
    genes = [
        row for row in _read_tsv(genes_path)
        if row["orthogroup_id"] in followup_orthogroups and row["species_id"] in proteomes
    ]
    seen_gene: dict[str, str] = {}
    for row in genes:
        previous = seen_gene.setdefault(row["gene_id"], row["orthogroup_id"])
        if previous != row["orthogroup_id"]:
            raise ValueError(f"Protein {row['gene_id']} assigned to conflicting orthogroups")

    sequences = {species: _read_fasta(path) for species, path in proteomes.items()}
    output_records: dict[str, list[tuple[str, str]]] = {species: [] for species in proteomes}
    rows: list[dict[str, str]] = []
    for row in sorted(genes, key=lambda item: (item["orthogroup_id"], item["species_id"], item["gene_id"])):
        species, gene_id = row["species_id"], row["gene_id"]
        fasta_id = gene_id.split("|", 1)[-1]
        sequence = sequences[species].get(fasta_id)
        if sequence is None:
            raise ValueError(f"Protein {gene_id} not found in {proteomes[species]}")
        output_records[species].append((gene_id, sequence))
        rows.append({
            "orthogroup_id": row["orthogroup_id"],
            "species_id": species,
            "gene_id": gene_id,
            "fasta_id": fasta_id,
            "sequence_length": str(len(sequence)),
            "source_fasta": str(proteomes[species]),
        })

    outputs = {"pden": Path(pden_output), "ptae": Path(ptae_output)}
    for species, output in outputs.items():
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", encoding="utf-8", newline="\n") as handle:
            for gene_id, sequence in output_records.get(species, []):
                handle.write(f">{gene_id}\n{sequence}\n")
    Path(manifest_path).parent.mkdir(parents=True, exist_ok=True)
    with Path(manifest_path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shortlist", type=Path, required=True)
    parser.add_argument("--genes", type=Path, required=True)
    parser.add_argument("--pden-proteome", type=Path, required=True)
    parser.add_argument("--ptae-proteome", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--pden-output", type=Path, required=True)
    parser.add_argument("--ptae-output", type=Path, required=True)
    args = parser.parse_args()
    rows = build_target_references(args.shortlist, args.genes, {"pden": args.pden_proteome, "ptae": args.ptae_proteome}, args.manifest, args.pden_output, args.ptae_output)
    counts = {species: len({row["orthogroup_id"] for row in rows if row["species_id"] == species}) for species in ("pden", "ptae")}
    print(f"Wrote {len(rows)} target proteins; orthogroups by species: {counts}")


if __name__ == "__main__":
    main()
