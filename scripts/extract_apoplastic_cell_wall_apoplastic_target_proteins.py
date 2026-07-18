#!/usr/bin/env python3
"""Extract all pine members of Tier A Route 2 apoplastic target orthogroups."""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
from collections import defaultdict
from pathlib import Path


FIELDS = [
    "orthogroup_id", "species_id", "gene_id", "fasta_id", "sequence_length",
    "expression_selected_member", "expression_evidence_tier", "source_fasta",
]
VALID_AA = set("ABCDEFGHIKLMNPQRSTVWXYZJUO*")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def scan_fasta(path: Path, requested: set[str]) -> dict[str, str]:
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    found: dict[str, str] = {}
    current = ""
    sequence: list[str] = []

    def store() -> None:
        if not current or current not in requested:
            return
        if current in found:
            raise ValueError(f"Duplicate requested FASTA ID {current} in {path}")
        value = "".join(sequence).upper()
        invalid = sorted(set(value) - VALID_AA)
        if invalid:
            raise ValueError(f"Invalid amino-acid characters for {current}: {''.join(invalid)}")
        found[current] = value

    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                store()
                current = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line.strip())
        store()
    return found


def extract_targets(
    shortlist_path: Path,
    genes_path: Path,
    species_path: Path,
    expression_path: Path,
    fasta_output: Path,
    manifest_output: Path,
) -> list[dict[str, str]]:
    orthogroups = {
        row["host_orthogroup_id"] for row in read_tsv(shortlist_path)
        if row.get("target_priority_tier") == "A"
    }
    expression_features: dict[str, set[str]] = defaultdict(set)
    expression_tier: dict[str, str] = {}
    for row in read_tsv(expression_path):
        if row.get("orthogroup_id") not in orthogroups:
            continue
        expression_tier[row["orthogroup_id"]] = row.get("evidence_tier", "")
        expression_features[row["orthogroup_id"]].update(
            item.strip() for item in row.get("feature_ids", "").split(";") if item.strip()
        )

    proteins = [row for row in read_tsv(genes_path) if row["orthogroup_id"] in orthogroups]
    sources = {
        row["species_id"]: Path(row["protein_fasta"])
        for row in read_tsv(species_path)
        if row.get("protein_fasta", "") not in {"", "NA"}
    }
    requested: dict[str, set[str]] = defaultdict(set)
    for row in proteins:
        if row["species_id"] not in sources:
            raise ValueError(f"No protein FASTA configured for {row['species_id']}")
        requested[row["species_id"]].add(row["gene_id"].split("|", 1)[-1])

    sequences: dict[str, dict[str, str]] = {}
    for species_id, ids in requested.items():
        if not sources[species_id].exists():
            raise FileNotFoundError(sources[species_id])
        sequences[species_id] = scan_fasta(sources[species_id], ids)
        missing = sorted(ids - set(sequences[species_id]))
        if missing:
            raise ValueError(f"Missing {len(missing)} proteins in {sources[species_id]}: {', '.join(missing[:5])}")

    rows: list[dict[str, str]] = []
    fasta_records: list[tuple[str, str, str, bool]] = []
    for protein in sorted(proteins, key=lambda row: (row["orthogroup_id"], row["species_id"], row["gene_id"])):
        orthogroup = protein["orthogroup_id"]
        species_id = protein["species_id"]
        fasta_id = protein["gene_id"].split("|", 1)[-1]
        sequence = sequences[species_id][fasta_id]
        selected = fasta_id in expression_features.get(orthogroup, set())
        rows.append({
            "orthogroup_id": orthogroup, "species_id": species_id, "gene_id": protein["gene_id"],
            "fasta_id": fasta_id, "sequence_length": str(len(sequence)),
            "expression_selected_member": str(selected).lower(),
            "expression_evidence_tier": expression_tier.get(orthogroup, "") if selected else "",
            "source_fasta": str(sources[species_id]),
        })
        fasta_records.append((protein["gene_id"], orthogroup, sequence, selected))

    fasta_output.parent.mkdir(parents=True, exist_ok=True)
    with fasta_output.open("w", encoding="utf-8", newline="\n") as handle:
        for gene_id, orthogroup, sequence, selected in fasta_records:
            handle.write(f">{gene_id} orthogroup={orthogroup} expression_selected={str(selected).lower()}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")
    manifest_output.parent.mkdir(parents=True, exist_ok=True)
    with manifest_output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Extracted %d proteins from %d species", len(rows), len(requested))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shortlist", type=Path, required=True)
    parser.add_argument("--genes", type=Path, required=True)
    parser.add_argument("--species", type=Path, required=True)
    parser.add_argument("--expression", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    extract_targets(args.shortlist, args.genes, args.species, args.expression, args.fasta, args.manifest)


if __name__ == "__main__":
    main()
