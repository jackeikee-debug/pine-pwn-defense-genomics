#!/usr/bin/env python3
"""Generate proposed mature-protein constructs for Route 2 Tier 1 assays."""

from __future__ import annotations

import argparse
import csv
import logging
import re
from pathlib import Path


FIELDS = [
    "construct_id", "source_protein_id", "protein_role", "host_family", "construct_type",
    "source_start", "source_end", "construct_length", "mutation", "design_status", "design_note",
]


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    current = ""
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                records[current] = ""
            else:
                records[current] += line.upper()
    return records


def safe_id(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def build_constructs(
    fasta_path: Path,
    audit_path: Path,
    output_fasta: Path,
    manifest_path: Path,
) -> list[dict[str, str]]:
    sequences = read_fasta(fasta_path)
    with audit_path.open(encoding="utf-8", newline="") as handle:
        audit = list(csv.DictReader(handle, delimiter="\t"))
    rows: list[dict[str, str]] = []
    for row in audit:
        protein_id = row["protein_id"]
        sequence = sequences[protein_id]
        signal_end = int(row["signal_peptide_end"])
        start = signal_end + 1
        mature = sequence[signal_end:]
        is_effector = row["protein_role"] == "candidate_effector"
        construct_type = "wild_type_effector" if is_effector else "wild_type_substrate"
        base = {
            "source_protein_id": protein_id, "protein_role": row["protein_role"],
            "host_family": row["host_family"], "source_start": str(start),
            "source_end": str(len(sequence)), "construct_length": str(len(mature)),
            "design_status": "proposed_not_experimentally_validated",
        }
        rows.append({
            **base, "construct_id": f"{safe_id(protein_id)}_residues_{start}_{len(sequence)}_WT",
            "construct_type": construct_type, "mutation": "none", "sequence": mature,
            "design_note": (
                "Predicted signal peptide removed; protease activation boundary remains unresolved"
                if is_effector else "Predicted signal peptide removed"
            ),
        })
        mutation = row.get("proposed_catalytic_dead_change", "")
        if is_effector and mutation.endswith("_proposed"):
            match = re.fullmatch(r"([A-Z])(\d+)([A-Z])_proposed", mutation)
            if not match:
                raise ValueError(f"Unparseable proposed mutation: {mutation}")
            reference, position_text, alternate = match.groups()
            position = int(position_text)
            if sequence[position - 1] != reference:
                raise ValueError(f"Mutation {mutation} does not match {protein_id} sequence")
            mutant_full = sequence[: position - 1] + alternate + sequence[position:]
            mutant_mature = mutant_full[signal_end:]
            rows.append({
                **base, "construct_id": f"{safe_id(protein_id)}_residues_{start}_{len(sequence)}_{reference}{position}{alternate}_proposed",
                "construct_type": "proposed_catalytic_dead_effector", "mutation": mutation,
                "sequence": mutant_mature,
                "design_note": "Computationally proposed HEXXH catalytic-glutamate mutant; verify against curated M8 evidence before synthesis",
            })
    rows.sort(key=lambda row: (row["protein_role"] != "candidate_effector", row["construct_type"], row["construct_id"]))
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with output_fasta.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(f'>{row["construct_id"]} source={row["source_protein_id"]} status={row["design_status"]}\n')
            sequence = row["sequence"]
            for index in range(0, len(sequence), 60):
                handle.write(sequence[index:index + 60] + "\n")
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row[field] for field in FIELDS} for row in rows)
    logging.info("Wrote %d proposed constructs", len(rows))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--output-fasta", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    build_constructs(args.fasta, args.audit, args.output_fasta, args.manifest)


if __name__ == "__main__":
    main()
