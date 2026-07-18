#!/usr/bin/env python3
"""Extract the standard secreted-soluble BXYJv5 protein set deterministically."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import logging
from pathlib import Path
from typing import TextIO


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def passing_ids(path: Path) -> set[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"protein_id", "standard_secreted_soluble_status"}
    if not rows or not required <= set(rows[0]):
        raise ValueError(f"Audit table must contain {sorted(required)}")
    selected = [
        row["protein_id"]
        for row in rows
        if row["standard_secreted_soluble_status"] == "pass"
    ]
    if len(selected) != len(set(selected)):
        raise ValueError("Duplicate passing protein identifiers in audit table")
    return set(selected)


def read_fasta(path: Path) -> list[tuple[str, str, str]]:
    records: list[tuple[str, str, str]] = []
    seen: set[str] = set()
    header = ""
    sequence_parts: list[str] = []

    def finish_record() -> None:
        nonlocal header, sequence_parts
        if not header:
            return
        protein_id = header[1:].split()[0]
        if protein_id in seen:
            raise ValueError(f"Duplicate FASTA identifier: {protein_id}")
        sequence = "".join(sequence_parts).replace(" ", "").upper()
        if not sequence:
            raise ValueError(f"Empty FASTA sequence: {protein_id}")
        seen.add(protein_id)
        records.append((protein_id, header, sequence))

    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish_record()
                header = line
                sequence_parts = []
            elif not header:
                raise ValueError("FASTA sequence encountered before the first header")
            else:
                sequence_parts.append(line)
    finish_record()
    return records


def extract_standard_secretome(
    source_fasta: Path,
    audit_table: Path,
    output_fasta: Path,
    manifest_path: Path,
) -> dict[str, int]:
    selected = passing_ids(audit_table)
    records = read_fasta(source_fasta)
    source_ids = {record[0] for record in records}
    missing = sorted(selected - source_ids)
    if missing:
        preview = ", ".join(missing[:10])
        raise ValueError(f"Passing identifiers missing from source FASTA ({len(missing)}): {preview}")

    chosen = [record for record in records if record[0] in selected]
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with output_fasta.open("w", encoding="utf-8", newline="\n") as handle:
        for _, header, sequence in chosen:
            handle.write(f"{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")

    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_rows = [
        {
            "protein_id": protein_id,
            "sequence_length": str(len(sequence)),
            "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            "source_fasta": source_fasta.as_posix(),
            "selection_rule": "standard_secreted_soluble_status=pass",
        }
        for protein_id, _, sequence in chosen
    ]
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    summary = {
        "source_records": len(records),
        "selected_records": len(chosen),
        "missing_records": len(missing),
    }
    logging.info("Selected %d of %d source proteins", len(chosen), len(records))
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-fasta", type=Path, required=True)
    parser.add_argument("--audit-table", type=Path, required=True)
    parser.add_argument("--output-fasta", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--log", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    handlers: list[logging.Handler] = [logging.StreamHandler()]
    if args.log:
        args.log.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(args.log, encoding="utf-8"))
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", handlers=handlers)
    extract_standard_secretome(
        args.source_fasta,
        args.audit_table,
        args.output_fasta,
        args.manifest,
    )


if __name__ == "__main__":
    main()
