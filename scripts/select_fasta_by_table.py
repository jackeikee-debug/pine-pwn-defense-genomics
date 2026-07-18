#!/usr/bin/env python3
"""Select FASTA records whose IDs appear in a TSV table."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


MANIFEST_FIELDS = ["record_id", "source_fasta", "selected_fasta"]


def read_fasta(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    header: str | None = None
    sequence_parts: list[str] = []
    with opener(path, "rt", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_parts)
                header = line[1:]
                sequence_parts = []
            else:
                sequence_parts.append(line)
    if header is not None:
        yield header, "".join(sequence_parts)


def read_selected_ids(
    table_path: Path,
    id_column: str,
    where_column: str | None,
    where_value: str | None,
) -> set[str]:
    selected: set[str] = set()
    with table_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if where_column and row.get(where_column) != where_value:
                continue
            selected.add(row[id_column])
    return selected


def select_fasta_by_table(
    fasta_path: Path,
    table_path: Path,
    id_column: str,
    output_fasta: Path,
    manifest_output: Path,
    where_column: str | None = None,
    where_value: str | None = None,
) -> list[dict[str, str]]:
    selected_ids = read_selected_ids(table_path, id_column, where_column, where_value)
    rows: list[dict[str, str]] = []
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with output_fasta.open("w", encoding="utf-8") as handle:
        for header, sequence in read_fasta(fasta_path):
            record_id = header.split()[0]
            if record_id not in selected_ids:
                continue
            handle.write(f">{record_id}\n{wrap_sequence(sequence)}\n")
            rows.append(
                {
                    "record_id": record_id,
                    "source_fasta": str(fasta_path),
                    "selected_fasta": str(output_fasta),
                }
            )
    write_tsv(manifest_output, MANIFEST_FIELDS, rows)
    return rows


def wrap_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[index : index + width] for index in range(0, len(sequence), width))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--table", type=Path, required=True)
    parser.add_argument("--id-column", required=True)
    parser.add_argument("--output-fasta", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--where-column")
    parser.add_argument("--where-value")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = select_fasta_by_table(
        fasta_path=args.fasta,
        table_path=args.table,
        id_column=args.id_column,
        output_fasta=args.output_fasta,
        manifest_output=args.manifest,
        where_column=args.where_column,
        where_value=args.where_value,
    )
    print(f"Selected FASTA records: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
