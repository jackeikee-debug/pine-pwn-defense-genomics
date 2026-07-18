#!/usr/bin/env python3
"""Clean protein FASTA records and standardize headers to species_id|gene_id."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path


def read_fasta(path: Path):
    header = None
    sequence_parts = []
    open_fn = gzip.open if path.suffix == ".gz" else open
    with open_fn(path, "rt", encoding="utf-8") as handle:
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


def extract_gene_id(header: str) -> str:
    token = header.split()[0]
    if "|" in token:
        parts = [part for part in token.split("|") if part]
        if len(parts) >= 2:
            return parts[1]
    return token


def wrap_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[index : index + width] for index in range(0, len(sequence), width))


def clean_sequence(sequence: str) -> tuple[str, bool, int]:
    sequence = "".join(sequence.split()).upper()
    had_terminal_stop = sequence.endswith("*")
    if had_terminal_stop:
        sequence = sequence[:-1]
    internal_stop_count = sequence.count("*")
    sequence = sequence.replace("*", "X")
    return sequence, had_terminal_stop, internal_stop_count


def write_summary(path: Path, summary: dict[str, int | str]) -> None:
    fields = [
        "species_id",
        "input_records",
        "proteins_written",
        "short_sequence_count",
        "duplicated_id_count",
        "terminal_stop_count",
        "internal_stop_count",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(fields) + "\n")
        handle.write("\t".join(str(summary[field]) for field in fields) + "\n")


def clean_fasta(
    input_path: str | Path,
    output_path: str | Path,
    species_id: str,
    min_length: int = 30,
    summary_path: str | Path | None = None,
    log_path: str | Path | None = None,
) -> dict[str, int | str]:
    input_path = Path(input_path)
    output_path = Path(output_path)
    seen_gene_ids: set[str] = set()
    records: list[tuple[str, str]] = []
    summary: dict[str, int | str] = {
        "species_id": species_id,
        "input_records": 0,
        "proteins_written": 0,
        "short_sequence_count": 0,
        "duplicated_id_count": 0,
        "terminal_stop_count": 0,
        "internal_stop_count": 0,
    }

    for header, raw_sequence in read_fasta(input_path):
        summary["input_records"] += 1
        gene_id = extract_gene_id(header)
        if gene_id in seen_gene_ids:
            summary["duplicated_id_count"] += 1
            continue
        seen_gene_ids.add(gene_id)

        sequence, had_terminal_stop, internal_stop_count = clean_sequence(raw_sequence)
        summary["terminal_stop_count"] += int(had_terminal_stop)
        summary["internal_stop_count"] += internal_stop_count
        if len(sequence) < min_length:
            summary["short_sequence_count"] += 1
            continue
        records.append((f"{species_id}|{gene_id}", sequence))

    summary["proteins_written"] = len(records)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n{wrap_sequence(sequence)}\n")

    if summary_path is not None:
        write_summary(Path(summary_path), summary)
    if log_path is not None:
        Path(log_path).parent.mkdir(parents=True, exist_ok=True)
        Path(log_path).write_text(
            f"Wrote {summary['proteins_written']} proteins from {summary['input_records']} input records.\n",
            encoding="utf-8",
        )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--species-id", required=True)
    parser.add_argument("--min-length", type=int, default=30)
    parser.add_argument("--summary")
    parser.add_argument("--log")
    args = parser.parse_args()
    clean_fasta(args.input, args.output, args.species_id, args.min_length, args.summary, args.log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
