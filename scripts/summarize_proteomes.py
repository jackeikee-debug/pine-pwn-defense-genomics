#!/usr/bin/env python3
"""Summarize processed proteome FASTA files."""

from __future__ import annotations

import argparse
from pathlib import Path


def read_fasta(path: Path):
    header = None
    sequence_parts = []
    with path.open("r", encoding="utf-8") as handle:
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


def median(values: list[int]) -> int:
    sorted_values = sorted(values)
    midpoint = len(sorted_values) // 2
    if len(sorted_values) % 2:
        return sorted_values[midpoint]
    return int(round((sorted_values[midpoint - 1] + sorted_values[midpoint]) / 2))


def protein_n50(lengths: list[int]) -> int:
    threshold = sum(lengths) / 2
    running_total = 0
    for length in sorted(lengths, reverse=True):
        running_total += length
        if running_total >= threshold:
            return length
    return 0


def summarize_fasta(path: Path) -> dict[str, int | str]:
    lengths: list[int] = []
    seen_ids: set[str] = set()
    duplicated_ids = 0
    missing_sequences = 0
    species_id = path.stem

    for header, sequence in read_fasta(path):
        gene_id = header.split()[0]
        if "|" in gene_id:
            species_id = gene_id.split("|", 1)[0]
        if gene_id in seen_ids:
            duplicated_ids += 1
        seen_ids.add(gene_id)
        clean_sequence = "".join(sequence.split())
        if not clean_sequence:
            missing_sequences += 1
            continue
        lengths.append(len(clean_sequence))

    return {
        "species_id": species_id,
        "protein_count": len(lengths),
        "median_length": median(lengths) if lengths else 0,
        "protein_n50": protein_n50(lengths) if lengths else 0,
        "missing_sequences": missing_sequences,
        "duplicated_ids": duplicated_ids,
    }


def write_summary_table(rows: list[dict[str, int | str]], output_path: Path) -> None:
    fields = [
        "species_id",
        "protein_count",
        "median_length",
        "protein_n50",
        "missing_sequences",
        "duplicated_ids",
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(fields) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[field]) for field in fields) + "\n")


def summarize_fastas(
    fasta_paths: list[str | Path],
    output_path: str | Path,
    log_path: str | Path | None = None,
):
    rows = [summarize_fasta(Path(path)) for path in sorted(fasta_paths)]
    rows.sort(key=lambda row: str(row["species_id"]))
    write_summary_table(rows, Path(output_path))
    if log_path is not None:
        Path(log_path).parent.mkdir(parents=True, exist_ok=True)
        Path(log_path).write_text(f"Summarized {len(rows)} proteomes.\n", encoding="utf-8")
    return rows


def summarize_directory(input_dir: str | Path, output_path: str | Path, log_path: str | Path | None = None):
    input_dir = Path(input_dir)
    return summarize_fastas(list(input_dir.glob("*.faa")), output_path, log_path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--input-dir")
    input_group.add_argument("--inputs", nargs="+")
    parser.add_argument("--output", required=True)
    parser.add_argument("--log")
    args = parser.parse_args()
    if args.inputs:
        summarize_fastas(args.inputs, args.output, args.log)
    else:
        summarize_directory(args.input_dir, args.output, args.log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
