#!/usr/bin/env python3
"""Normalize whitespace in a database FASTA while preserving identifiers."""

from __future__ import annotations

import argparse
import csv
import hashlib
import logging
import re
from pathlib import Path


VALID_SEQUENCE = re.compile(r"^[A-Za-z*?.-]+$")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalize_fasta(
    source_path: Path,
    output_path: Path,
    summary_path: Path,
    input_encoding: str = "utf-8",
) -> dict[str, int | str]:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    records = 0
    sequence_lines = 0
    changed_lines = 0
    removed = 0
    saw_header = False

    with source_path.open(encoding=input_encoding) as source, output_path.open(
        "w", encoding="utf-8", newline="\n"
    ) as output:
        for line_number, raw_line in enumerate(source, start=1):
            line = raw_line.rstrip("\r\n")
            if line.startswith(">"):
                records += 1
                saw_header = True
                output.write(line + "\n")
                continue
            if not line.strip():
                continue
            if not saw_header:
                raise ValueError(f"Sequence encountered before FASTA header on line {line_number}")
            sequence_lines += 1
            normalized = "".join(line.split())
            whitespace_removed = len(line) - len(normalized)
            if whitespace_removed:
                changed_lines += 1
                removed += whitespace_removed
            if not VALID_SEQUENCE.fullmatch(normalized):
                raise ValueError(f"Invalid sequence character on line {line_number}")
            output.write(normalized + "\n")

    if not records:
        raise ValueError("Input FASTA contains no records")
    summary: dict[str, int | str] = {
        "records": records,
        "sequence_lines": sequence_lines,
        "sequence_lines_changed": changed_lines,
        "whitespace_characters_removed": removed,
        "input_encoding": input_encoding,
        "output_encoding": "utf-8",
        "input_sha256": sha256_file(source_path),
        "output_sha256": sha256_file(output_path),
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary), delimiter="\t")
        writer.writeheader()
        writer.writerow(summary)
    logging.info(
        "Normalized %d records; removed %d whitespace characters from %d sequence lines",
        records,
        removed,
        changed_lines,
    )
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--input-encoding", default="utf-8")
    parser.add_argument("--log", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    handlers: list[logging.Handler] = [logging.StreamHandler()]
    if args.log:
        args.log.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(args.log, encoding="utf-8"))
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=handlers,
    )
    normalize_fasta(args.input, args.output, args.summary, args.input_encoding)


if __name__ == "__main__":
    main()
