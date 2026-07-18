#!/usr/bin/env python3
"""Build a gene-level map for the NCBI GIIE Pinus strobus TSA reference."""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import re
from collections import Counter
from pathlib import Path
from typing import TextIO


FIELDS = ["transcript_id", "gene_id"]
TRINITY_GENE = re.compile(r"\b(?P<gene_id>c\d+_g\d+)(?:_i\d+)?(?:,|\s|$)")


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def build_giie_tx2gene(fasta_path: Path, output_path: Path) -> list[dict[str, str]]:
    """Map GIIE TSA accession IDs to their original Trinity gene clusters."""
    rows: list[dict[str, str]] = []
    seen_transcripts: set[str] = set()
    with open_text(fasta_path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            transcript_id = header.split(maxsplit=1)[0] if header else ""
            if not transcript_id:
                raise ValueError(f"Empty FASTA identifier at line {line_number}")
            if transcript_id in seen_transcripts:
                raise ValueError(f"Duplicate transcript ID: {transcript_id}")
            match = TRINITY_GENE.search(header)
            if match is None:
                raise ValueError(
                    f"Missing Trinity gene cluster in FASTA header at line {line_number}: {header}"
                )
            seen_transcripts.add(transcript_id)
            rows.append({"transcript_id": transcript_id, "gene_id": match.group("gene_id")})
    if not rows:
        raise ValueError(f"No FASTA identifiers found in {fasta_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Wrote %d transcript-to-gene rows across %d Trinity gene clusters", len(rows), len({row['gene_id'] for row in rows}))
    return rows


def write_markdown(path: Path, rows: list[dict[str, str]], fasta_path: Path, output_path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    gene_counts = Counter(row["gene_id"] for row in rows)
    lines = [
        "# P. strobus GIIE TSA Transcript-to-Gene Map",
        "",
        f"- Source FASTA: `{fasta_path}`",
        f"- Output map: `{output_path}`",
        f"- Transcript accessions: {len(rows)}",
        f"- Trinity gene clusters: {len(gene_counts)}",
        f"- Gene clusters with multiple isoforms: {sum(count > 1 for count in gene_counts.values())}",
        "",
        "NCBI TSA accession identifiers are mapped to the original Trinity `c<number>_g<number>` cluster found in each FASTA header. This supports transcript-to-gene aggregation but is not a genome-coordinate annotation.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = build_parser().parse_args(argv)
    rows = build_giie_tx2gene(args.fasta, args.output)
    if args.markdown:
        write_markdown(args.markdown, rows, args.fasta, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
