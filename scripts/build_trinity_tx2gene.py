#!/usr/bin/env python3
"""Build a transcript-to-gene map from Trinity FASTA identifiers."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter
from pathlib import Path


FIELDS = ["transcript_id", "gene_id"]
ISOFORM_SUFFIX = re.compile(r"^(?P<gene_id>.+)_i\d+$")


def trinity_gene_id(transcript_id: str) -> str:
    match = ISOFORM_SUFFIX.match(transcript_id)
    return match.group("gene_id") if match else transcript_id


def read_fasta_ids(path: Path) -> list[str]:
    identifiers = []
    seen = set()
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.startswith(">"):
                continue
            transcript_id = line[1:].strip().split()[0]
            if not transcript_id:
                raise ValueError(f"Empty FASTA identifier at line {line_number}")
            if transcript_id in seen:
                raise ValueError(f"Duplicate transcript ID: {transcript_id}")
            seen.add(transcript_id)
            identifiers.append(transcript_id)
    if not identifiers:
        raise ValueError(f"No FASTA identifiers found in {path}")
    return identifiers


def build_trinity_tx2gene(
    fasta_path: Path,
    output_path: Path,
    markdown_path: Path | None = None,
) -> list[dict[str, str]]:
    rows = [
        {"transcript_id": transcript_id, "gene_id": trinity_gene_id(transcript_id)}
        for transcript_id in read_fasta_ids(fasta_path)
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if markdown_path is not None:
        write_markdown(markdown_path, rows, fasta_path, output_path)
    return rows


def write_markdown(
    path: Path,
    rows: list[dict[str, str]],
    fasta_path: Path,
    output_path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    gene_counts = Counter(row["gene_id"] for row in rows)
    multi_isoform_genes = sum(count > 1 for count in gene_counts.values())
    lines = [
        "# P. densiflora Trinity Transcript-to-Gene Map",
        "",
        f"- Source FASTA: `{fasta_path}`",
        f"- Output map: `{output_path}`",
        f"- Transcript IDs: {len(rows)}",
        f"- Gene IDs: {len(gene_counts)}",
        f"- Genes with multiple isoforms: {multi_isoform_genes}",
        "",
        "Gene identifiers are derived only by removing a terminal Trinity `_i<number>` isoform suffix. The map supports Salmon gene-level aggregation and tximport; it is not a genome annotation crosswalk.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = build_trinity_tx2gene(args.fasta, args.output, args.markdown)
    gene_count = len({row["gene_id"] for row in rows})
    print(f"Trinity tx2gene map written: {len(rows)} transcripts; {gene_count} genes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
