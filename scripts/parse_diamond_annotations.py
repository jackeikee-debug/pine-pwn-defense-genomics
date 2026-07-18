#!/usr/bin/env python3
"""Convert DIAMOND SwissProt tabular hits into project annotation tables."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


ANNOTATION_FIELDS = ["gene_id", "annotation_source", "annotation_text"]
HIT_FIELDS = [
    "gene_id",
    "subject_id",
    "pident",
    "alignment_length",
    "evalue",
    "bitscore",
    "annotation_text",
    "annotation_source",
]


def parse_diamond_annotations(
    diamond_path: Path,
    annotation_output: Path,
    hits_output: Path,
    min_pident: float,
    min_bitscore: float,
    source_name: str,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    best_hits: dict[str, dict[str, str]] = {}
    gene_order: list[str] = []

    if diamond_path.exists():
        with diamond_path.open("r", encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.rstrip("\n")
                if not line:
                    continue
                fields = line.split("\t")
                if len(fields) < 7:
                    continue
                gene_id, subject_id, pident, length, evalue, bitscore, annotation_text = fields[:7]
                if float(pident) < min_pident or float(bitscore) < min_bitscore:
                    continue
                hit = {
                    "gene_id": gene_id,
                    "subject_id": subject_id,
                    "pident": pident,
                    "alignment_length": length,
                    "evalue": evalue,
                    "bitscore": bitscore,
                    "annotation_text": annotation_text,
                    "annotation_source": source_name,
                }
                current = best_hits.get(gene_id)
                if current is None:
                    gene_order.append(gene_id)
                if current is None or float(bitscore) > float(current["bitscore"]):
                    best_hits[gene_id] = hit

    hit_rows = [best_hits[gene_id] for gene_id in gene_order]
    annotation_rows = [
        {
            "gene_id": row["gene_id"],
            "annotation_source": row["annotation_source"],
            "annotation_text": row["annotation_text"],
        }
        for row in hit_rows
    ]

    write_tsv(annotation_output, ANNOTATION_FIELDS, annotation_rows)
    write_tsv(hits_output, HIT_FIELDS, hit_rows)
    return annotation_rows, hit_rows


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diamond", type=Path, required=True)
    parser.add_argument("--annotation-output", type=Path, required=True)
    parser.add_argument("--hits-output", type=Path, required=True)
    parser.add_argument("--min-pident", type=float, default=35.0)
    parser.add_argument("--min-bitscore", type=float, default=80.0)
    parser.add_argument("--source-name", default="SwissProt_DIAMOND")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    annotation_rows, hit_rows = parse_diamond_annotations(
        diamond_path=args.diamond,
        annotation_output=args.annotation_output,
        hits_output=args.hits_output,
        min_pident=args.min_pident,
        min_bitscore=args.min_bitscore,
        source_name=args.source_name,
    )
    print(f"Protein function annotations written: {len(annotation_rows)}")
    print(f"SwissProt DIAMOND hits written: {len(hit_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
