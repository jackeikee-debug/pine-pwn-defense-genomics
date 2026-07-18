#!/usr/bin/env python3
"""Build an identity tx2gene map for P. massoniana annotated CDS features."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


FIELDS = ["transcript_id", "gene_id"]


def build_identity_tx2gene(fasta_path: Path, output_path: Path) -> list[dict[str, str]]:
    """Map each distinct annotated CDS feature ID to itself for tximport."""
    rows: list[dict[str, str]] = []
    seen: set[str] = set()
    with fasta_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.startswith(">"):
                continue
            feature_id = line[1:].strip().split()[0]
            if not feature_id:
                raise ValueError(f"Empty FASTA identifier at line {line_number}")
            if feature_id in seen:
                raise ValueError(f"Duplicate CDS feature ID: {feature_id}")
            seen.add(feature_id)
            rows.append({"transcript_id": feature_id, "gene_id": feature_id})
    if not rows:
        raise ValueError(f"No FASTA identifiers found in {fasta_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Wrote identity tx2gene mapping for %d annotated CDS features", len(rows))
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = parse_args()
    build_identity_tx2gene(args.fasta, args.output)


if __name__ == "__main__":
    main()
