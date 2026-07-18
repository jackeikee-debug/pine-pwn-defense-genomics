#!/usr/bin/env python3
"""Predict candidate secreted nematode proteins with transparent heuristics."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


FIELDS = [
    "protein_id",
    "annotation_text",
    "length",
    "cysteine_fraction",
    "n_terminal_hydrophobic_count",
    "signal_peptide_heuristic",
    "transmembrane_domain_heuristic_count",
    "small_secreted_candidate",
    "cysteine_rich_candidate",
    "secretome_candidate",
    "evidence",
]
HYDROPHOBIC = set("AILMFWV")


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


def protein_id_from_header(header: str) -> str:
    return header.split()[0]


def annotation_from_header(header: str) -> str:
    parts = header.split(maxsplit=1)
    return parts[1] if len(parts) > 1 else ""


def has_signal_peptide_heuristic(sequence: str) -> bool:
    if len(sequence) < 25 or not sequence.startswith("M"):
        return False
    n_region = sequence[1:20]
    hydrophobic_count = sum(1 for amino_acid in n_region if amino_acid in HYDROPHOBIC)
    polar_tail = sequence[15:30]
    return hydrophobic_count >= 8 and any(amino_acid in "ASTCG" for amino_acid in polar_tail)


def count_tm_domains_after_signal(sequence: str, window: int = 19, min_hydrophobic: int = 15) -> int:
    count = 0
    index = 30
    while index + window <= len(sequence):
        segment = sequence[index : index + window]
        if sum(1 for amino_acid in segment if amino_acid in HYDROPHOBIC) >= min_hydrophobic:
            count += 1
            index += window
        else:
            index += 1
    return count


def predict_secretome(
    proteins_path: Path,
    output_path: Path,
    min_length: int,
    small_secreted_max_length: int,
    cysteine_rich_threshold: float,
    max_transmembrane_domains: int,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for header, raw_sequence in read_fasta(proteins_path):
        sequence = raw_sequence.upper().replace("*", "")
        length = len(sequence)
        if length < min_length:
            continue
        cysteine_fraction = sequence.count("C") / length if length else 0.0
        signal_peptide = has_signal_peptide_heuristic(sequence)
        tm_count = count_tm_domains_after_signal(sequence)
        secretome_candidate = signal_peptide and tm_count <= max_transmembrane_domains
        small_secreted = secretome_candidate and length <= small_secreted_max_length
        cysteine_rich = secretome_candidate and cysteine_fraction >= cysteine_rich_threshold
        n_terminal_hydrophobic_count = sum(1 for amino_acid in sequence[1:20] if amino_acid in HYDROPHOBIC)
        evidence = []
        if signal_peptide:
            evidence.append("n_terminal_signal_peptide_heuristic")
        if tm_count > max_transmembrane_domains:
            evidence.append("post_signal_transmembrane_filtered")
        if small_secreted:
            evidence.append("small_secreted")
        if cysteine_rich:
            evidence.append("cysteine_rich")
        rows.append(
            {
                "protein_id": protein_id_from_header(header),
                "annotation_text": annotation_from_header(header),
                "length": str(length),
                "cysteine_fraction": f"{cysteine_fraction:.4f}",
                "n_terminal_hydrophobic_count": str(n_terminal_hydrophobic_count),
                "signal_peptide_heuristic": yes_no(signal_peptide),
                "transmembrane_domain_heuristic_count": str(tm_count),
                "small_secreted_candidate": yes_no(small_secreted),
                "cysteine_rich_candidate": yes_no(cysteine_rich),
                "secretome_candidate": yes_no(secretome_candidate),
                "evidence": ";".join(evidence) if evidence else "none",
            }
        )

    write_tsv(output_path, FIELDS, rows)
    return rows


def yes_no(value: bool) -> str:
    return "yes" if value else "no"


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proteins", type=Path, required=True)
    parser.add_argument("--signalp")
    parser.add_argument("--tm-domains")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--min-length", type=int, default=30)
    parser.add_argument("--small-secreted-max-length", type=int, default=300)
    parser.add_argument("--cysteine-rich-threshold", type=float, default=0.03)
    parser.add_argument("--max-transmembrane-domains", type=int, default=0)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = predict_secretome(
        proteins_path=args.proteins,
        output_path=args.output,
        min_length=args.min_length,
        small_secreted_max_length=args.small_secreted_max_length,
        cysteine_rich_threshold=args.cysteine_rich_threshold,
        max_transmembrane_domains=args.max_transmembrane_domains,
    )
    secreted = sum(1 for row in rows if row["secretome_candidate"] == "yes")
    print(f"Proteins evaluated: {len(rows)}")
    print(f"Secretome candidates written: {secreted}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
