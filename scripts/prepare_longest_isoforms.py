#!/usr/bin/env python3
"""Select the longest protein isoform per gene from a FASTA file."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


ISOFORM_SUFFIX = re.compile(r"^(.+)\.\d+$")


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


def wrap_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[index : index + width] for index in range(0, len(sequence), width))


def gene_id_from_header(header: str) -> str:
    token = header.split()[0]
    protein_id = token.split("|", 1)[1] if "|" in token else token
    match = ISOFORM_SUFFIX.match(protein_id)
    return match.group(1) if match else protein_id


def select_longest_isoforms(
    input_path: str | Path,
    output_path: str | Path,
    summary_path: str | Path | None = None,
    log_path: str | Path | None = None,
) -> dict[str, int]:
    input_path = Path(input_path)
    output_path = Path(output_path)
    selected: dict[str, tuple[str, str]] = {}
    isoform_counts: dict[str, int] = {}
    input_records = 0

    for header, sequence in read_fasta(input_path):
        input_records += 1
        gene_id = gene_id_from_header(header)
        isoform_counts[gene_id] = isoform_counts.get(gene_id, 0) + 1
        current = selected.get(gene_id)
        if current is None or len(sequence) > len(current[1]):
            selected[gene_id] = (header, sequence)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        for header, sequence in selected.values():
            handle.write(f">{header}\n{wrap_sequence(sequence)}\n")

    summary = {
        "input_records": input_records,
        "genes_written": len(selected),
        "isoforms_removed": input_records - len(selected),
        "multi_isoform_genes": sum(1 for count in isoform_counts.values() if count > 1),
    }

    if summary_path is not None:
        summary_path = Path(summary_path)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        fields = ["input_records", "genes_written", "isoforms_removed", "multi_isoform_genes"]
        with summary_path.open("w", encoding="utf-8") as handle:
            handle.write("\t".join(fields) + "\n")
            handle.write("\t".join(str(summary[field]) for field in fields) + "\n")
    if log_path is not None:
        log_path = Path(log_path)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(
            f"Wrote {summary['genes_written']} longest isoforms from {summary['input_records']} input records.\n",
            encoding="utf-8",
        )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary")
    parser.add_argument("--log")
    parser.add_argument("--mapping-table")
    args = parser.parse_args()
    if args.mapping_table:
        raise SystemExit("Mapping-table isoform selection is not implemented yet; use terminal numeric suffixes.")
    select_longest_isoforms(args.input, args.output, args.summary, args.log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
