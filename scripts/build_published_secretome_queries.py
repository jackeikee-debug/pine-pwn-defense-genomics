#!/usr/bin/env python3
"""Extract published PWN secretome identifiers and their legacy sequences."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import logging
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


RECORD_FIELDS = [
    "source_id",
    "citation_key",
    "doi",
    "pmcid",
    "evidence_type",
    "external_protein_id",
    "source_file",
    "source_sheet",
    "source_row_number",
]
MANIFEST_FIELDS = [
    "external_protein_id",
    "source_record_count",
    "source_ids",
    "sequence_status",
    "sequence_source",
    "sequence_length",
    "sequence_sha256",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(
        "r", encoding="utf-8"
    )


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier: str | None = None
    chunks: list[str] = []
    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    if identifier in records:
                        raise ValueError(f"Duplicate FASTA identifier in {path}: {identifier}")
                    records[identifier] = "".join(chunks)
                identifier = line[1:].split()[0]
                chunks = []
            elif identifier is None:
                raise ValueError(f"Sequence before first FASTA header in {path}")
            else:
                chunks.append(line)
    if identifier is not None:
        if identifier in records:
            raise ValueError(f"Duplicate FASTA identifier in {path}: {identifier}")
        records[identifier] = "".join(chunks)
    return records


def resolve_path(raw_path: str, config_path: Path) -> Path:
    path = Path(raw_path)
    if path.is_absolute():
        return path
    cwd_path = Path.cwd() / path
    return cwd_path if cwd_path.exists() else config_path.parent / path


def extract_source_records(config_path: Path) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    for source in read_tsv(config_path):
        source_path = resolve_path(source["source_file"], config_path)
        header_row = int(source["header_row"])
        frame = pd.read_excel(source_path, sheet_name=source["sheet"], header=header_row)
        if source["id_column"] not in frame.columns:
            raise ValueError(
                f"Missing ID column {source['id_column']!r} in {source_path} sheet {source['sheet']}"
            )
        for index, value in frame[source["id_column"]].items():
            if pd.isna(value):
                continue
            identifier = str(value).strip()
            if not identifier or identifier.lower() == "nan":
                continue
            records.append(
                {
                    "source_id": source["source_id"],
                    "citation_key": source["citation_key"],
                    "doi": source["doi"],
                    "pmcid": source["pmcid"],
                    "evidence_type": source["evidence_type"],
                    "external_protein_id": identifier,
                    "source_file": source["source_file"],
                    "source_sheet": source["sheet"],
                    "source_row_number": str(header_row + 2 + int(index)),
                }
            )
    return records


def build_queries(
    sources_path: Path,
    wbps_proteome_path: Path,
    legacy_secretome_fasta_path: Path,
    records_output_path: Path,
    query_fasta_path: Path,
    manifest_output_path: Path,
) -> list[dict[str, str]]:
    source_records = extract_source_records(sources_path)
    write_tsv(records_output_path, source_records, RECORD_FIELDS)

    wbps = read_fasta(wbps_proteome_path)
    legacy = read_fasta(legacy_secretome_fasta_path)
    conflicting = {
        identifier for identifier in wbps.keys() & legacy.keys() if wbps[identifier] != legacy[identifier]
    }
    if conflicting:
        raise ValueError(f"Conflicting legacy sequences: {', '.join(sorted(conflicting)[:10])}")

    record_counts = Counter(row["external_protein_id"] for row in source_records)
    sources_by_id: dict[str, set[str]] = defaultdict(set)
    for row in source_records:
        sources_by_id[row["external_protein_id"]].add(row["source_id"])

    manifest: list[dict[str, str]] = []
    available: dict[str, str] = {}
    for identifier in sorted(record_counts):
        if identifier in legacy:
            sequence = legacy[identifier]
            sequence_source = str(legacy_secretome_fasta_path)
        elif identifier in wbps:
            sequence = wbps[identifier]
            sequence_source = str(wbps_proteome_path)
        else:
            sequence = ""
            sequence_source = "NA"
        if sequence:
            available[identifier] = sequence
        manifest.append(
            {
                "external_protein_id": identifier,
                "source_record_count": str(record_counts[identifier]),
                "source_ids": ";".join(sorted(sources_by_id[identifier])),
                "sequence_status": "available" if sequence else "missing",
                "sequence_source": sequence_source,
                "sequence_length": str(len(sequence)) if sequence else "NA",
                "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest() if sequence else "NA",
            }
        )

    query_fasta_path.parent.mkdir(parents=True, exist_ok=True)
    with query_fasta_path.open("w", encoding="utf-8", newline="\n") as handle:
        for identifier in sorted(available):
            handle.write(f">{identifier}\n")
            sequence = available[identifier]
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")
    write_tsv(manifest_output_path, manifest, MANIFEST_FIELDS)
    logging.info(
        "Extracted %d source rows, %d unique identifiers, and %d available sequences",
        len(source_records),
        len(manifest),
        len(available),
    )
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--wbps-proteome", type=Path, required=True)
    parser.add_argument("--legacy-secretome-fasta", type=Path, required=True)
    parser.add_argument("--records-output", type=Path, required=True)
    parser.add_argument("--query-fasta", type=Path, required=True)
    parser.add_argument("--manifest-output", type=Path, required=True)
    parser.add_argument("--log", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    handlers: list[logging.Handler] = [logging.StreamHandler()]
    if args.log:
        args.log.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(args.log, encoding="utf-8"))
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", handlers=handlers)
    build_queries(
        args.sources,
        args.wbps_proteome,
        args.legacy_secretome_fasta,
        args.records_output,
        args.query_fasta,
        args.manifest_output,
    )


if __name__ == "__main__":
    main()
