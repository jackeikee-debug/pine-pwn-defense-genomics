#!/usr/bin/env python3
"""Build a manuscript-facing proteome source and quality audit table."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


FIELDS = [
    "species_id", "species_name", "region", "panel_role", "source", "accession", "version",
    "protein_fasta", "protein_count", "median_length", "protein_n50", "missing_sequences",
    "duplicated_ids", "busco_lineage", "busco_complete_single_copy_pct",
    "busco_complete_duplicated_pct", "busco_fragmented_pct", "busco_missing_pct",
    "busco_version", "busco_status", "quality_claim_ceiling",
]
BUSCO_FIELDS = [
    "species_id", "lineage", "complete_single_copy_pct", "complete_duplicated_pct",
    "fragmented_pct", "missing_pct", "busco_version",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_busco_short_summary(species_id: str, path: Path) -> dict[str, str]:
    text = path.read_text(encoding="utf-8")
    version_match = re.search(r"BUSCO version is:\s*([^\s]+)", text)
    lineage_match = re.search(r"lineage dataset is:\s*([^\s(]+)", text, flags=re.IGNORECASE)
    score_match = re.search(
        r"C:[\d.]+%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:\d+",
        text,
    )
    if not version_match or not lineage_match or not score_match:
        raise ValueError(f"Could not parse BUSCO short summary: {path}")
    return {
        "species_id": species_id,
        "lineage": lineage_match.group(1),
        "complete_single_copy_pct": score_match.group(1),
        "complete_duplicated_pct": score_match.group(2),
        "fragmented_pct": score_match.group(3),
        "missing_pct": score_match.group(4),
        "busco_version": version_match.group(1),
    }


def write_busco_table(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=BUSCO_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: row["species_id"]))


def build_proteome_quality_audit(
    metadata_path: Path,
    summary_path: Path,
    busco_path: Path | None,
    output_path: Path,
    report_path: Path,
) -> tuple[list[dict[str, str]], dict[str, int]]:
    metadata = {row["species_id"]: row for row in read_tsv(metadata_path)}
    summaries = read_tsv(summary_path)
    busco = {row["species_id"]: row for row in read_tsv(busco_path)} if busco_path else {}
    rows: list[dict[str, str]] = []

    for summary in summaries:
        species_id = summary["species_id"]
        if species_id not in metadata:
            raise ValueError(f"No species metadata for {species_id}")
        meta = metadata[species_id]
        busco_row = busco.get(species_id)
        busco_status = "complete" if busco_row else "not_run"
        rows.append({
            "species_id": species_id,
            "species_name": meta.get("species_name", ""),
            "region": meta.get("region", ""),
            "panel_role": meta.get("group", ""),
            "source": meta.get("source", ""),
            "accession": meta.get("accession", ""),
            "version": meta.get("version", ""),
            "protein_fasta": meta.get("protein_fasta", ""),
            "protein_count": summary.get("protein_count", ""),
            "median_length": summary.get("median_length", ""),
            "protein_n50": summary.get("protein_n50", ""),
            "missing_sequences": summary.get("missing_sequences", ""),
            "duplicated_ids": summary.get("duplicated_ids", ""),
            "busco_lineage": "" if busco_row is None else busco_row.get("lineage", ""),
            "busco_complete_single_copy_pct": "" if busco_row is None else busco_row.get("complete_single_copy_pct", ""),
            "busco_complete_duplicated_pct": "" if busco_row is None else busco_row.get("complete_duplicated_pct", ""),
            "busco_fragmented_pct": "" if busco_row is None else busco_row.get("fragmented_pct", ""),
            "busco_missing_pct": "" if busco_row is None else busco_row.get("missing_pct", ""),
            "busco_version": "" if busco_row is None else busco_row.get("busco_version", ""),
            "busco_status": busco_status,
            "quality_claim_ceiling": "sequence_statistics_plus_busco" if busco_row else "descriptive_sequence_statistics_only",
        })

    rows.sort(key=lambda row: row["species_id"])
    counts = {
        "panel_species": len(rows),
        "busco_complete": sum(row["busco_status"] == "complete" for row in rows),
        "busco_not_run": sum(row["busco_status"] == "not_run" for row in rows),
    }
    write_tsv(output_path, rows)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        "# Host proteome quality audit\n\n"
        f"- Proteome sets in the orthology panel: {counts['panel_species']}.\n"
        f"- Sets with parsed BUSCO summaries: {counts['busco_complete']}.\n"
        f"- Sets without BUSCO summaries: {counts['busco_not_run']}.\n\n"
        + (
            "BUSCO has not been run for the current panel. Protein count, median length, protein N50, missing sequences, "
            "and duplicate identifiers are descriptive sequence-set statistics only; they do not demonstrate equivalent completeness.\n"
            if counts["busco_not_run"] else
            "BUSCO summaries were parsed for every panel proteome. Source and annotation-version differences still require interpretation.\n"
        ),
        encoding="utf-8",
    )
    return rows, counts


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--busco", type=Path)
    parser.add_argument(
        "--busco-summary", action="append", default=[], metavar="SPECIES_ID=PATH",
        help="Parse a BUSCO short summary; may be repeated.",
    )
    parser.add_argument("--busco-table-output", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    busco_path = args.busco
    if args.busco_summary:
        if args.busco is not None or args.busco_table_output is None:
            raise ValueError("--busco-summary requires --busco-table-output and cannot be combined with --busco")
        parsed = []
        for specification in args.busco_summary:
            if "=" not in specification:
                raise ValueError(f"Expected SPECIES_ID=PATH, found: {specification}")
            species_id, summary_path = specification.split("=", 1)
            parsed.append(parse_busco_short_summary(species_id, Path(summary_path)))
        write_busco_table(args.busco_table_output, parsed)
        busco_path = args.busco_table_output
    _, counts = build_proteome_quality_audit(
        args.metadata, args.summary, busco_path, args.output, args.report
    )
    print(f"Audited {counts['panel_species']} proteomes; BUSCO available for {counts['busco_complete']}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
