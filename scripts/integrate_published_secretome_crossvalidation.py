#!/usr/bin/env python3
"""Map published PWN secretome sequences to BXYJv5 and integrate support."""

from __future__ import annotations

import argparse
import csv
import logging
from collections import Counter, defaultdict
from pathlib import Path


MAPPING_FIELDS = [
    "external_protein_id",
    "mapping_status",
    "bxyjv5_protein_id",
    "accepted_top_hit_count",
    "best_pident",
    "best_query_coverage",
    "best_subject_coverage",
    "best_evalue",
    "best_bitscore",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_alignments(path: Path) -> dict[str, list[dict[str, str]]]:
    by_query: dict[str, list[dict[str, str]]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 10:
                raise ValueError(f"DIAMOND row {line_number} must contain 10 fields")
            qseqid, sseqid, pident, length, qlen, slen, evalue, bitscore, qcov, scov = fields
            by_query[qseqid].append(
                {
                    "qseqid": qseqid,
                    "sseqid": sseqid,
                    "pident": pident,
                    "length": length,
                    "qlen": qlen,
                    "slen": slen,
                    "evalue": evalue,
                    "bitscore": bitscore,
                    "qcov": qcov,
                    "scov": scov,
                }
            )
    return by_query


def classify_mappings(
    manifest_path: Path,
    alignments_path: Path,
    min_identity: float = 90.0,
    min_query_coverage: float = 80.0,
    min_subject_coverage: float = 80.0,
    max_evalue: float = 1e-20,
) -> list[dict[str, str]]:
    manifest = read_tsv(manifest_path)
    alignments = parse_alignments(alignments_path)
    rows: list[dict[str, str]] = []
    for source in manifest:
        identifier = source["external_protein_id"]
        hits = alignments.get(identifier, [])
        if source.get("sequence_status") == "missing":
            status = "source_sequence_missing"
            accepted: list[dict[str, str]] = []
            best = None
        elif not hits:
            status = "no_alignment"
            accepted = []
            best = None
        else:
            best = max(hits, key=lambda hit: float(hit["bitscore"]))
            accepted = [
                hit
                for hit in hits
                if float(hit["pident"]) >= min_identity
                and float(hit["qcov"]) >= min_query_coverage
                and float(hit["scov"]) >= min_subject_coverage
                and float(hit["evalue"]) <= max_evalue
            ]
            if not accepted:
                status = "below_threshold"
            else:
                top_score = max(float(hit["bitscore"]) for hit in accepted)
                accepted = [hit for hit in accepted if abs(float(hit["bitscore"]) - top_score) <= 1e-9]
                best = accepted[0]
                if len({hit["sseqid"] for hit in accepted}) > 1:
                    status = "ambiguous_top_hit"
                elif all(float(best[key]) == 100.0 for key in ("pident", "qcov", "scov")):
                    status = "exact_sequence"
                else:
                    status = "high_confidence_sequence"
        accepted_target = (
            best["sseqid"]
            if best is not None and status in {"exact_sequence", "high_confidence_sequence"}
            else "NA"
        )
        rows.append(
            {
                "external_protein_id": identifier,
                "mapping_status": status,
                "bxyjv5_protein_id": accepted_target,
                "accepted_top_hit_count": str(len({hit["sseqid"] for hit in accepted})),
                "best_pident": best["pident"] if best else "NA",
                "best_query_coverage": best["qcov"] if best else "NA",
                "best_subject_coverage": best["scov"] if best else "NA",
                "best_evalue": best["evalue"] if best else "NA",
                "best_bitscore": best["bitscore"] if best else "NA",
            }
        )
    return rows


def integrate_support(
    source_records_path: Path,
    mappings_path: Path,
    functional_path: Path,
    source_mappings_output_path: Path,
    output_path: Path,
    report_path: Path,
) -> list[dict[str, str]]:
    source_records = read_tsv(source_records_path)
    mappings = {row["external_protein_id"]: row for row in read_tsv(mappings_path)}
    functional_rows = read_tsv(functional_path)
    accepted_statuses = {"exact_sequence", "high_confidence_sequence"}

    source_mappings: list[dict[str, str]] = []
    sources_by_protein: dict[str, set[str]] = defaultdict(set)
    external_ids_by_protein: dict[str, set[str]] = defaultdict(set)
    for source in source_records:
        mapping = mappings.get(
            source["external_protein_id"],
            {
                "mapping_status": "not_evaluated",
                "bxyjv5_protein_id": "NA",
            },
        )
        joined = {**source, **mapping}
        source_mappings.append(joined)
        if mapping["mapping_status"] in accepted_statuses:
            protein_id = mapping["bxyjv5_protein_id"]
            sources_by_protein[protein_id].add(source["source_id"])
            external_ids_by_protein[protein_id].add(source["external_protein_id"])

    source_mapping_fields = list(source_records[0]) + [
        field for field in MAPPING_FIELDS if field != "external_protein_id"
    ]
    write_tsv(source_mappings_output_path, source_mappings, source_mapping_fields)

    output_rows: list[dict[str, str]] = []
    for row in functional_rows:
        protein_id = row["protein_id"]
        published = bool(sources_by_protein[protein_id])
        channel_count = int(row.get("functional_support_count", "0")) + int(published)
        row = dict(row)
        row.update(
            {
                "published_secretome_support": "yes" if published else "no",
                "published_secretome_source_count": str(len(sources_by_protein[protein_id])),
                "published_secretome_sources": ";".join(sorted(sources_by_protein[protein_id])) or "NA",
                "published_secretome_external_ids": ";".join(
                    sorted(external_ids_by_protein[protein_id])
                )
                or "NA",
                "functional_support_count": str(channel_count),
                "evidence_tier": (
                    "multi_source_functional_support"
                    if channel_count >= 2
                    else "single_source_functional_support"
                    if channel_count == 1
                    else "sequence_defined_only"
                ),
                "claim_ceiling": (
                    "functionally or externally supported candidate secreted-soluble protein"
                    if channel_count
                    else "candidate secreted-soluble protein"
                ),
            }
        )
        output_rows.append(row)
    write_tsv(output_path, output_rows, list(output_rows[0]))

    status_counts = Counter(row.get("mapping_status", "not_evaluated") for row in mappings.values())
    source_counts = Counter(row["source_id"] for row in source_records)
    functional_ids = {row["protein_id"] for row in functional_rows}
    standard_targets_by_source: dict[str, set[str]] = defaultdict(set)
    for row in source_mappings:
        if (
            row["mapping_status"] in accepted_statuses
            and row["bxyjv5_protein_id"] in functional_ids
        ):
            standard_targets_by_source[row["source_id"]].add(row["bxyjv5_protein_id"])
    supported = sum(row["published_secretome_support"] == "yes" for row in output_rows)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        "# Published PWN secretome cross-validation\n\n"
        + "## Source records\n\n"
        + "\n".join(
            f"- {source}: {count} records; "
            f"{len(standard_targets_by_source[source])} mapped standard candidates"
            for source, count in sorted(source_counts.items())
        )
        + "\n\n## Sequence mapping\n\n"
        + "\n".join(f"- {status}: {count}" for status, count in sorted(status_counts.items()))
        + f"\n\n- Standard candidates supported by a published secretome: {supported}\n"
        + f"- Standard candidates without mapped published-secretome support: {len(output_rows) - supported}\n\n"
        "Mappings required a unique best BXYJv5 hit with amino-acid identity >= 90%, "
        "query coverage >= 80%, subject coverage >= 80%, and E-value <= 1e-20. "
        "Tied best hits, missing legacy sequences, and sub-threshold alignments were retained "
        "but did not count as external support. Published proteomic detection supports prior "
        "secretion evidence only; it does not validate effector activity, infection-stage delivery, "
        "direct host targeting, or causality.\n",
        encoding="utf-8",
    )
    logging.info("Integrated published-secretome support for %d standard candidates", supported)
    return output_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--alignments", type=Path, required=True)
    parser.add_argument("--source-records", type=Path, required=True)
    parser.add_argument("--functional", type=Path, required=True)
    parser.add_argument("--mapping-output", type=Path, required=True)
    parser.add_argument("--source-mappings-output", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--min-identity", type=float, default=90.0)
    parser.add_argument("--min-query-coverage", type=float, default=80.0)
    parser.add_argument("--min-subject-coverage", type=float, default=80.0)
    parser.add_argument("--max-evalue", type=float, default=1e-20)
    parser.add_argument("--log", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    handlers: list[logging.Handler] = [logging.StreamHandler()]
    if args.log:
        args.log.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(args.log, encoding="utf-8"))
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", handlers=handlers)
    mappings = classify_mappings(
        args.manifest,
        args.alignments,
        args.min_identity,
        args.min_query_coverage,
        args.min_subject_coverage,
        args.max_evalue,
    )
    write_tsv(args.mapping_output, mappings, MAPPING_FIELDS)
    integrate_support(
        args.source_records,
        args.mapping_output,
        args.functional,
        args.source_mappings_output,
        args.output,
        args.report,
    )


if __name__ == "__main__":
    main()
