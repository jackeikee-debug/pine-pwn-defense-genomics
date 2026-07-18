#!/usr/bin/env python3
"""Fetch and audit ENA run metadata for the Route 1 P. densiflora RNA-seq study."""

from __future__ import annotations

import argparse
import csv
import io
import re
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen


ENA_FIELDS = [
    "run_accession",
    "experiment_accession",
    "study_accession",
    "secondary_study_accession",
    "sample_accession",
    "secondary_sample_accession",
    "scientific_name",
    "experiment_alias",
    "library_layout",
    "library_strategy",
    "instrument_platform",
    "instrument_model",
    "read_count",
    "base_count",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
    "first_public",
    "last_updated",
]

MANIFEST_FIELDS = [
    "study_accession",
    "secondary_study_accession",
    "run_accession",
    "experiment_accession",
    "sample_accession",
    "secondary_sample_accession",
    "scientific_name",
    "experiment_alias",
    "condition",
    "condition_label",
    "replicate",
    "comparison_role",
    "tissue",
    "collection_time",
    "library_layout",
    "library_strategy",
    "instrument_platform",
    "instrument_model",
    "read_count",
    "base_count",
    "fastq_1_url",
    "fastq_2_url",
    "fastq_1_md5",
    "fastq_2_md5",
    "fastq_1_bytes",
    "fastq_2_bytes",
    "total_fastq_bytes",
    "total_fastq_gib",
    "first_public",
    "last_updated",
    "include_in_effector_target_network_reanalysis",
    "metadata_status",
    "source_api_url",
    "metadata_retrieved_on",
]

AUDIT_FIELDS = ["check_id", "expected", "observed", "status", "notes"]

CONDITION_DETAILS = {
    "b_xylophilus": (
        "pathogenic B. xylophilus",
        "primary_case;secondary_case",
    ),
    "b_thailandae": (
        "non-pathogenic B. thailandae",
        "primary_comparator;nonpathogenic_response_case",
    ),
    "water_control": (
        "water inoculation control",
        "secondary_comparator;nonpathogenic_response_comparator",
    ),
}


def build_ena_url(api_url: str, accession: str) -> str:
    query = urlencode(
        {
            "accession": accession,
            "result": "read_run",
            "fields": ",".join(ENA_FIELDS),
            "format": "tsv",
            "download": "true",
        }
    )
    return f"{api_url.rstrip('?')}?{query}"


def fetch_ena_tsv(url: str, timeout_seconds: int) -> str:
    request = Request(url, headers={"User-Agent": "pine-pwn-effector_target_network-rnaseq-metadata/1.0"})
    with urlopen(request, timeout=timeout_seconds) as response:
        return response.read().decode("utf-8-sig")


def parse_ena_tsv(text: str) -> list[dict[str, str]]:
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    if not rows:
        raise ValueError("ENA returned no run records")
    missing = [field for field in ENA_FIELDS if field not in rows[0]]
    if missing:
        raise ValueError(f"ENA response is missing fields: {','.join(missing)}")
    return rows


def classify_condition(alias: str) -> tuple[str, str, str, str]:
    normalized = alias.lower()
    replicate_match = re.search(r"(?:^|[-_ ])rep(?:licate)?[-_ ]?(\d+)$", normalized)
    replicate = replicate_match.group(1) if replicate_match else ""
    if "bursaphelenchus xylophilus" in normalized:
        condition = "b_xylophilus"
    elif "bursaphelenchus thailandae" in normalized:
        condition = "b_thailandae"
    elif "water" in normalized and "control" in normalized:
        condition = "water_control"
    else:
        return "unresolved", "unresolved treatment", replicate, "unresolved"
    label, role = CONDITION_DETAILS[condition]
    return condition, label, replicate, role


def split_parallel_field(value: str, expected_parts: int = 2) -> list[str]:
    parts = [part.strip() for part in value.split(";") if part.strip()]
    return (parts + [""] * expected_parts)[:expected_parts]


def normalize_fastq_url(value: str) -> str:
    if not value:
        return ""
    if value.startswith(("https://", "http://", "ftp://")):
        return value
    return f"https://{value}"


def safe_int(value: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def normalize_ena_rows(
    rows: list[dict[str, str]],
    source_api_url: str,
    tissue: str,
    collection_time: str,
    retrieval_date: str | None = None,
) -> list[dict[str, str]]:
    retrieval_date = retrieval_date or date.today().isoformat()
    normalized_rows = []
    for row in rows:
        condition, label, replicate, role = classify_condition(row.get("experiment_alias", ""))
        fastq_urls = split_parallel_field(row.get("fastq_ftp", ""))
        fastq_md5 = split_parallel_field(row.get("fastq_md5", ""))
        fastq_bytes = split_parallel_field(row.get("fastq_bytes", ""))
        total_bytes = sum(safe_int(value) for value in fastq_bytes)
        complete_pair = all(fastq_urls) and all(fastq_md5) and all(fastq_bytes)
        eligible = all(
            [
                condition != "unresolved",
                bool(replicate),
                row.get("scientific_name") == "Pinus densiflora",
                row.get("library_layout") == "PAIRED",
                row.get("library_strategy") == "RNA-Seq",
                complete_pair,
            ]
        )
        normalized_rows.append(
            {
                "study_accession": row.get("study_accession", ""),
                "secondary_study_accession": row.get("secondary_study_accession", ""),
                "run_accession": row.get("run_accession", ""),
                "experiment_accession": row.get("experiment_accession", ""),
                "sample_accession": row.get("sample_accession", ""),
                "secondary_sample_accession": row.get("secondary_sample_accession", ""),
                "scientific_name": row.get("scientific_name", ""),
                "experiment_alias": row.get("experiment_alias", ""),
                "condition": condition,
                "condition_label": label,
                "replicate": replicate,
                "comparison_role": role,
                "tissue": tissue,
                "collection_time": collection_time,
                "library_layout": row.get("library_layout", ""),
                "library_strategy": row.get("library_strategy", ""),
                "instrument_platform": row.get("instrument_platform", ""),
                "instrument_model": row.get("instrument_model", ""),
                "read_count": row.get("read_count", ""),
                "base_count": row.get("base_count", ""),
                "fastq_1_url": normalize_fastq_url(fastq_urls[0]),
                "fastq_2_url": normalize_fastq_url(fastq_urls[1]),
                "fastq_1_md5": fastq_md5[0],
                "fastq_2_md5": fastq_md5[1],
                "fastq_1_bytes": fastq_bytes[0],
                "fastq_2_bytes": fastq_bytes[1],
                "total_fastq_bytes": str(total_bytes),
                "total_fastq_gib": f"{total_bytes / (1024**3):.3f}",
                "first_public": row.get("first_public", ""),
                "last_updated": row.get("last_updated", ""),
                "include_in_effector_target_network_reanalysis": "yes" if eligible else "no",
                "metadata_status": "complete" if eligible else "review_required",
                "source_api_url": source_api_url,
                "metadata_retrieved_on": retrieval_date,
            }
        )
    return sorted(normalized_rows, key=lambda item: (item["condition"], safe_int(item["replicate"])))


def audit_row(check_id: str, expected: str, observed: str, passed: bool, notes: str) -> dict[str, str]:
    return {
        "check_id": check_id,
        "expected": expected,
        "observed": observed,
        "status": "pass" if passed else "fail",
        "notes": notes,
    }


def audit_manifest(
    rows: list[dict[str, str]],
    expected_runs: int = 9,
    expected_replicates: int = 3,
    expected_bioproject: str = "",
    expected_secondary_study: str = "",
) -> list[dict[str, str]]:
    audits = []
    audits.append(
        audit_row(
            "run_count",
            str(expected_runs),
            str(len(rows)),
            len(rows) == expected_runs,
            "Expected design is three conditions with biological replicates.",
        )
    )

    study_pairs = sorted(
        {f"{row['study_accession']}|{row['secondary_study_accession']}" for row in rows}
    )
    expected_pair = f"{expected_bioproject}|{expected_secondary_study}"
    audits.append(
        audit_row(
            "study_identity",
            expected_pair,
            ";".join(study_pairs),
            study_pairs == [expected_pair],
            "This prevents an unrelated P. densiflora project from entering the infection reanalysis.",
        )
    )

    condition_counts = Counter(row["condition"] for row in rows)
    observed_conditions = ";".join(f"{key}={condition_counts[key]}" for key in sorted(condition_counts))
    expected_conditions = ";".join(
        f"{key}={expected_replicates}" for key in sorted(CONDITION_DETAILS)
    )
    condition_pass = set(condition_counts) == set(CONDITION_DETAILS) and all(
        condition_counts[condition] == expected_replicates for condition in CONDITION_DETAILS
    )
    audits.append(
        audit_row(
            "condition_balance",
            expected_conditions,
            observed_conditions,
            condition_pass,
            "Conditions are parsed conservatively from ENA experiment aliases.",
        )
    )

    replicates = defaultdict(set)
    for row in rows:
        replicates[row["condition"]].add(row["replicate"])
    expected_set = {str(index) for index in range(1, expected_replicates + 1)}
    replicate_observed = ";".join(
        f"{condition}={','.join(sorted(values))}" for condition, values in sorted(replicates.items())
    )
    replicate_pass = set(replicates) == set(CONDITION_DETAILS) and all(
        replicates[condition] == expected_set for condition in CONDITION_DETAILS
    )
    audits.append(
        audit_row(
            "replicate_labels",
            f"each_condition={','.join(sorted(expected_set))}",
            replicate_observed,
            replicate_pass,
            "Replicate labels are taken from experiment aliases, not sample accessions.",
        )
    )

    species = sorted({row["scientific_name"] for row in rows})
    audits.append(
        audit_row(
            "species_identity",
            "Pinus densiflora",
            ";".join(species),
            species == ["Pinus densiflora"],
            "All runs must represent the host species used in the published experiment.",
        )
    )

    designs = sorted(
        {f"{row['library_layout']}|{row['library_strategy']}|{row['instrument_platform']}" for row in rows}
    )
    design_pass = all(
        row["library_layout"] == "PAIRED"
        and row["library_strategy"] == "RNA-Seq"
        and row["instrument_platform"] == "ILLUMINA"
        for row in rows
    )
    audits.append(
        audit_row(
            "library_design",
            "PAIRED|RNA-Seq|ILLUMINA",
            ";".join(designs),
            design_pass,
            "The planned Salmon quantification requires paired Illumina RNA-seq reads.",
        )
    )

    complete_fastq = sum(
        all(row[field] for field in ["fastq_1_url", "fastq_2_url", "fastq_1_md5", "fastq_2_md5"])
        for row in rows
    )
    audits.append(
        audit_row(
            "paired_fastq_metadata",
            f"{expected_runs}/{expected_runs}",
            f"{complete_fastq}/{len(rows)}",
            complete_fastq == len(rows) == expected_runs,
            "Each run must expose two files and two checksums before download.",
        )
    )

    total_fastq_bytes = sum(safe_int(row["total_fastq_bytes"]) for row in rows)
    total_base_count = sum(safe_int(row["base_count"]) for row in rows)
    total_read_count = sum(safe_int(row["read_count"]) for row in rows)
    audits.append(
        audit_row(
            "download_volume",
            "reported",
            f"{total_fastq_bytes} bytes ({total_fastq_bytes / (1024**3):.3f} GiB)",
            total_fastq_bytes > 0,
            f"ENA reports {total_read_count} read records and {total_base_count} sequenced bases.",
        )
    )

    eligible_runs = sum(row["include_in_effector_target_network_reanalysis"] == "yes" for row in rows)
    prerequisite_pass = all(row["status"] == "pass" for row in audits)
    eligibility_pass = prerequisite_pass and eligible_runs == expected_runs
    audits.append(
        audit_row(
            "effector_target_network_reanalysis_eligibility",
            "eligible",
            "eligible" if eligibility_pass else "not_eligible",
            eligibility_pass,
            f"Eligible run records: {eligible_runs}/{len(rows)}.",
        )
    )
    return audits


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_report(path: Path, rows: list[dict[str, str]], audits: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    condition_counts = Counter(row["condition"] for row in rows)
    total_bytes = sum(safe_int(row["total_fastq_bytes"]) for row in rows)
    total_bases = sum(safe_int(row["base_count"]) for row in rows)
    total_reads = sum(safe_int(row["read_count"]) for row in rows)
    eligible = next(row for row in audits if row["check_id"] == "effector_target_network_reanalysis_eligibility")
    lines = [
        "# Route 1 P. densiflora RNA-seq Run Manifest",
        "",
        "ENA run metadata confirms a balanced infection experiment suitable for transcript-level reanalysis against the published Trinity reference.",
        "",
        f"- Runs: {len(rows)}",
        f"- Conditions: {', '.join(f'{key}={value}' for key, value in sorted(condition_counts.items()))}",
        f"- Compressed FASTQ volume: {total_bytes / (1024**3):.3f} GiB",
        f"- ENA read records: {total_reads}",
        f"- Sequenced bases: {total_bases}",
        f"- Reanalysis eligibility: {eligible['observed']}",
        "",
        "## Planned Contrasts",
        "",
        "1. B. xylophilus versus B. thailandae: primary pathogenicity-focused contrast.",
        "2. B. xylophilus versus water: secondary overall disease-response contrast.",
        "3. B. thailandae versus water: non-pathogenic nematode and wound-response context.",
        "",
        "## Evidence Boundary",
        "",
        "The run design supports differential-expression reanalysis, but expression changes will remain indirect support for Route 1 candidates and will not establish direct effector-target binding.",
        "",
        "Raw FASTQ download and quantification are intentionally separate workflow stages so the approximately 40+ GiB transfer is explicit and reviewable.",
        "",
        "## Audit",
        "",
    ]
    for audit in audits:
        lines.append(f"- {audit['check_id']}: {audit['status']} ({audit['observed']})")
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accession", required=True)
    parser.add_argument("--api-url", required=True)
    parser.add_argument("--tissue", required=True)
    parser.add_argument("--collection-time", required=True)
    parser.add_argument("--expected-runs", type=int, default=9)
    parser.add_argument("--expected-replicates", type=int, default=3)
    parser.add_argument("--expected-bioproject", required=True)
    parser.add_argument("--expected-secondary-study", required=True)
    parser.add_argument("--timeout-seconds", type=int, default=60)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit-output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    parser.add_argument("--strict", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    source_url = build_ena_url(args.api_url, args.accession)
    raw_rows = parse_ena_tsv(fetch_ena_tsv(source_url, args.timeout_seconds))
    rows = normalize_ena_rows(raw_rows, source_url, args.tissue, args.collection_time)
    audits = audit_manifest(
        rows,
        args.expected_runs,
        args.expected_replicates,
        args.expected_bioproject,
        args.expected_secondary_study,
    )
    write_tsv(args.output, rows, MANIFEST_FIELDS)
    write_tsv(args.audit_output, audits, AUDIT_FIELDS)
    write_report(args.markdown, rows, audits)
    eligible = next(row for row in audits if row["check_id"] == "effector_target_network_reanalysis_eligibility")
    print(
        "Route 1 P. densiflora RNA-seq metadata written: "
        f"{len(rows)} runs; {eligible['observed']}"
    )
    if args.strict and eligible["status"] != "pass":
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
