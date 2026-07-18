#!/usr/bin/env python3
"""Build traceable FastQC and Salmon audit tables for manuscript datasets."""

from __future__ import annotations

import argparse
import csv
import logging
import re
from collections import Counter, defaultdict
from pathlib import Path


MODULES = (
    "Basic Statistics",
    "Per base sequence quality",
    "Per sequence quality scores",
    "Per base sequence content",
    "Per sequence GC content",
    "Per base N content",
    "Sequence Length Distribution",
    "Sequence Duplication Levels",
    "Overrepresented sequences",
    "Adapter Content",
)


def slug(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def parse_fastqc(path: Path) -> dict[str, str]:
    basic: dict[str, str] = {}
    statuses: dict[str, str] = {}
    deduplicated = None
    current_module = None
    in_basic = False

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.rstrip("\n")
        if line.startswith(">>") and not line.startswith(">>END_MODULE"):
            parts = line[2:].split("\t")
            current_module = parts[0]
            statuses[current_module] = parts[1] if len(parts) > 1 else "unknown"
            in_basic = current_module == "Basic Statistics"
            continue
        if line.startswith(">>END_MODULE"):
            current_module = None
            in_basic = False
            continue
        if current_module == "Sequence Duplication Levels" and line.startswith(
            "#Total Deduplicated Percentage\t"
        ):
            deduplicated = float(line.split("\t", 1)[1])
        elif in_basic and line and not line.startswith("#"):
            parts = line.split("\t", 1)
            if len(parts) == 2:
                basic[parts[0]] = parts[1]

    accession_match = re.search(r"(SRR\d+)", path.name)
    if not accession_match:
        accession_match = re.search(r"(SRR\d+)", basic.get("Filename", ""))
    if not accession_match:
        raise ValueError(f"Could not identify an SRR accession in {path}")

    row = {
        "run_accession": accession_match.group(1),
        "fastqc_file": path.as_posix(),
        "source_filename": basic.get("Filename", "NA"),
        "total_sequences": basic.get("Total Sequences", "NA"),
        "sequence_length": basic.get("Sequence length", "NA"),
        "gc_pct": basic.get("%GC", "NA"),
        "duplication_pct": "NA" if deduplicated is None else f"{100.0 - deduplicated:.2f}",
    }
    for module in MODULES:
        row[slug(module)] = statuses.get(module, "not_reported")
    return row


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"No rows available for {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def format_range(values: list[float], digits: int = 2) -> str:
    return f"{min(values):.{digits}f}-{max(values):.{digits}f}"


def build_audit(
    dataset_dirs: dict[str, Path],
    provenance_path: Path,
    run_output: Path,
    summary_output: Path,
    report_path: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    provenance = read_tsv(provenance_path)
    provenance_by_key = {(row["dataset_id"], row["run_accession"]): row for row in provenance}
    run_rows: list[dict[str, str]] = []

    for dataset_id, directory in dataset_dirs.items():
        fastqc_paths = sorted(directory.glob("*.fastqc.txt"))
        if not fastqc_paths:
            raise ValueError(f"No FastQC text files found for {dataset_id}: {directory}")
        for path in fastqc_paths:
            parsed = parse_fastqc(path)
            key = (dataset_id, parsed["run_accession"])
            if key not in provenance_by_key:
                raise ValueError(f"Missing quantification provenance for {dataset_id}/{parsed['run_accession']}")
            source = provenance_by_key[key]
            run_rows.append(
                {
                    "dataset_id": dataset_id,
                    "species": source["species"],
                    "sample_id": source["sample_id"],
                    "run_accession": parsed["run_accession"],
                    "fastqc_file": parsed["fastqc_file"],
                    "source_filename": parsed["source_filename"],
                    "platform": source["platform"],
                    "instrument_model": source["instrument_model"],
                    "layout": source["layout"],
                    "total_sequences": parsed["total_sequences"],
                    "sequence_length": parsed["sequence_length"],
                    "gc_pct": parsed["gc_pct"],
                    "duplication_pct": parsed["duplication_pct"],
                    "salmon_mapping_rate_pct": source["salmon_mapping_rate_pct"],
                    **{slug(module): parsed[slug(module)] for module in MODULES},
                }
            )

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in run_rows:
        grouped[row["dataset_id"]].append(row)

    summaries: list[dict[str, str]] = []
    for dataset_id in sorted(grouped):
        rows = grouped[dataset_id]
        unique_runs = sorted({row["run_accession"] for row in rows})
        mapping_by_run = {
            row["run_accession"]: float(row["salmon_mapping_rate_pct"]) for row in rows
        }
        mapping_values = list(mapping_by_run.values())
        duplication_values = [
            float(row["duplication_pct"]) for row in rows if row["duplication_pct"] != "NA"
        ]
        gc_values = [float(row["gc_pct"]) for row in rows if row["gc_pct"] != "NA"]
        lengths = sorted({row["sequence_length"] for row in rows})
        summary = {
            "dataset_id": dataset_id,
            "species": rows[0]["species"],
            "platform": rows[0]["platform"],
            "instrument_model": rows[0]["instrument_model"],
            "layout": rows[0]["layout"],
            "library_count": str(len(unique_runs)),
            "fastqc_file_count": str(len(rows)),
            "read_length": ";".join(lengths),
            "gc_range_pct": format_range(gc_values, 0),
            "duplication_range_pct": format_range(duplication_values),
            "mapping_rate_mean_pct": f"{sum(mapping_values) / len(mapping_values):.2f}",
            "mapping_rate_range_pct": format_range(mapping_values),
        }
        for module in MODULES:
            key = slug(module)
            counts = Counter(row[key] for row in rows)
            summary[f"{key}_counts"] = ";".join(
                f"{status}:{counts[status]}" for status in ("pass", "warn", "fail", "not_reported") if counts[status]
            )
        summaries.append(summary)

    write_tsv(run_output, run_rows)
    write_tsv(summary_output, summaries)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# RNA-seq upstream quality audit",
        "",
        "FastQC results are counted per FASTQ file; library counts are unique SRA runs. Salmon mapping rates are counted once per library.",
        "",
    ]
    for row in summaries:
        lines.append(
            f"- `{row['dataset_id']}`: {row['library_count']} libraries, {row['fastqc_file_count']} FastQC files, "
            f"read length {row['read_length']} nt, Salmon mapping {row['mapping_rate_mean_pct']}% "
            f"(range {row['mapping_rate_range_pct']}%)."
        )
    lines.extend(
        [
            "",
            "For the *P. densiflora* DESeq2 contrasts, no explicit count prefilter was applied; DESeq2 independent filtering remained enabled. FastQC warnings or failures are retained as audit evidence and are not silently converted to exclusions.",
        ]
    )
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    logging.info("Wrote %d FastQC records and %d dataset summaries", len(run_rows), len(summaries))
    return run_rows, summaries


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pmas-fastqc-dir", type=Path, required=True)
    parser.add_argument("--pden-fastqc-dir", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    parser.add_argument("--run-output", type=Path, required=True)
    parser.add_argument("--summary-output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--log", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    handlers: list[logging.Handler] = [logging.StreamHandler()]
    if args.log:
        args.log.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(args.log, encoding="utf-8"))
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", handlers=handlers)
    build_audit(
        {"pmas_1dpi": args.pmas_fastqc_dir, "pden_full": args.pden_fastqc_dir},
        args.provenance,
        args.run_output,
        args.summary_output,
        args.report,
    )


if __name__ == "__main__":
    main()
