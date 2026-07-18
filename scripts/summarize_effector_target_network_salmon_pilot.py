#!/usr/bin/env python3
"""Summarize Route 1 candidate abundance in a single Salmon pilot sample."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


FIELDS = [
    "sample_id",
    "condition",
    "transcript_id",
    "orthogroup_ids",
    "effector_target_network_seed_symbols",
    "mapping_status",
    "best_effector_target_network_candidate_pd_gene",
    "isoforms_in_reference",
    "isoforms_with_positive_reads",
    "summed_isoform_tpm",
    "summed_isoform_numreads",
    "gene_tpm",
    "gene_numreads",
    "pilot_abundance_status",
    "claim_ceiling",
]

CLAIM_CEILING = (
    "Single-sample Salmon abundance confirms quantifiability only; it is not differential-expression, "
    "infection-response, resistance, or direct effector-target evidence"
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def trinity_gene_id(transcript_id: str) -> str:
    return re.sub(r"_i\d+$", "", transcript_id)


def safe_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def format_number(value: float) -> str:
    return f"{value:.10g}"


def summarize_effector_target_network_salmon_pilot(
    candidates_path: Path,
    transcript_quant_path: Path,
    gene_quant_path: Path,
    output_path: Path,
    markdown_path: Path,
    sample_id: str,
    condition: str,
) -> list[dict[str, str]]:
    isoforms = defaultdict(list)
    for row in read_tsv(transcript_quant_path):
        isoforms[trinity_gene_id(row["Name"])].append(row)
    genes = {row["Name"]: row for row in read_tsv(gene_quant_path)}

    rows = []
    for candidate in read_tsv(candidates_path):
        transcript_id = candidate["transcript_id"]
        matched_isoforms = isoforms.get(transcript_id, [])
        gene_row = genes.get(transcript_id)
        positive_isoforms = sum(safe_float(row.get("NumReads", "0")) > 0 for row in matched_isoforms)
        summed_tpm = sum(safe_float(row.get("TPM", "0")) for row in matched_isoforms)
        summed_reads = sum(safe_float(row.get("NumReads", "0")) for row in matched_isoforms)
        gene_tpm = safe_float(gene_row.get("TPM", "0")) if gene_row else 0.0
        gene_reads = safe_float(gene_row.get("NumReads", "0")) if gene_row else 0.0
        if not matched_isoforms or gene_row is None:
            status = "missing_from_pilot_quantification"
        elif gene_reads > 0:
            status = "detected_in_pilot_sample"
        else:
            status = "zero_reads_in_pilot_sample"
        rows.append(
            {
                "sample_id": sample_id,
                "condition": condition,
                "transcript_id": transcript_id,
                "orthogroup_ids": candidate.get("orthogroup_ids", ""),
                "effector_target_network_seed_symbols": candidate.get("effector_target_network_seed_symbols", ""),
                "mapping_status": candidate.get("mapping_status", ""),
                "best_effector_target_network_candidate_pd_gene": candidate.get("best_effector_target_network_candidate_pd_gene", ""),
                "isoforms_in_reference": str(len(matched_isoforms)),
                "isoforms_with_positive_reads": str(positive_isoforms),
                "summed_isoform_tpm": format_number(summed_tpm),
                "summed_isoform_numreads": format_number(summed_reads),
                "gene_tpm": format_number(gene_tpm),
                "gene_numreads": format_number(gene_reads),
                "pilot_abundance_status": status,
                "claim_ceiling": CLAIM_CEILING,
            }
        )
    write_tsv(output_path, rows)
    write_markdown(markdown_path, rows, transcript_quant_path, gene_quant_path)
    return rows


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(
    path: Path,
    rows: list[dict[str, str]],
    transcript_quant_path: Path,
    gene_quant_path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    status_counts = Counter(row["pilot_abundance_status"] for row in rows)
    current_pd_rows = [row for row in rows if row["mapping_status"] == "current_pd_effector_target_network_candidate_hit"]
    current_pd_detected = sum(row["pilot_abundance_status"] == "detected_in_pilot_sample" for row in current_pd_rows)
    lines = [
        "# Route 1 Salmon Pilot Candidate Abundance",
        "",
        f"- Sample: {rows[0]['sample_id'] if rows else 'NA'}",
        f"- Condition: {rows[0]['condition'] if rows else 'NA'}",
        f"- Transcript quantification: `{transcript_quant_path}`",
        f"- Gene quantification: `{gene_quant_path}`",
        f"- Candidate Trinity genes: {len(rows)}",
        f"- Detected candidate genes: {status_counts.get('detected_in_pilot_sample', 0)}",
        f"- Zero-read candidate genes: {status_counts.get('zero_reads_in_pilot_sample', 0)}",
        f"- Missing candidate genes: {status_counts.get('missing_from_pilot_quantification', 0)}",
        f"- Current-Pd Route 1 bridges detected: {current_pd_detected}/{len(current_pd_rows)}",
        "",
        "Evidence boundary: this is single-sample abundance from a water-control pilot. It verifies that candidate IDs can pass through the Galaxy/Salmon bridge, but it cannot support differential-expression or infection-response claims.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--transcript-quant", type=Path, required=True)
    parser.add_argument("--gene-quant", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--condition", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = summarize_effector_target_network_salmon_pilot(
        args.candidates,
        args.transcript_quant,
        args.gene_quant,
        args.output,
        args.markdown,
        args.sample_id,
        args.condition,
    )
    detected = sum(row["pilot_abundance_status"] == "detected_in_pilot_sample" for row in rows)
    print(f"Route 1 Salmon pilot summarized: {detected}/{len(rows)} candidate genes detected")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
