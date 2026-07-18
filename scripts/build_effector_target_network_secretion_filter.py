#!/usr/bin/env python3
"""Integrate DeepSig and TMHMM evidence for Route 1 effector representatives."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cluster-summary", type=Path, required=True)
    parser.add_argument("--priority-table", type=Path, required=True)
    parser.add_argument("--deepsig-gff3", type=Path, required=True)
    parser.add_argument("--tmhmm-tsv", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--representative-fasta", type=Path)
    parser.add_argument("--retained-fasta", type=Path)
    parser.add_argument("--min-signal-score", type=float, default=0.80)
    parser.add_argument("--galaxy-history-id", default="not_recorded")
    parser.add_argument("--deepsig-dataset-id", default="not_recorded")
    parser.add_argument("--tmhmm-dataset-id", default="not_recorded")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_deepsig(path: Path) -> dict[str, tuple[float, int, int]]:
    calls: dict[str, tuple[float, int, int]] = {}
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            if not raw_line.strip() or raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "Signal peptide":
                continue
            calls[fields[0]] = (float(fields[5]), int(fields[3]), int(fields[4]))
    return calls


def parse_tmhmm(path: Path) -> dict[str, int]:
    calls: dict[str, int] = {}
    for row in read_tsv(path):
        calls[row["#ID"]] = int(row["PredHel"])
    return calls


def write_retained_fasta(source: Path, target: Path, retained_ids: set[str]) -> None:
    records: dict[str, list[str]] = {}
    current_id: str | None = None
    with source.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                current_id = line[1:].split()[0]
                records[current_id] = [line]
            elif current_id is not None:
                records[current_id].append(line)
    missing = retained_ids - records.keys()
    if missing:
        raise ValueError(f"Retained FASTA IDs not found: {sorted(missing)}")
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="\n") as handle:
        for protein_id in sorted(retained_ids):
            handle.write("\n".join(records[protein_id]) + "\n")


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    clusters = read_tsv(args.cluster_summary)
    priority = {row["protein_id"]: row for row in read_tsv(args.priority_table)}
    deepsig = parse_deepsig(args.deepsig_gff3)
    tmhmm = parse_tmhmm(args.tmhmm_tsv)

    output_rows: list[dict[str, str]] = []
    for cluster in clusters:
        representative = cluster["representative_id"]
        if representative not in priority or representative not in tmhmm:
            raise ValueError(f"Missing representative input for {representative}")
        score, start, end = deepsig.get(representative, (None, None, None))
        predicted_helices = tmhmm[representative]
        signal_supported = score is not None and score >= args.min_signal_score
        soluble_supported = predicted_helices == 0
        status = "prioritize_for_structure_and_target_hypotheses" if signal_supported and soluble_supported else "deprioritize_pending_sequence_review"
        output_rows.append(
            {
                "representative_id": representative,
                "cluster_size": cluster["cluster_size"],
                "effector_classes": cluster["effector_classes"],
                "swissprot_annotation": priority[representative]["swissprot_annotation"],
                "deepsig_signal_peptide_score": "" if score is None else f"{score:.2f}",
                "deepsig_signal_peptide_start": "" if start is None else str(start),
                "deepsig_signal_peptide_end": "" if end is None else str(end),
                "tmhmm_predicted_helices": str(predicted_helices),
                "signal_peptide_supported": str(signal_supported).lower(),
                "no_predicted_transmembrane_helix": str(soluble_supported).lower(),
                "secretion_filter_status": status,
                "member_inference_scope": "Representative prediction; cluster members require individual prediction before sequence-level claims",
                "claim_ceiling": "Sequence-prediction support only; does not validate secretion, effector activity, host binding, or host targeting",
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    if bool(args.representative_fasta) != bool(args.retained_fasta):
        raise ValueError("--representative-fasta and --retained-fasta must be supplied together")
    if args.representative_fasta and args.retained_fasta:
        retained_ids = {
            row["representative_id"]
            for row in output_rows
            if row["secretion_filter_status"].startswith("prioritize")
        }
        write_retained_fasta(args.representative_fasta, args.retained_fasta, retained_ids)

    total = len(output_rows)
    signal = sum(row["signal_peptide_supported"] == "true" for row in output_rows)
    retained = sum(row["secretion_filter_status"].startswith("prioritize") for row in output_rows)
    by_class: dict[str, int] = {}
    for row in output_rows:
        if row["secretion_filter_status"].startswith("prioritize"):
            for effector_class in row["effector_classes"].split(";"):
                by_class[effector_class] = by_class.get(effector_class, 0) + 1
    report = [
        "# Route 1 Candidate Secretion Filter",
        "",
        "## Inputs and rule",
        "",
        f"- {total} MMseqs2 sequence-cluster representatives were evaluated with DeepSig (eukaryote mode) and TMHMM 2.0 through Galaxy.",
        f"- Galaxy history: `{args.galaxy_history_id}`; DeepSig output dataset: `{args.deepsig_dataset_id}`; TMHMM output dataset: `{args.tmhmm_dataset_id}`.",
        f"- A representative is retained when DeepSig reports an N-terminal signal peptide with score >= {args.min_signal_score:.2f} and TMHMM predicts zero transmembrane helices.",
        f"- {signal}/{total} representatives have DeepSig signal-peptide support; {retained}/{total} pass the combined soluble-secreted filter.",
        "",
        "## Retained representative classes",
        "",
    ]
    report.extend(f"- {name}: {count}" for name, count in sorted(by_class.items()))
    report.extend(
        [
            "",
            "## Interpretation boundary",
            "",
            "These are protein-sequence predictions, not experimental validation of secretion or effector function. The filter does not establish a physical interaction with pine proteins. A clustered member inherits no sequence-level call from its representative; it remains a prioritization aid until individually predicted.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(report) + "\n", encoding="utf-8")
    logging.info("Wrote %d representative calls; retained %d", total, retained)


if __name__ == "__main__":
    main()
