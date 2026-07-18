#!/usr/bin/env python3
"""Integrate heuristic, DeepSig, TMHMM, and functional secretome evidence."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


FIELDS = [
    "protein_id",
    "annotation_text",
    "heuristic_secretome_candidate",
    "signal_peptide_tool",
    "signal_peptide_score",
    "signal_peptide_start",
    "signal_peptide_end",
    "signal_peptide_status",
    "transmembrane_tool",
    "tmhmm_predicted_helices",
    "signal_overlapping_tm_count",
    "mature_protein_tm_count",
    "post_signal_tm_status",
    "standard_prediction_status",
    "standard_secreted_soluble_status",
    "domain_support",
    "merops_support",
    "dbcan_support",
    "published_secretome_support",
    "evidence_tier",
    "claim_ceiling",
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


def parse_deepsig(path: Path) -> tuple[set[str], dict[str, tuple[float, int, int]]]:
    evaluated: set[str] = set()
    signals: dict[str, tuple[float, int, int]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 8:
                raise ValueError(f"Invalid DeepSig GFF3 row: {line}")
            protein_id = fields[0]
            evaluated.add(protein_id)
            if fields[2].lower().replace("_", " ") != "signal peptide":
                continue
            signals[protein_id] = (float(fields[5]), int(fields[3]), int(fields[4]))
    return evaluated, signals


def parse_tmhmm(path: Path) -> dict[str, dict[str, str]]:
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames:
            reader.fieldnames = [field.lstrip("#") for field in reader.fieldnames]
        return {row["ID"]: row for row in reader}


def topology_helices(topology: str) -> list[tuple[int, int]]:
    return [(int(start), int(end)) for start, end in re.findall(r"(\d+)-(\d+)", topology)]


def truthy(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes", "y", "supported", "present"}


def functional_by_id(path: Path | None) -> dict[str, dict[str, str]]:
    if path is None:
        return {}
    return {row["protein_id"]: row for row in read_tsv(path)}


def build_secretome_audit(
    heuristic_path: Path,
    deepsig_path: Path,
    tmhmm_path: Path,
    functional_path: Path | None,
    output_path: Path,
    report_path: Path,
    min_signal_score: float = 0.80,
    signal_overlap_tolerance: int = 5,
) -> tuple[list[dict[str, str]], dict[str, int]]:
    heuristic_rows = read_tsv(heuristic_path)
    evaluated, signals = parse_deepsig(deepsig_path)
    tmhmm = parse_tmhmm(tmhmm_path)
    functional = functional_by_id(functional_path)
    rows: list[dict[str, str]] = []

    for source in heuristic_rows:
        protein_id = source["protein_id"]
        signal = signals.get(protein_id)
        signal_supported = signal is not None and signal[0] >= min_signal_score
        tm = tmhmm.get(protein_id)
        complete = protein_id in evaluated and tm is not None
        helices = topology_helices(tm.get("Topology", "")) if tm else []
        if tm and len(helices) != int(tm.get("PredHel", "0")):
            raise ValueError(f"Could not reconcile TMHMM topology for {protein_id}: {tm.get('Topology', '')}")
        signal_end = signal[2] if signal_supported and signal else None
        overlap_count = sum(
            end <= signal_end + signal_overlap_tolerance
            for _, end in helices
        ) if signal_end is not None else 0
        mature_tm_count = len(helices) - overlap_count
        standard_pass = complete and signal_supported and mature_tm_count == 0

        function = functional.get(protein_id, {})
        support = {
            field: function.get(field, "not_evaluated")
            for field in (
                "domain_support", "merops_support", "dbcan_support", "published_secretome_support"
            )
        }
        has_functional_support = any(truthy(value) for value in support.values())
        if standard_pass and has_functional_support:
            tier = "functionally_supported_candidate"
        elif standard_pass:
            tier = "high_confidence_secreted_protein"
        else:
            tier = "exploratory_sequence_candidate"

        row = {
            "protein_id": protein_id,
            "annotation_text": source.get("annotation_text", ""),
            "heuristic_secretome_candidate": source.get("secretome_candidate", ""),
            "signal_peptide_tool": "DeepSig_1.2.5",
            "signal_peptide_score": "" if signal is None else f"{signal[0]:.2f}",
            "signal_peptide_start": "" if signal is None else str(signal[1]),
            "signal_peptide_end": "" if signal is None else str(signal[2]),
            "signal_peptide_status": (
                "missing" if protein_id not in evaluated else "supported" if signal_supported else "not_supported"
            ),
            "transmembrane_tool": "TMHMM_2.0_wrapper_0.0.17",
            "tmhmm_predicted_helices": "" if tm is None else tm.get("PredHel", ""),
            "signal_overlapping_tm_count": "" if tm is None else str(overlap_count),
            "mature_protein_tm_count": "" if tm is None else str(mature_tm_count),
            "post_signal_tm_status": (
                "missing" if tm is None else "absent" if mature_tm_count == 0 else "present"
            ),
            "standard_prediction_status": "complete" if complete else "incomplete",
            "standard_secreted_soluble_status": "pass" if standard_pass else "fail" if complete else "incomplete",
            **support,
            "evidence_tier": tier,
            "claim_ceiling": "sequence_prediction_only_not_validated_secretion_or_effector_activity",
        }
        rows.append(row)

    rows.sort(key=lambda row: row["protein_id"])
    summary = {
        "proteins": len(rows),
        "standard_prediction_complete": sum(row["standard_prediction_status"] == "complete" for row in rows),
        "standard_secreted_soluble": sum(row["standard_secreted_soluble_status"] == "pass" for row in rows),
        "functionally_supported": sum(row["evidence_tier"] == "functionally_supported_candidate" for row in rows),
    }
    write_tsv(output_path, rows)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        "# Full-proteome secretome evidence audit\n\n"
        f"- Proteins in the source set: {summary['proteins']}.\n"
        f"- Proteins with complete DeepSig/TMHMM records: {summary['standard_prediction_complete']}.\n"
        f"- Proteins passing the signal-peptide plus no mature-protein TM rule: {summary['standard_secreted_soluble']}.\n"
        f"- Passing proteins with at least one recorded functional support channel: {summary['functionally_supported']}.\n\n"
        "A TMHMM helix ending no more than "
        f"{signal_overlap_tolerance} residues after the DeepSig signal-peptide end was treated as signal-overlapping. "
        "All categories are sequence-based priorities and are not experimental validation of secretion, effector activity, host binding, or cleavage.\n",
        encoding="utf-8",
    )
    return rows, summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--heuristic", type=Path, required=True)
    parser.add_argument("--deepsig", type=Path, required=True)
    parser.add_argument("--tmhmm", type=Path, required=True)
    parser.add_argument("--functional", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--min-signal-score", type=float, default=0.80)
    parser.add_argument("--signal-overlap-tolerance", type=int, default=5)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    _, summary = build_secretome_audit(
        args.heuristic, args.deepsig, args.tmhmm, args.functional,
        args.output, args.report, args.min_signal_score, args.signal_overlap_tolerance,
    )
    print(f"Audited {summary['proteins']} proteins; {summary['standard_secreted_soluble']} passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
