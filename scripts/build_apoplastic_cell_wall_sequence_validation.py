#!/usr/bin/env python3
"""Integrate family annotation, DeepSig, and TMHMM for Route 2 host targets."""

from __future__ import annotations

import argparse
import csv
import logging
import re
from collections import defaultdict
from pathlib import Path


FIELDS = [
    "orthogroup_id", "expected_host_family", "species_id", "gene_id", "sequence_length",
    "expression_selected_member", "expression_evidence_tier", "annotation_source", "annotation_text",
    "family_annotation_status", "deepsig_signal_score", "deepsig_signal_start", "deepsig_signal_end",
    "tmhmm_predicted_helices", "tmhmm_topology", "signal_overlapping_tmhmm_helices",
    "mature_protein_tm_helices", "sequence_localization_status", "claim_ceiling",
]
SUMMARY_FIELDS = [
    "orthogroup_id", "expected_host_family", "protein_count", "family_annotation_supported_count",
    "family_related_annotation_count", "annotation_missing_count", "signal_supported_count",
    "secreted_soluble_supported_count", "post_signal_tm_supported_count", "n_terminal_ambiguous_count", "expression_selected_gene_id",
    "expression_selected_family_status", "expression_selected_localization_status",
]
CLAIM = (
    "Sequence prediction and homology annotation only; does not validate in vivo pine localization, "
    "M8 binding, proteolytic cleavage, or causal hydraulic failure"
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_deepsig(path: Path) -> dict[str, tuple[float, int, int]]:
    calls: dict[str, tuple[float, int, int]] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) == 9 and fields[2] == "Signal peptide":
                calls[fields[0]] = (float(fields[5]), int(fields[3]), int(fields[4]))
    return calls


def parse_tmhmm(path: Path) -> dict[str, dict[str, str]]:
    return {row["#ID"]: row for row in read_tsv(path)}


def family_status(expected: str, annotation: str) -> str:
    text = annotation.lower()
    if not text:
        return "annotation_missing"
    if expected == "XTH" and ("xth" in text or "xyloglucan endotrans" in text):
        return "family_annotation_supported"
    if expected == "PMEI" and ("pmei" in text or "pectinesterase inhibitor" in text):
        return "family_annotation_supported"
    if expected == "PMEI" and "21 kda protein" in text:
        return "family_related_annotation_only"
    return "family_annotation_conflict"


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build_validation(
    manifest_path: Path,
    shortlist_path: Path,
    annotations_path: Path,
    deepsig_path: Path,
    tmhmm_path: Path,
    output_path: Path,
    summary_path: Path,
    report_path: Path,
    min_signal_score: float = 0.80,
    signal_overlap_tolerance: int = 5,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    families = {row["host_orthogroup_id"]: row["host_family"] for row in read_tsv(shortlist_path)}
    annotations = {row["gene_id"]: row for row in read_tsv(annotations_path)}
    deepsig = parse_deepsig(deepsig_path)
    tmhmm = parse_tmhmm(tmhmm_path)
    rows: list[dict[str, str]] = []

    for protein in read_tsv(manifest_path):
        gene_id = protein["gene_id"]
        if gene_id not in tmhmm:
            raise ValueError(f"Missing TMHMM result for {gene_id}")
        expected = families[protein["orthogroup_id"]]
        annotation = annotations.get(gene_id, {})
        signal = deepsig.get(gene_id)
        signal_supported = signal is not None and signal[0] >= min_signal_score
        topology = tmhmm[gene_id].get("Topology", "")
        helices = [(int(start), int(end)) for start, end in re.findall(r"(\d+)-(\d+)", topology)]
        if len(helices) != int(tmhmm[gene_id].get("PredHel", "0")):
            raise ValueError(f"Could not reconcile TMHMM topology for {gene_id}: {topology}")
        overlap = 0
        mature = len(helices)
        if signal_supported:
            overlap = sum(start <= signal[2] + signal_overlap_tolerance for start, _ in helices)
            mature -= overlap
        if signal_supported and mature == 0:
            localization = "supported_secreted_soluble"
        elif signal_supported:
            localization = "signal_supported_with_mature_tm"
        elif helices and all(start <= 25 and end <= 45 for start, end in helices):
            localization = "n_terminal_signal_or_anchor_ambiguous"
        elif helices:
            localization = "tm_supported_without_signal"
        else:
            localization = "no_signal_peptide_support"
        rows.append({
            "orthogroup_id": protein["orthogroup_id"], "expected_host_family": expected,
            "species_id": protein["species_id"], "gene_id": gene_id,
            "sequence_length": protein["sequence_length"],
            "expression_selected_member": protein["expression_selected_member"],
            "expression_evidence_tier": protein["expression_evidence_tier"],
            "annotation_source": annotation.get("annotation_source", ""),
            "annotation_text": annotation.get("annotation_text", ""),
            "family_annotation_status": family_status(expected, annotation.get("annotation_text", "")),
            "deepsig_signal_score": "" if signal is None else f"{signal[0]:.2f}",
            "deepsig_signal_start": "" if signal is None else str(signal[1]),
            "deepsig_signal_end": "" if signal is None else str(signal[2]),
            "tmhmm_predicted_helices": tmhmm[gene_id].get("PredHel", ""),
            "tmhmm_topology": topology, "signal_overlapping_tmhmm_helices": str(overlap),
            "mature_protein_tm_helices": str(mature), "sequence_localization_status": localization,
            "claim_ceiling": CLAIM,
        })
    rows.sort(key=lambda row: (row["orthogroup_id"], row["species_id"], row["gene_id"]))

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["orthogroup_id"]].append(row)
    summaries: list[dict[str, str]] = []
    for orthogroup, members in sorted(grouped.items()):
        selected = [row for row in members if row["expression_selected_member"] == "true"]
        summaries.append({
            "orthogroup_id": orthogroup, "expected_host_family": members[0]["expected_host_family"],
            "protein_count": str(len(members)),
            "family_annotation_supported_count": str(sum(row["family_annotation_status"] == "family_annotation_supported" for row in members)),
            "family_related_annotation_count": str(sum(row["family_annotation_status"] == "family_related_annotation_only" for row in members)),
            "annotation_missing_count": str(sum(row["family_annotation_status"] == "annotation_missing" for row in members)),
            "signal_supported_count": str(sum(row["deepsig_signal_score"] != "" for row in members)),
            "secreted_soluble_supported_count": str(sum(row["sequence_localization_status"] == "supported_secreted_soluble" for row in members)),
            "post_signal_tm_supported_count": str(sum(row["sequence_localization_status"] == "signal_supported_with_mature_tm" for row in members)),
            "n_terminal_ambiguous_count": str(sum(row["sequence_localization_status"] == "n_terminal_signal_or_anchor_ambiguous" for row in members)),
            "expression_selected_gene_id": ";".join(row["gene_id"] for row in selected),
            "expression_selected_family_status": ";".join(row["family_annotation_status"] for row in selected),
            "expression_selected_localization_status": ";".join(row["sequence_localization_status"] for row in selected),
        })
    write_tsv(output_path, FIELDS, rows)
    write_tsv(summary_path, SUMMARY_FIELDS, summaries)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        "# Route 2 target sequence validation\n\n"
        f"- Proteins evaluated: {len(rows)} across {len(summaries)} orthogroups.\n"
        f"- DeepSig support threshold: {min_signal_score:.2f}.\n"
        f"- A TMHMM helix ending no more than {signal_overlap_tolerance} residues after the DeepSig signal-peptide end is treated as signal-overlapping, not a mature-protein transmembrane helix.\n\n"
        "Family annotations and localization predictions are sequence-level support only. They do not establish "
        "that M8 binds or cleaves any pine protein.\n",
        encoding="utf-8",
    )
    logging.info("Validated %d proteins across %d orthogroups", len(rows), len(summaries))
    return rows, summaries


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--shortlist", type=Path, required=True)
    parser.add_argument("--annotations", type=Path, required=True)
    parser.add_argument("--deepsig", type=Path, required=True)
    parser.add_argument("--tmhmm", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--min-signal-score", type=float, default=0.80)
    parser.add_argument("--signal-overlap-tolerance", type=int, default=5)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    build_validation(
        args.manifest, args.shortlist, args.annotations, args.deepsig, args.tmhmm,
        args.output, args.summary, args.report, args.min_signal_score, args.signal_overlap_tolerance,
    )


if __name__ == "__main__":
    main()
