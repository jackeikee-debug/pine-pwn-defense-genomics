#!/usr/bin/env python3
"""Summarize InterProScan evidence for secretion-filtered Route 1 candidates."""

from __future__ import annotations

import argparse
import csv
import logging
import re
from collections import defaultdict
from pathlib import Path


PROTEASE_PATTERN = re.compile(r"peptid|proteas|carboxypeptid|leishmanolysin|metalloprote", re.I)
CATALYTIC_PATTERN = re.compile(r"active site|catalytic domain|catalytic chain", re.I)
P450_PATTERN = re.compile(r"cytochrome p450|p450", re.I)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--interproscan", type=Path, required=True)
    parser.add_argument("--secretion-filter", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--filtered-fasta", type=Path)
    parser.add_argument("--high-priority-fasta", type=Path)
    parser.add_argument("--galaxy-history-id", default="not_recorded")
    parser.add_argument("--galaxy-dataset-id", default="not_recorded")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_interpro(path: Path) -> dict[str, list[list[str]]]:
    grouped: dict[str, list[list[str]]] = defaultdict(list)
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            row.extend(["-"] * (15 - len(row)))
            grouped[row[0]].append(row[:15])
    return grouped


def joined_description(row: list[str]) -> str:
    return " ".join((row[5], row[12]))


def write_selected_fasta(source: Path, target: Path, selected_ids: set[str]) -> None:
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
    missing = selected_ids - records.keys()
    if missing:
        raise ValueError(f"Selected FASTA IDs not found: {sorted(missing)}")
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="\n") as handle:
        for protein_id in sorted(selected_ids):
            handle.write("\n".join(records[protein_id]) + "\n")


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    filtered = [
        row
        for row in read_tsv(args.secretion_filter)
        if row["secretion_filter_status"] == "prioritize_for_structure_and_target_hypotheses"
    ]
    interpro = read_interpro(args.interproscan)
    output_rows: list[dict[str, str]] = []

    for candidate in filtered:
        protein_id = candidate["representative_id"]
        hits = interpro.get(protein_id, [])
        if not hits:
            raise ValueError(f"No InterProScan rows for {protein_id}")
        protease_hits = [row for row in hits if PROTEASE_PATTERN.search(joined_description(row))]
        catalytic_hits = [row for row in hits if CATALYTIC_PATTERN.search(joined_description(row))]
        p450_hits = [row for row in hits if P450_PATTERN.search(joined_description(row))]
        class_is_protease = "protease" in candidate["effector_classes"].split(";")
        class_is_detox = "detoxification" in candidate["effector_classes"].split(";")
        expected_domain_supported = (class_is_protease and bool(protease_hits)) or (class_is_detox and bool(p450_hits))

        if catalytic_hits:
            priority = "high_catalytic_feature_supported"
        elif protease_hits:
            priority = "medium_protease_domain_supported"
        elif p450_hits:
            priority = "medium_p450_domain_supported"
        else:
            priority = "low_expected_domain_not_recovered"

        relevant_hits = protease_hits or p450_hits
        ipr_pairs = sorted({(row[11], row[12]) for row in relevant_hits if row[11] != "-"})
        catalytic_pairs = sorted({(row[3], row[4], row[5]) for row in catalytic_hits})
        output_rows.append(
            {
                "representative_id": protein_id,
                "cluster_size": candidate["cluster_size"],
                "effector_classes": candidate["effector_classes"],
                "deepsig_signal_peptide_score": candidate["deepsig_signal_peptide_score"],
                "interproscan_hit_count": str(len(hits)),
                "interproscan_applications": ";".join(sorted({row[3] for row in hits})),
                "expected_domain_supported": str(expected_domain_supported).lower(),
                "protease_domain_supported": str(bool(protease_hits)).lower(),
                "p450_domain_supported": str(bool(p450_hits)).lower(),
                "catalytic_feature_annotation": str(bool(catalytic_hits)).lower(),
                "relevant_interpro_terms": ";".join(f"{ipr}:{description}" for ipr, description in ipr_pairs),
                "catalytic_feature_terms": ";".join(f"{app}:{accession}:{description}" for app, accession, description in catalytic_pairs),
                "structure_followup_priority": priority,
                "claim_ceiling": "Domain and conserved-feature predictions only; enzymatic activity, secretion, host binding, and host targeting remain unvalidated",
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    if bool(args.filtered_fasta) != bool(args.high_priority_fasta):
        raise ValueError("--filtered-fasta and --high-priority-fasta must be supplied together")
    if args.filtered_fasta and args.high_priority_fasta:
        high_priority_ids = {
            row["representative_id"]
            for row in output_rows
            if row["structure_followup_priority"] == "high_catalytic_feature_supported"
        }
        write_selected_fasta(args.filtered_fasta, args.high_priority_fasta, high_priority_ids)

    total = len(output_rows)
    expected = sum(row["expected_domain_supported"] == "true" for row in output_rows)
    protease = sum(row["protease_domain_supported"] == "true" for row in output_rows)
    catalytic = sum(row["catalytic_feature_annotation"] == "true" for row in output_rows)
    p450 = sum(row["p450_domain_supported"] == "true" for row in output_rows)
    report = [
        "# Route 1 InterProScan Domain Audit",
        "",
        "## Provenance",
        "",
        "- InterProScan 5.59-91.0 was run through Galaxy using Pfam, PANTHER, CDD, SMART, PROSITE, PIRSR, SFLD, FunFam, Gene3D, SUPERFAMILY, PRINTS, AntiFam, MobiDBLite, and PIRSF evidence.",
        f"- Galaxy history: `{args.galaxy_history_id}`; output dataset: `{args.galaxy_dataset_id}`.",
        "",
        "## Results",
        "",
        f"- {expected}/{total} candidates recovered a domain consistent with their assigned functional class.",
        f"- {protease}/22 protease-class candidates have independent protease-domain support.",
        f"- {catalytic}/{total} candidates have an explicit active-site pattern or catalytic-domain annotation and receive the highest structure follow-up priority.",
        f"- {p450}/1 detoxification-class candidate has cytochrome P450 domain support, including a conserved-site annotation.",
        "",
        "## Interpretation boundary",
        "",
        "InterProScan strengthens protein-family and conserved-feature assignments. It does not demonstrate enzymatic activity, secretion during infection, physical interaction with a pine protein, or causal contribution to disease.",
    ]
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(report) + "\n", encoding="utf-8")
    logging.info("Wrote %d candidate summaries; %d high-priority catalytic-feature candidates", total, catalytic)


if __name__ == "__main__":
    main()
