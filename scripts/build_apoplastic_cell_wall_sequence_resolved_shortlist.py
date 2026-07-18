#!/usr/bin/env python3
"""Resolve Route 2 target priorities with exact pine sequence evidence."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


FIELDS = [
    "effector_id", "host_orthogroup_id", "host_symbols", "host_family",
    "expression_selected_gene_id", "expression_selected_family_status",
    "expression_selected_localization_status", "prior_target_priority_score",
    "prior_target_priority_tier", "sequence_resolved_tier", "recommended_next_step",
    "direct_interaction_evidence", "claim_ceiling",
]
CLAIM = (
    "Sequence-resolved substrate nomination only; secretion prediction and family annotation do not "
    "demonstrate M8 binding, cleavage, or causal contribution to pine wilt disease"
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_shortlist(
    targets_path: Path,
    summary_path: Path,
    output_path: Path,
    report_path: Path,
) -> list[dict[str, str]]:
    summaries = {row["orthogroup_id"]: row for row in read_tsv(summary_path)}
    rows: list[dict[str, str]] = []
    for target in read_tsv(targets_path):
        summary = summaries[target["host_orthogroup_id"]]
        localization = summary["expression_selected_localization_status"]
        family = summary["expression_selected_family_status"]
        if family == "family_annotation_supported" and localization == "supported_secreted_soluble":
            tier = "Tier 1_sequence_supported"
            next_step = "confirm intact pine protein and M8 catalytic motif; prioritize targeted cleavage assay"
        elif family == "family_annotation_supported" and localization == "n_terminal_signal_or_anchor_ambiguous":
            tier = "Tier 2_localization_ambiguous"
            next_step = "resolve N-terminal localization with independent SignalP or experimental secretion evidence"
        else:
            tier = "Tier 3_sequence_gap"
            next_step = "repair annotation or sequence evidence before interaction follow-up"
        rows.append({
            "effector_id": target["effector_id"], "host_orthogroup_id": target["host_orthogroup_id"],
            "host_symbols": target["host_symbols"], "host_family": target["host_family"],
            "expression_selected_gene_id": summary["expression_selected_gene_id"],
            "expression_selected_family_status": family,
            "expression_selected_localization_status": localization,
            "prior_target_priority_score": target["target_priority_score"],
            "prior_target_priority_tier": target["target_priority_tier"],
            "sequence_resolved_tier": tier, "recommended_next_step": next_step,
            "direct_interaction_evidence": "none", "claim_ceiling": CLAIM,
        })
    order = {"Tier 1_sequence_supported": 1, "Tier 2_localization_ambiguous": 2, "Tier 3_sequence_gap": 3}
    rows.sort(key=lambda row: (order[row["sequence_resolved_tier"]], -int(row["prior_target_priority_score"]), row["host_orthogroup_id"]))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    tier1 = [row for row in rows if row["sequence_resolved_tier"] == "Tier 1_sequence_supported"]
    tier2 = [row for row in rows if row["sequence_resolved_tier"] == "Tier 2_localization_ambiguous"]
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        "# Route 2 sequence-resolved shortlist\n\n"
        f"- Tier 1 sequence-supported hypotheses: {len(tier1)}\n"
        f"- Tier 2 localization-ambiguous hypotheses: {len(tier2)}\n\n"
        "Tier 1 means that the exact expression-selected pine member has compatible family annotation and "
        "DeepSig/TMHMM secretion support. It is a follow-up priority, not direct substrate evidence.\n",
        encoding="utf-8",
    )
    logging.info("Wrote %d sequence-resolved targets", len(rows))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    build_shortlist(args.targets, args.summary, args.output, args.report)


if __name__ == "__main__":
    main()
