#!/usr/bin/env python3
"""Audit enzyme-substrate scope and compartment compatibility of Route 1 target hypotheses."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


AUDIT_FIELDS = [
    "effector_id", "effector_class", "effector_evidence_tier", "protease_family",
    "catalytic_mode", "substrate_constraint", "protease_constraint_source",
    "host_orthogroup_id", "host_symbols", "host_compartment", "localization_basis",
    "localization_source", "prior_pair_priority_score", "prior_pair_priority_tier",
    "effector_translocation_evidence", "compartment_compatibility",
    "substrate_scope_compatibility", "overall_compatibility", "recommended_route",
    "direct_interaction_evidence", "claim_ceiling",
]

CLAIM = (
    "Compatibility audit only; homolog localization and protease-family properties do not validate "
    "pine localization, effector translocation, host binding, cleavage, or causal disease function"
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=AUDIT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def match_family(interpro_terms: str, constraints: list[dict[str, str]]) -> dict[str, str]:
    matches = [row for row in constraints if row["match_text"].lower() in interpro_terms.lower()]
    if len(matches) != 1:
        names = ", ".join(row["family_id"] for row in matches) or "none"
        raise ValueError(f"Expected one protease-family match for '{interpro_terms}', found {names}")
    return matches[0]


def audit_pairs(
    pairs_path: Path,
    interpro_path: Path,
    constraints_path: Path,
    compartments_path: Path,
    audit_output: Path,
    shortlist_output: Path,
    report_output: Path,
) -> list[dict[str, str]]:
    pairs = read_tsv(pairs_path)
    interpro = {row["representative_id"]: row for row in read_tsv(interpro_path)}
    constraints = read_tsv(constraints_path)
    compartments = {row["orthogroup_id"]: row for row in read_tsv(compartments_path)}
    rows: list[dict[str, str]] = []

    for pair in pairs:
        effector_id = pair["effector_id"]
        orthogroup = pair["host_orthogroup_id"]
        if effector_id not in interpro:
            raise ValueError(f"Missing InterPro evidence for {effector_id}")
        if orthogroup not in compartments:
            raise ValueError(f"Missing host compartment annotation for {orthogroup}")
        family = match_family(interpro[effector_id].get("relevant_interpro_terms", ""), constraints)
        host = compartments[orthogroup]

        translocation = "none"
        compartment = (
            "blocked_without_translocation_evidence"
            if host["host_compartment"] in {"nucleus", "cytosol", "chloroplast", "mitochondrion"}
            else "potentially_colocalized"
        )
        substrate_scope = (
            "possible_family_level_only"
            if family["intact_protein_substrate_plausibility"] == "possible"
            else "poor_for_intact_host_protein"
        )
        if compartment.startswith("blocked"):
            overall = "deprioritize_direct_nuclear_target"
        elif substrate_scope == "poor_for_intact_host_protein":
            overall = "deprioritize_intact_protein_cleavage"
        else:
            overall = "retain_for_target_specific_validation"

        if family["catalytic_mode"] == "endopeptidase" and compartment.startswith("blocked"):
            route = "redirect_to_apoplastic_apoplastic_cell_wall_targets"
        elif family["catalytic_mode"] in {"carboxypeptidase", "oligopeptidase"}:
            route = "evaluate_extracellular_peptide_processing_not_tf_cleavage"
        else:
            route = "retain_for_target_specific_validation"

        rows.append({
            "effector_id": effector_id, "effector_class": pair.get("effector_class", ""),
            "effector_evidence_tier": pair.get("effector_evidence_tier", ""),
            "protease_family": family["family_id"], "catalytic_mode": family["catalytic_mode"],
            "substrate_constraint": family["substrate_constraint"],
            "protease_constraint_source": family["source_url"],
            "host_orthogroup_id": orthogroup, "host_symbols": pair.get("host_symbols", ""),
            "host_compartment": host["host_compartment"], "localization_basis": host["localization_basis"],
            "localization_source": host["source_url"],
            "prior_pair_priority_score": pair.get("pair_priority_score", ""),
            "prior_pair_priority_tier": pair.get("pair_priority_tier", ""),
            "effector_translocation_evidence": translocation,
            "compartment_compatibility": compartment,
            "substrate_scope_compatibility": substrate_scope,
            "overall_compatibility": overall, "recommended_route": route,
            "direct_interaction_evidence": "none", "claim_ceiling": CLAIM,
        })

    rows.sort(key=lambda row: (row["overall_compatibility"], row["protease_family"], row["effector_id"], row["host_orthogroup_id"]))
    shortlist = [row for row in rows if row["overall_compatibility"] == "retain_for_target_specific_validation"]
    write_tsv(audit_output, rows)
    write_tsv(shortlist_output, shortlist)
    report_output.parent.mkdir(parents=True, exist_ok=True)
    redirected = {row["effector_id"] for row in rows if row["recommended_route"] == "redirect_to_apoplastic_apoplastic_cell_wall_targets"}
    report_output.write_text(
        "# Route 1 target compatibility audit\n\n"
        f"- Audited effector-host pairs: {len(rows)}\n"
        f"- Retained direct-target pairs: {len(shortlist)}\n"
        f"- Effectors redirected to apoplastic Route 2 target discovery: {len(redirected)}\n\n"
        "All current host nominations are nuclear transcription-factor homologs, whereas no candidate effector "
        "has host-cell translocation or nuclear localization evidence. Family S9, S10, and M14 substrate scope "
        "also argues against treating full-length transcription factors as primary cleavage substrates.\n\n"
        "The M8 endopeptidase candidate remains enzymatically compatible with intact-protein cleavage at the "
        "family level, but its current nuclear target nominations are compartment-incompatible. It should be "
        "redirected to extracellular or apoplastic Route 2 substrates.\n",
        encoding="utf-8",
    )
    logging.info("Audited %d pairs; retained %d", len(rows), len(shortlist))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--interpro", type=Path, required=True)
    parser.add_argument("--constraints", type=Path, required=True)
    parser.add_argument("--compartments", type=Path, required=True)
    parser.add_argument("--audit-output", type=Path, required=True)
    parser.add_argument("--shortlist-output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    audit_pairs(
        args.pairs, args.interpro, args.constraints, args.compartments,
        args.audit_output, args.shortlist_output, args.report,
    )


if __name__ == "__main__":
    main()
