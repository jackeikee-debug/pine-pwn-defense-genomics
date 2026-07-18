#!/usr/bin/env python3
"""Nominate compartment-compatible Route 2 substrates for redirected M8 effectors."""

from __future__ import annotations

import argparse
import csv
import logging
import math
import re
from pathlib import Path


FIELDS = [
    "effector_id", "protease_family", "catalytic_mode", "host_orthogroup_id", "host_symbols",
    "host_family", "host_compartment", "module_ids", "target_rationale", "localization_source",
    "evidence_tier", "pmas_simes_padj", "significant_member_count", "max_abs_interaction_lfc",
    "concordant_species_count", "concordant_species", "discordant_species_count", "discordant_species",
    "ppi_supported_symbol_count", "max_weighted_degree", "target_priority_score", "target_priority_tier",
    "nomination_basis", "direct_interaction_evidence", "interaction_confidence", "recommended_follow_up",
    "claim_ceiling",
]

CLAIM = (
    "Compartment-compatible M8-substrate hypothesis only; homolog localization, expression, and network "
    "support do not demonstrate pine colocalization, binding, proteolysis, or causal hydraulic failure"
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def number(value: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return 0.0
    return parsed if math.isfinite(parsed) else 0.0


def match_target_family(symbols: str, families: list[dict[str, str]]) -> dict[str, str] | None:
    symbol_list = [symbol.strip() for symbol in symbols.replace(",", ";").split(";") if symbol.strip()]
    matches = [family for family in families if any(re.search(family["symbol_regex"], symbol) for symbol in symbol_list)]
    if not matches:
        return None
    if len(matches) > 1:
        raise ValueError(f"Multiple apoplastic families match '{symbols}': {', '.join(row['family_id'] for row in matches)}")
    return matches[0]


def score_candidate(row: dict[str, str]) -> int:
    score = 3 if row.get("evidence_tier") == "Tier A" else (2 if row.get("evidence_tier") == "Tier B" else 1)
    score += 2 if 0 < number(row.get("pmas_simes_padj", "")) < 0.05 else 0
    concordant = int(number(row.get("concordant_species_count", "")))
    score += 2 if concordant >= 2 else (1 if concordant == 1 else 0)
    score += 1 if number(row.get("ppi_supported_symbol_count", "")) > 0 else 0
    score += 2  # Passed the curated apoplastic-family filter.
    score -= 1 if number(row.get("discordant_species_count", "")) > 0 else 0
    return score


def build_hypotheses(
    compatibility_audit_path: Path,
    candidates_path: Path,
    families_path: Path,
    output_path: Path,
    shortlist_path: Path,
    report_path: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    redirected: dict[str, dict[str, str]] = {}
    for row in read_tsv(compatibility_audit_path):
        if row.get("recommended_route") == "redirect_to_apoplastic_apoplastic_cell_wall_targets":
            redirected[row["effector_id"]] = row
    families = read_tsv(families_path)
    rows: list[dict[str, str]] = []

    for candidate in read_tsv(candidates_path):
        if "hydraulic_xylem_collapse" not in candidate.get("mechanism_axes", ""):
            continue
        family = match_target_family(candidate.get("ppi_symbols", ""), families)
        if family is None:
            continue
        score = score_candidate(candidate)
        tier = "A" if score >= 8 else ("B" if score >= 6 else "C")
        for effector_id, effector in sorted(redirected.items()):
            rows.append({
                "effector_id": effector_id, "protease_family": effector.get("protease_family", ""),
                "catalytic_mode": effector.get("catalytic_mode", ""),
                "host_orthogroup_id": candidate["orthogroup_id"], "host_symbols": candidate.get("ppi_symbols", ""),
                "host_family": family["family_id"], "host_compartment": family["host_compartment"],
                "module_ids": candidate.get("module_ids", ""), "target_rationale": family["target_rationale"],
                "localization_source": family["source_url"], "evidence_tier": candidate.get("evidence_tier", ""),
                "pmas_simes_padj": candidate.get("pmas_simes_padj", ""),
                "significant_member_count": candidate.get("significant_member_count", ""),
                "max_abs_interaction_lfc": candidate.get("max_abs_interaction_lfc", ""),
                "concordant_species_count": candidate.get("concordant_species_count", ""),
                "concordant_species": candidate.get("concordant_species", ""),
                "discordant_species_count": candidate.get("discordant_species_count", ""),
                "discordant_species": candidate.get("discordant_species", ""),
                "ppi_supported_symbol_count": candidate.get("ppi_supported_symbol_count", ""),
                "max_weighted_degree": candidate.get("max_weighted_degree", ""),
                "target_priority_score": str(score), "target_priority_tier": tier,
                "nomination_basis": "apoplastic-family localization; infection expression; cross-species response; plant network support",
                "direct_interaction_evidence": "none", "interaction_confidence": "hypothesis_only",
                "recommended_follow_up": "verify pine protein secretion and intact M8 catalytic motif, then perform targeted cleavage assay",
                "claim_ceiling": CLAIM,
            })
    rows.sort(key=lambda row: (-int(row["target_priority_score"]), row["host_orthogroup_id"], row["effector_id"]))
    shortlist = [row for row in rows if row["target_priority_tier"] == "A"]
    write_tsv(output_path, rows)
    write_tsv(shortlist_path, shortlist)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        "# Route 2 apoplastic target hypotheses\n\n"
        f"- Redirected endopeptidase effectors: {len(redirected)}\n"
        f"- Compartment-compatible target hypotheses: {len(rows)}\n"
        f"- Tier A follow-up hypotheses: {len(shortlist)}\n\n"
        "Candidates are selected because their Arabidopsis homolog families are secreted, apoplastic, or "
        "cell-wall associated and because the corresponding pine orthogroups carry infection-expression or "
        "cross-species evidence. These are target nominations, not demonstrated M8 substrates.\n",
        encoding="utf-8",
    )
    logging.info("Wrote %d Route 2 hypotheses and %d Tier A candidates", len(rows), len(shortlist))
    return rows, shortlist


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compatibility-audit", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--families", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--shortlist", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    build_hypotheses(
        args.compatibility_audit, args.candidates, args.families,
        args.output, args.shortlist, args.report,
    )


if __name__ == "__main__":
    main()
