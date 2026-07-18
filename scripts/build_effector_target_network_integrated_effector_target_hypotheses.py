#!/usr/bin/env python3
"""Integrate Route 1 effector evidence and nominate class-compatible host targets."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


EFFECTOR_FIELDS = [
    "effector_id", "effector_class", "signal_peptide_score", "signal_peptide_supported",
    "no_predicted_transmembrane_helix", "expected_domain_supported", "protease_domain_supported",
    "catalytic_feature_annotation", "interpro_terms", "catalytic_feature_terms", "ranking_score",
    "ptm", "mean_ca_plddt", "active_site_mean_plddt", "global_quality_tier", "best_cif_path",
    "effector_evidence_score", "effector_evidence_tier", "evidence_summary", "claim_ceiling",
]

PAIR_FIELDS = [
    "effector_id", "effector_class", "effector_evidence_tier", "host_orthogroup_id",
    "host_symbols", "host_module", "host_network_role", "effector_target_network_role", "top_hub_symbol",
    "ppi_supported_symbol_count", "max_partner_count", "max_weighted_degree",
    "max_betweenness_centrality", "expression_evidence_profile", "mechanism_hypothesis_status",
    "host_network_score", "host_expression_score", "pair_priority_score", "pair_priority_tier",
    "nomination_basis", "direct_interaction_evidence", "interaction_confidence",
    "recommended_follow_up", "claim_ceiling",
]

CLAIM = (
    "Class-compatible effector-host mechanism nomination only; no direct physical interaction, "
    "pine target validation, proteolytic cleavage, or causal resistance evidence"
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def truthy(value: str) -> bool:
    return value.strip().lower() in {"true", "yes", "1"}


def number(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def effector_score(secretion: dict[str, str], domain: dict[str, str], structure: dict[str, str]) -> int:
    score = 0
    score += 2 if truthy(secretion.get("signal_peptide_supported", "")) else 0
    score += 1 if truthy(secretion.get("no_predicted_transmembrane_helix", "")) else 0
    score += 1 if truthy(domain.get("expected_domain_supported", "")) else 0
    score += 1 if truthy(domain.get("catalytic_feature_annotation", "")) else 0
    quality = structure.get("global_quality_tier", "")
    score += 2 if quality == "high_global_confidence" else (1 if quality == "moderate_global_confidence" else 0)
    score += 1 if number(structure.get("active_site_mean_plddt", "")) >= 90 else 0
    return score


def effector_tier(score: int) -> str:
    return "A" if score >= 8 else ("B" if score >= 6 else "C")


def host_scores(row: dict[str, str]) -> tuple[int, int]:
    network = 0
    network += 1 if int(number(row.get("ppi_supported_symbol_count", ""))) > 0 else 0
    network += 1 if number(row.get("max_weighted_degree", "")) >= 10 else 0
    network += 1 if number(row.get("max_betweenness_centrality", "")) > 0 else 0
    fdr_fields = [
        "pden_fdr_supported_transcript_count", "pmas_fdr_supported_feature_count",
        "pmas_interaction_fdr_supported_feature_count", "pthun_two_week_significant_gene_count",
        "pthun_four_week_significant_gene_count",
    ]
    expression = 2 if any(number(row.get(field, "")) > 0 for field in fdr_fields) else 0
    if expression == 0 and "family_level_literature" in row.get("expression_evidence_profile", ""):
        expression = 1
    return network, expression


def pair_tier(score: int, effector_level: str, network: int, expression: int) -> str:
    if score >= 11 and effector_level == "A" and network >= 1 and expression >= 2:
        return "A"
    if score >= 8 and network >= 1:
        return "B"
    return "C"


def build_tables(
    secretion_path: Path,
    interpro_path: Path,
    structures_path: Path,
    mechanisms_path: Path,
    effector_output: Path,
    pair_output: Path,
    report_output: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    secretion = {row["representative_id"]: row for row in read_tsv(secretion_path)}
    domains = {row["representative_id"]: row for row in read_tsv(interpro_path)}
    structures = {row["protein_id"]: row for row in read_tsv(structures_path)}
    missing = sorted(set(structures) - set(secretion) | (set(structures) - set(domains)))
    if missing:
        raise ValueError(f"Structure candidates missing upstream evidence: {', '.join(missing)}")

    effector_rows: list[dict[str, str]] = []
    effector_by_id: dict[str, dict[str, str]] = {}
    for protein_id in sorted(structures):
        sec, dom, struct = secretion[protein_id], domains[protein_id], structures[protein_id]
        score = effector_score(sec, dom, struct)
        tier = effector_tier(score)
        row = {
            "effector_id": protein_id,
            "effector_class": dom.get("effector_classes", ""),
            "signal_peptide_score": sec.get("deepsig_signal_peptide_score", ""),
            "signal_peptide_supported": sec.get("signal_peptide_supported", ""),
            "no_predicted_transmembrane_helix": sec.get("no_predicted_transmembrane_helix", ""),
            "expected_domain_supported": dom.get("expected_domain_supported", ""),
            "protease_domain_supported": dom.get("protease_domain_supported", ""),
            "catalytic_feature_annotation": dom.get("catalytic_feature_annotation", ""),
            "interpro_terms": dom.get("relevant_interpro_terms", ""),
            "catalytic_feature_terms": dom.get("catalytic_feature_terms", ""),
            "ranking_score": struct.get("ranking_score", ""),
            "ptm": struct.get("ptm", ""),
            "mean_ca_plddt": struct.get("mean_ca_plddt", ""),
            "active_site_mean_plddt": struct.get("active_site_mean_plddt", ""),
            "global_quality_tier": struct.get("global_quality_tier", ""),
            "best_cif_path": struct.get("best_cif_path", ""),
            "effector_evidence_score": str(score),
            "effector_evidence_tier": tier,
            "evidence_summary": "sequence secretion prediction; domain annotation; monomer structure prediction",
            "claim_ceiling": "Candidate secreted effector with structural support; activity and host targeting unvalidated",
        }
        effector_rows.append(row)
        effector_by_id[protein_id] = row

    pair_rows: list[dict[str, str]] = []
    mechanisms = sorted(read_tsv(mechanisms_path), key=lambda row: (row["orthogroup_id"], row["effector_class"]))
    for effector in effector_rows:
        classes = {value.strip() for value in effector["effector_class"].replace(",", ";").split(";") if value.strip()}
        for host in mechanisms:
            if host.get("effector_class", "") not in classes:
                continue
            network, expression = host_scores(host)
            score = int(effector["effector_evidence_score"]) + network + expression
            tier = pair_tier(score, effector["effector_evidence_tier"], network, expression)
            pair_rows.append({
                "effector_id": effector["effector_id"], "effector_class": host["effector_class"],
                "effector_evidence_tier": effector["effector_evidence_tier"],
                "host_orthogroup_id": host["orthogroup_id"], "host_symbols": host.get("host_symbols", ""),
                "host_module": host.get("host_module", ""), "host_network_role": host.get("host_network_role", ""),
                "effector_target_network_role": host.get("effector_target_network_role", ""), "top_hub_symbol": host.get("top_hub_symbol", ""),
                "ppi_supported_symbol_count": host.get("ppi_supported_symbol_count", ""),
                "max_partner_count": host.get("max_partner_count", ""),
                "max_weighted_degree": host.get("max_weighted_degree", ""),
                "max_betweenness_centrality": host.get("max_betweenness_centrality", ""),
                "expression_evidence_profile": host.get("expression_evidence_profile", ""),
                "mechanism_hypothesis_status": host.get("mechanism_hypothesis_status", ""),
                "host_network_score": str(network), "host_expression_score": str(expression),
                "pair_priority_score": str(score), "pair_priority_tier": tier,
                "nomination_basis": "effector-class compatibility; host interolog centrality; infection-expression overlay",
                "direct_interaction_evidence": "none", "interaction_confidence": "hypothesis_only",
                "recommended_follow_up": "verify colocalization and cleavage motif; then test a small number of pairs experimentally",
                "claim_ceiling": CLAIM,
            })
    pair_rows.sort(key=lambda row: (-int(row["pair_priority_score"]), row["effector_id"], row["host_orthogroup_id"]))

    write_tsv(effector_output, EFFECTOR_FIELDS, effector_rows)
    write_tsv(pair_output, PAIR_FIELDS, pair_rows)
    report_output.parent.mkdir(parents=True, exist_ok=True)
    counts = {tier: sum(row["pair_priority_tier"] == tier for row in pair_rows) for tier in "ABC"}
    report_output.write_text(
        "# Route 1 integrated effector-target hypotheses\n\n"
        f"- Structure-prioritized effectors: {len(effector_rows)}\n"
        f"- Class-compatible host nominations: {len(pair_rows)}\n"
        f"- Pair priority tiers: A={counts['A']}, B={counts['B']}, C={counts['C']}\n\n"
        "Priority tiers rank follow-up value; they are not interaction-confidence levels. Every pair remains "
        "a class-compatible mechanism hypothesis without direct physical-interaction evidence.\n\n"
        "The current structure-prioritized set contains proteases. Accordingly, pair generation excludes "
        "detoxification-class links such as the earlier SOD module unless a matching effector class is present.\n",
        encoding="utf-8",
    )
    logging.info("Wrote %d effectors and %d target hypotheses", len(effector_rows), len(pair_rows))
    return effector_rows, pair_rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--secretion", type=Path, required=True)
    parser.add_argument("--interpro", type=Path, required=True)
    parser.add_argument("--structures", type=Path, required=True)
    parser.add_argument("--mechanisms", type=Path, required=True)
    parser.add_argument("--effector-output", type=Path, required=True)
    parser.add_argument("--pair-output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    build_tables(
        args.secretion, args.interpro, args.structures, args.mechanisms,
        args.effector_output, args.pair_output, args.report,
    )


if __name__ == "__main__":
    main()
