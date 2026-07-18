#!/usr/bin/env python3
"""Integrate route 1 effector-class anchors, host hubs, and expression evidence."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


FIELDS = [
    "orthogroup_id", "effector_target_network_role", "host_module", "host_network_role", "effector_class",
    "effector_count", "effector_ids", "host_symbols", "ppi_supported_symbol_count",
    "top_hub_symbol", "max_partner_count", "max_weighted_degree", "max_betweenness_centrality",
    "pden_current_candidate_transcript_count", "pden_fdr_supported_transcript_count",
    "pmas_fdr_supported_feature_count", "pmas_interaction_fdr_supported_feature_count",
    "pthun_two_week_significant_gene_count", "pthun_two_week_direction",
    "pthun_four_week_significant_gene_count", "pthun_four_week_direction",
    "expression_evidence_profile", "mechanism_hypothesis_status", "edge_type", "claim_ceiling",
]
CLAIM = "Effector-class to host-hub mechanism hypothesis only; not a direct physical interaction, validated target, or causal resistance mechanism"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(value: str) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def split_multi(value: str) -> list[str]:
    return [item.strip() for item in value.replace(",", ";").split(";") if item.strip()]


def read_de(path: Path) -> dict[str, dict[str, str]]:
    return {row["gene_id"]: row for row in read_tsv(path)}


def summarize_de(genes: set[str], de: dict[str, dict[str, str]]) -> tuple[int, str]:
    significant = []
    for gene in genes:
        if gene not in de:
            continue
        padj = number(de[gene].get("padj", ""))
        if padj is not None and padj < 0.05:
            significant.append(de[gene])
    signs = {"up" if (number(row.get("log2FoldChange", "")) or 0) > 0 else "down" for row in significant if (number(row.get("log2FoldChange", "")) or 0) != 0}
    direction = "mixed" if len(signs) > 1 else (next(iter(signs)) if signs else "not_significant")
    return len(significant), direction


def build_hypotheses(
    seed_path: Path,
    centrality_path: Path,
    pden_path: Path,
    pmas_path: Path,
    stable_genes_path: Path,
    pthun_two_week_path: Path,
    pthun_four_week_path: Path,
    output_path: Path,
    report_path: Path,
) -> list[dict[str, str]]:
    seeds = read_tsv(seed_path)
    centrality: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(centrality_path):
        centrality[row["orthogroup_id"]].append(row)
    pden: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(pden_path):
        for orthogroup in split_multi(row.get("orthogroup_ids", "")):
            pden[orthogroup].append(row)
    pmas: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(pmas_path):
        pmas[row["orthogroup_id"]].append(row)
    ptae_genes: dict[str, set[str]] = defaultdict(set)
    for row in read_tsv(stable_genes_path):
        if row.get("species_id") == "ptae":
            ptae_genes[row["orthogroup_id"]].add(row["gene_id"].split("|", 1)[-1])
    de_two, de_four = read_de(pthun_two_week_path), read_de(pthun_four_week_path)

    rows = []
    seen: set[tuple[str, str]] = set()
    for seed in seeds:
        key = (seed["orthogroup_id"], seed["effector_class"])
        if key in seen:
            continue
        seen.add(key)
        orthogroup = seed["orthogroup_id"]
        network = centrality.get(orthogroup, [])
        supported = [row for row in network if row.get("ppi_support_status") == "string_neighbors_found"]
        top = max(network, key=lambda row: (number(row.get("weighted_degree_score_sum", "")) or 0, row.get("query_symbol", "")), default={})
        pden_rows = pden.get(orthogroup, [])
        pden_current = [row for row in pden_rows if row.get("mapping_status") == "current_pd_effector_target_network_candidate_hit"]
        pden_fdr = [row for row in pden_current if row.get("expression_evidence_tier") != "no_fdr_expression_support"]
        pmas_rows = pmas.get(orthogroup, [])
        pmas_fdr = [row for row in pmas_rows if row.get("resistant_fdr_status") == "fdr_lt_0.05" or row.get("susceptible_fdr_status") == "fdr_lt_0.05"]
        pmas_interaction = [row for row in pmas_rows if row.get("interaction_fdr_status") == "fdr_lt_0.05"]
        two_count, two_direction = summarize_de(ptae_genes.get(orthogroup, set()), de_two)
        four_count, four_direction = summarize_de(ptae_genes.get(orthogroup, set()), de_four)
        profile = []
        if pden_fdr:
            profile.append("pden_current_candidate_fdr")
        elif pden_current:
            profile.append("pden_current_candidate_no_fdr")
        if pmas_fdr:
            profile.append("pmas_inoculation_fdr")
        if pmas_interaction:
            profile.append("pmas_genotype_by_inoculum_fdr")
        if two_count or four_count:
            profile.append("pthun_surrogate_reference_fdr")
        if seed.get("expression_support_level") == "family_level_literature_support":
            profile.append("family_level_literature")
        exact_or_proxy_expression = bool(pden_fdr or pmas_fdr or pmas_interaction or two_count or four_count)
        if supported and exact_or_proxy_expression:
            status = "network_and_expression_supported_hypothesis"
        elif supported:
            status = "network_supported_expression_gap"
        elif exact_or_proxy_expression:
            status = "expression_supported_network_gap"
        else:
            status = "network_and_expression_gap"
        rows.append({
            "orthogroup_id": orthogroup,
            "effector_target_network_role": seed.get("effector_target_network_role", ""),
            "host_module": seed.get("host_module", ""),
            "host_network_role": seed.get("host_network_role", ""),
            "effector_class": seed.get("effector_class", ""),
            "effector_count": seed.get("effector_count", ""),
            "effector_ids": seed.get("effector_ids", ""),
            "host_symbols": ";".join(sorted({row.get("query_symbol", "") for row in network if row.get("query_symbol")})),
            "ppi_supported_symbol_count": str(len(supported)),
            "top_hub_symbol": top.get("query_symbol", ""),
            "max_partner_count": str(max((int(row.get("partner_count", "0") or 0) for row in network), default=0)),
            "max_weighted_degree": top.get("weighted_degree_score_sum", "0.000"),
            "max_betweenness_centrality": f"{max((number(row.get('betweenness_centrality', '')) or 0 for row in network), default=0):.6f}",
            "pden_current_candidate_transcript_count": str(len(pden_current)),
            "pden_fdr_supported_transcript_count": str(len(pden_fdr)),
            "pmas_fdr_supported_feature_count": str(len(pmas_fdr)),
            "pmas_interaction_fdr_supported_feature_count": str(len(pmas_interaction)),
            "pthun_two_week_significant_gene_count": str(two_count),
            "pthun_two_week_direction": two_direction,
            "pthun_four_week_significant_gene_count": str(four_count),
            "pthun_four_week_direction": four_direction,
            "expression_evidence_profile": ";".join(profile) or "none",
            "mechanism_hypothesis_status": status,
            "edge_type": "effector_class_to_host_hub_functional_hypothesis",
            "claim_ceiling": CLAIM,
        })
    order = {"primary": 0, "supporting": 1}
    rows.sort(key=lambda row: (order.get(row["effector_target_network_role"], 2), row["orthogroup_id"], row["effector_class"]))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    supported_rows = [row for row in rows if row["mechanism_hypothesis_status"] == "network_and_expression_supported_hypothesis"]
    lines = [
        "# Route 1 mechanism hypotheses", "",
        f"- Effector-class to host-hub hypothesis rows: {len(rows)}",
        f"- Rows with predicted network plus expression support: {len(supported_rows)}", "",
        "These are functional hypotheses anchored by effector class, host orthogroup annotation, Arabidopsis STRING neighborhoods, and infection-expression evidence. They are not predicted physical docking pairs or validated effector targets.", "",
        "| effector class | host OG | symbols | status | expression evidence |",
        "|---|---|---|---|---|",
    ]
    for row in rows:
        lines.append(f"| {row['effector_class']} | {row['orthogroup_id']} | {row['host_symbols'] or 'none'} | {row['mechanism_hypothesis_status']} | {row['expression_evidence_profile']} |")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["seed", "centrality", "pden", "pmas", "stable-genes", "pthun-two-week", "pthun-four-week", "output", "report"]:
        parser.add_argument(f"--{name}", type=Path, required=True)
    args = parser.parse_args()
    rows = build_hypotheses(args.seed, args.centrality, args.pden, args.pmas, args.stable_genes, args.pthun_two_week, args.pthun_four_week, args.output, args.report)
    print(f"Wrote {len(rows)} route 1 mechanism hypothesis rows")


if __name__ == "__main__":
    main()
