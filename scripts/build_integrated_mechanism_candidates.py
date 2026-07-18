#!/usr/bin/env python3
"""Integrate OG, expression, and STRING evidence without a composite score."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


FIELDS = [
    "orthogroup_id", "evidence_tier", "module_ids", "mechanism_axes",
    "pmas_simes_padj", "significant_member_count", "max_abs_interaction_lfc",
    "concordant_species_count", "concordant_species", "discordant_species_count",
    "discordant_species", "mixed_response_species", "ppi_symbols", "ppi_supported_symbol_count",
    "max_partner_count", "max_weighted_degree", "max_betweenness_centrality",
    "evidence_profile", "claim_ceiling",
]
CLAIM = "Integrated candidate evidence only; no direct effector target, pine-specific PPI, causal resistance, coevolution, or regional effect"


def _read(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _num(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def build_integrated_candidates(shortlist_path: Path, validation_path: Path, centrality_path: Path, output_path: Path, report_path: Path) -> list[dict[str, str]]:
    validations: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in _read(validation_path): validations[row["orthogroup_id"]].append(row)
    networks: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in _read(centrality_path): networks[row["orthogroup_id"]].append(row)
    rows = []
    for candidate in _read(shortlist_path):
        if candidate["evidence_tier"] not in {"Tier A", "Tier B"}: continue
        orthogroup = candidate["orthogroup_id"]
        evidence = validations.get(orthogroup, [])
        concordant = sorted({r["species"] for r in evidence if r["evidence_classification"] == "directionally_concordant"})
        discordant = sorted({r["species"] for r in evidence if r["evidence_classification"] == "directionally_discordant"})
        mixed = sorted({r["species"] for r in evidence if r["evidence_classification"] == "mixed_cross_species_response"})
        network = networks.get(orthogroup, [])
        supported = [r for r in network if r["ppi_support_status"] == "string_neighbors_found"]
        profile = []
        if concordant: profile.append("cross_species_concordance")
        if discordant or mixed: profile.append("context_dependent_expression")
        if supported: profile.append("predicted_network_support")
        if not network: profile.append("ppi_annotation_gap")
        rows.append({
            "orthogroup_id": orthogroup, "evidence_tier": candidate["evidence_tier"],
            "module_ids": candidate.get("module_ids", ""), "mechanism_axes": candidate.get("mechanism_axes", ""),
            "pmas_simes_padj": candidate.get("pmas_simes_padj", ""),
            "significant_member_count": candidate.get("significant_member_count", ""),
            "max_abs_interaction_lfc": candidate.get("max_abs_interaction_lfc", ""),
            "concordant_species_count": str(len(concordant)), "concordant_species": ";".join(concordant),
            "discordant_species_count": str(len(discordant)), "discordant_species": ";".join(discordant),
            "mixed_response_species": ";".join(mixed),
            "ppi_symbols": ";".join(sorted({r["query_symbol"] for r in network})),
            "ppi_supported_symbol_count": str(len(supported)),
            "max_partner_count": str(max((_num(r["partner_count"]) for r in network), default=0)).removesuffix(".0"),
            "max_weighted_degree": f"{max((_num(r['weighted_degree_score_sum']) for r in network), default=0):.3f}",
            "max_betweenness_centrality": f"{max((_num(r['betweenness_centrality']) for r in network), default=0):.6g}",
            "evidence_profile": ";".join(profile), "claim_ceiling": CLAIM,
        })
    tier_order = {"Tier A": 0, "Tier B": 1}
    rows.sort(key=lambda r: (tier_order[r["evidence_tier"]], -int(r["concordant_species_count"]), -_num(r["max_betweenness_centrality"]), r["orthogroup_id"]))
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with Path(output_path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n"); writer.writeheader(); writer.writerows(rows)

    apoplastic_cell_wall = [r for r in rows if "hydraulic_xylem_collapse" in r["mechanism_axes"]]
    apoplastic_cell_wall_ppi = [r for r in apoplastic_cell_wall if int(r["ppi_supported_symbol_count"]) > 0]
    apoplastic_cell_wall_cross = [r for r in apoplastic_cell_wall if int(r["concordant_species_count"]) > 0]
    apoplastic_cell_wall_modules = {"hydraulic_xylem", "wound_periderm", "phenylpropanoid_lignin", "ros_detoxification"}
    modules = sorted({m for r in apoplastic_cell_wall for m in r["module_ids"].split(";") if m in apoplastic_cell_wall_modules})
    lines = [
        "# Integrated Mechanism Candidates", "", f"- Tier A/B orthogroups: {len(rows)}",
        f"- Candidates with predicted STRING neighbors: {sum(int(r['ppi_supported_symbol_count']) > 0 for r in rows)}",
        f"- Candidates with concordant evidence in at least one independent pine material: {sum(int(r['concordant_species_count']) > 0 for r in rows)}",
        "", "## Route 2 multi-node support", "",
        f"Route 2 contains {len(apoplastic_cell_wall)} candidate orthogroups across {len(modules)} modules ({'; '.join(modules)}). {len(apoplastic_cell_wall_ppi)} have predicted STRING neighbors and {len(apoplastic_cell_wall_cross)} have directionally concordant expression in at least one independent pine material.",
        "", "This pattern is compatible with distributed cell-wall, wound, lignin, and ROS network disruption rather than dependence on one gene. It remains predicted functional-association and expression evidence, not direct effector targeting or causal hydraulic failure.",
    ]
    Path(report_path).parent.mkdir(parents=True, exist_ok=True); Path(report_path).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["shortlist", "validation", "centrality", "output", "report"]: parser.add_argument(f"--{name}", type=Path, required=True)
    args = parser.parse_args(); rows = build_integrated_candidates(args.shortlist, args.validation, args.centrality, args.output, args.report)
    print(f"Wrote {len(rows)} integrated mechanism candidates")


if __name__ == "__main__": main()
