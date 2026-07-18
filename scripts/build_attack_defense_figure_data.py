#!/usr/bin/env python3
"""Build validated plotting tables for the PWN-pine attack-defense figure."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path


CONTRASTS = (
    ("resistant_pwn_vs_water", "R"),
    ("susceptible_pwn_vs_water", "S"),
    ("genotype_by_inoculum", "GxI"),
)
RESPONSE_NODES = (
    ("r_response_greater", "R-response greater"),
    ("s_response_greater", "S-response greater"),
    ("interaction_unresolved", "unresolved / padj unavailable"),
)
SUPPORTED_EDGE_TYPES = {"functional_prior", "candidate_membership"}
BUBBLE_FIELDS = [
    "display_order",
    "orthogroup_id",
    "candidate_label",
    "contrast_id",
    "contrast_label",
    "log2_fold_change",
    "whole_transcriptome_padj",
    "whole_transcriptome_state",
    "candidate_set_simes_padj",
    "priority_class",
    "cross_species_state",
    "network_support",
    "structure_support",
    "sequence_support",
]
NODE_FIELDS = [
    "node_id",
    "node_label",
    "node_layer",
    "display_order",
    "color_group",
    "evidence_boundary",
]
EDGE_FIELDS = [
    "source_id",
    "target_id",
    "source_layer",
    "target_layer",
    "weight",
    "edge_class",
    "color_group",
    "statistical_state",
    "evidence_boundary",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def finite_float(value: object) -> float | None:
    try:
        number = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def split_values(value: object) -> list[str]:
    return [item.strip() for item in str(value or "").split(";") if item.strip()]


def display_order(row: dict[str, str], identifier: str) -> int:
    try:
        return int(row["display_order"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"Invalid display_order for {identifier}") from exc


def candidate_label(candidate: dict[str, str]) -> str:
    orthogroup_id = candidate.get("orthogroup_id", "")
    if not orthogroup_id:
        raise ValueError("Candidate orthogroup_id must be non-empty")
    symbols = split_values(candidate.get("mapped_symbols", ""))
    return f"{orthogroup_id} | {symbols[0]}" if symbols else orthogroup_id


def significance_state(candidate: dict[str, str], contrast_id: str) -> str:
    state = candidate.get(f"{contrast_id}_whole_transcriptome_state", "")
    if not state:
        state = whole_transcriptome_state(
            candidate.get(f"{contrast_id}_whole_transcriptome_padj", "")
        )
    if state not in {"significant", "not_significant", "unavailable"}:
        raise ValueError(f"Invalid {contrast_id} state for {candidate.get('orthogroup_id', '')}: {state}")
    return state


def validated_fdr(value: object, identifier: str, contrast_id: str, state: str) -> float | None:
    value_text = str(value or "").strip()
    if state == "unavailable" and not value_text:
        return None
    fdr = finite_float(value)
    if fdr is None:
        raise ValueError(f"A finite FDR is required for {contrast_id}: {identifier}")
    if not 0 <= fdr <= 1:
        raise ValueError(f"FDR outside [0, 1] for {contrast_id}: {identifier}")
    return fdr


def classify_response(candidate: dict[str, str]) -> str:
    lfc = finite_float(candidate.get("genotype_by_inoculum_representative_log2_fold_change"))
    state = candidate.get("genotype_by_inoculum_whole_transcriptome_state", "unavailable")
    fdr = validated_fdr(
        candidate.get("genotype_by_inoculum_whole_transcriptome_padj"),
        candidate.get("orthogroup_id", ""),
        "genotype_by_inoculum",
        state,
    )
    if lfc is None or fdr is None or state == "unavailable" or lfc == 0:
        return "interaction_unresolved"
    return "s_response_greater" if lfc > 0 else "r_response_greater"


def whole_transcriptome_state(value: object) -> str:
    value_text = str(value or "").strip()
    if not value_text:
        return "unavailable"
    padj = finite_float(value)
    if padj is None:
        raise ValueError(f"Whole-transcriptome padj must be finite or blank: {value_text}")
    return "significant" if padj < 0.05 else "not_significant"


def attach_whole_transcriptome_statistics(
    candidates: list[dict[str, str]],
    member_audit: list[dict[str, str]],
) -> list[dict[str, str]]:
    selected = {
        (row.get("orthogroup_id", ""), row.get("contrast", "")): row
        for row in member_audit
        if row.get("species") == "Pinus massoniana" and row.get("is_displayed_member") == "yes"
    }
    enriched = []
    for candidate in candidates:
        row = dict(candidate)
        for contrast_id, _ in CONTRASTS:
            audit = selected.get((row.get("orthogroup_id", ""), contrast_id))
            if audit is None:
                raise ValueError(
                    f"Missing displayed P. massoniana member audit for {row.get('orthogroup_id', '')}: {contrast_id}"
                )
            padj = audit.get("whole_transcriptome_padj", "")
            row[f"{contrast_id}_whole_transcriptome_padj"] = padj
            row[f"{contrast_id}_whole_transcriptome_state"] = whole_transcriptome_state(padj)
        enriched.append(row)
    return enriched


def build_bubble_rows(candidates: list[dict[str, str]]) -> list[dict[str, str]]:
    validate_candidates(candidates)
    rows: list[dict[str, str]] = []
    for candidate in sorted(candidates, key=lambda row: display_order(row, row.get("orthogroup_id", ""))):
        for contrast_id, contrast_label in CONTRASTS:
            rows.append(
                {
                    "display_order": candidate["display_order"],
                    "orthogroup_id": candidate["orthogroup_id"],
                    "candidate_label": candidate_label(candidate),
                    "contrast_id": contrast_id,
                    "contrast_label": contrast_label,
                    "log2_fold_change": candidate.get(f"{contrast_id}_representative_log2_fold_change", ""),
                    "whole_transcriptome_padj": candidate.get(
                        f"{contrast_id}_whole_transcriptome_padj", ""
                    ),
                    "whole_transcriptome_state": significance_state(candidate, contrast_id),
                    "candidate_set_simes_padj": candidate.get(f"{contrast_id}_simes_padj", ""),
                    "priority_class": candidate.get("priority_class", ""),
                    "cross_species_state": candidate.get("cross_species_state", ""),
                    "network_support": candidate.get("network_support", ""),
                    "structure_support": candidate.get("structure_support", ""),
                    "sequence_support": candidate.get("sequence_support", ""),
                }
            )
    return rows


def reject_prohibited_boundaries(rows: list[dict[str, str]], boundary_field: str = "evidence_boundary") -> None:
    for row in rows:
        boundary = row.get(boundary_field, "").lower()
        if "direct_interaction" in boundary or "validated_binding" in boundary:
            raise ValueError(f"Unsupported evidence boundary: {row.get(boundary_field, '')}")


def validate_candidates(candidates: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    candidate_by_id = {row.get("orthogroup_id", ""): row for row in candidates}
    if "" in candidate_by_id or len(candidate_by_id) != len(candidates):
        raise ValueError("Candidate orthogroup IDs must be non-empty and unique")
    reject_prohibited_boundaries(candidates, boundary_field="claim_boundary")
    display_orders: set[int] = set()
    for candidate_id, candidate in candidate_by_id.items():
        order = display_order(candidate, candidate_id)
        if order in display_orders:
            raise ValueError(f"Candidate display_order values must be unique: {order}")
        display_orders.add(order)
        for contrast_id, _ in CONTRASTS:
            state = significance_state(candidate, contrast_id)
            validated_fdr(
                candidate.get(f"{contrast_id}_whole_transcriptome_padj"),
                candidate_id,
                contrast_id,
                state,
            )
    return candidate_by_id


def sector_index(sectors: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    indexed: dict[str, dict[str, str]] = {}
    for sector in sectors:
        sector_id = sector.get("sector_id", "")
        if not sector_id or sector_id in indexed:
            raise ValueError(f"Sector IDs must be non-empty and unique: {sector_id}")
        if sector.get("direct_interaction_evidence", "none") != "none":
            raise ValueError(f"Sector reports direct interaction evidence: {sector_id}")
        display_order(sector, sector_id)
        indexed[sector_id] = sector
    return indexed


def require_sector_type(sectors: dict[str, dict[str, str]], sector_id: str, sector_type: str) -> None:
    sector = sectors.get(sector_id)
    if sector is None or sector.get("sector_type") != sector_type:
        raise ValueError(f"Expected {sector_type} sector: {sector_id}")


def node_row(
    node_id: str,
    node_label: str,
    node_layer: str,
    display_order_value: int,
    color_group: str,
    evidence_boundary: str,
) -> dict[str, str]:
    return {
        "node_id": node_id,
        "node_label": node_label,
        "node_layer": node_layer,
        "display_order": str(display_order_value),
        "color_group": color_group,
        "evidence_boundary": evidence_boundary,
    }


def edge_row(
    source_id: str,
    target_id: str,
    source_layer: str,
    target_layer: str,
    weight: float,
    edge_class: str,
    color_group: str,
    evidence_boundary: str,
    statistical_state: str = "not_applicable",
) -> dict[str, str]:
    return {
        "source_id": source_id,
        "target_id": target_id,
        "source_layer": source_layer,
        "target_layer": target_layer,
        "weight": format(weight, ".12g"),
        "edge_class": edge_class,
        "color_group": color_group,
        "statistical_state": statistical_state,
        "evidence_boundary": evidence_boundary,
    }


def build_alluvial_tables(
    candidates: list[dict[str, str]],
    mechanism_edges: list[dict[str, str]],
    sectors: list[dict[str, str]],
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    candidate_by_id = validate_candidates(candidates)

    sector_by_id = sector_index(sectors)
    reject_prohibited_boundaries(mechanism_edges)
    seen_edges: set[tuple[str, str, str]] = set()
    functional_by_module: dict[str, list[dict[str, str]]] = defaultdict(list)
    membership_by_candidate: dict[str, list[dict[str, str]]] = defaultdict(list)
    for edge in mechanism_edges:
        edge_type = edge.get("edge_type", "")
        source_id = edge.get("source_id", "")
        target_id = edge.get("target_id", "")
        key = (source_id, target_id, edge_type)
        if edge_type not in SUPPORTED_EDGE_TYPES:
            raise ValueError(f"Unsupported edge type: {edge_type}")
        if not source_id or not target_id or key in seen_edges:
            raise ValueError(f"Edges must be non-empty and unique: {key}")
        seen_edges.add(key)
        if edge_type == "functional_prior":
            require_sector_type(sector_by_id, source_id, "effector_class")
            require_sector_type(sector_by_id, target_id, "host_module")
            functional_by_module[target_id].append(edge)
        else:
            require_sector_type(sector_by_id, source_id, "host_module")
            require_sector_type(sector_by_id, target_id, "candidate")
            if target_id not in candidate_by_id:
                raise ValueError(f"Candidate membership references unknown candidate: {target_id}")
            membership_by_candidate[target_id].append(edge)

    for candidate_id, candidate in candidate_by_id.items():
        memberships = membership_by_candidate.get(candidate_id, [])
        if not memberships:
            raise ValueError(f"Candidate has no module membership: {candidate_id}")
        declared_modules = set(split_values(candidate.get("module_ids", "")))
        membership_modules = {edge["source_id"] for edge in memberships}
        if declared_modules != membership_modules:
            raise ValueError(f"Candidate module membership disagrees with frozen candidate table: {candidate_id}")
        for module_id in membership_modules:
            if not functional_by_module.get(module_id):
                raise ValueError(f"Host module has no functional-prior source: {module_id}")

    candidate_mass_by_module: dict[str, float] = defaultdict(float)
    membership_edges: list[dict[str, str]] = []
    for candidate in sorted(candidates, key=lambda row: display_order(row, row["orthogroup_id"])):
        candidate_id = candidate["orthogroup_id"]
        memberships = sorted(
            membership_by_candidate[candidate_id],
            key=lambda edge: display_order(sector_by_id[edge["source_id"]], edge["source_id"]),
        )
        weight = 1.0 / len(memberships)
        for membership in memberships:
            module_id = membership["source_id"]
            candidate_mass_by_module[module_id] += weight
            membership_edges.append(
                edge_row(
                    module_id,
                    candidate_id,
                    "host_module",
                    "candidate",
                    weight,
                    "candidate_membership",
                    module_id,
                    membership["evidence_boundary"],
                )
            )

    functional_edges: list[dict[str, str]] = []
    for module_id in sorted(candidate_mass_by_module, key=lambda node: display_order(sector_by_id[node], node)):
        priors = sorted(
            functional_by_module[module_id],
            key=lambda edge: display_order(sector_by_id[edge["source_id"]], edge["source_id"]),
        )
        weight = candidate_mass_by_module[module_id] / len(priors)
        for prior in priors:
            functional_edges.append(
                edge_row(
                    prior["source_id"],
                    module_id,
                    "effector_class",
                    "host_module",
                    weight,
                    "functional_prior",
                    "functional_prior",
                    prior["evidence_boundary"],
                )
            )

    response_edges = [
        edge_row(
            candidate["orthogroup_id"],
            classify_response(candidate),
            "candidate",
            "response",
            1.0,
            "response_summary",
            classify_response(candidate),
            "transcriptomic_interaction_summary",
            significance_state(candidate, "genotype_by_inoculum"),
        )
        for candidate in sorted(candidates, key=lambda row: display_order(row, row["orthogroup_id"]))
    ]

    functional_sources = {edge["source_id"] for edges in functional_by_module.values() for edge in edges}
    active_modules = set(candidate_mass_by_module)
    nodes: list[dict[str, str]] = []
    for node_id in sorted(functional_sources, key=lambda node: display_order(sector_by_id[node], node)):
        sector = sector_by_id[node_id]
        nodes.append(node_row(node_id, sector["label"], "effector_class", display_order(sector, node_id), "functional_prior", "class_level_functional_context"))
    for node_id in sorted(active_modules, key=lambda node: display_order(sector_by_id[node], node)):
        sector = sector_by_id[node_id]
        nodes.append(node_row(node_id, sector["label"], "host_module", display_order(sector, node_id), node_id, "module_annotation_context"))
    for candidate in sorted(candidates, key=lambda row: display_order(row, row["orthogroup_id"])):
        candidate_id = candidate["orthogroup_id"]
        nodes.append(node_row(candidate_id, candidate_label(candidate), "candidate", display_order(candidate, candidate_id), "candidate", "candidate_evidence_synthesis"))
    for order, (node_id, label) in enumerate(RESPONSE_NODES, start=1):
        nodes.append(node_row(node_id, label, "response", order, node_id, "transcriptomic_interaction_summary"))

    all_edges = functional_edges + membership_edges + response_edges
    validate_alluvial_tables(nodes, all_edges, candidate_by_id)
    return nodes, all_edges


def validate_bubble_rows(rows: list[dict[str, str]], expected_candidate_count: int | None = None) -> None:
    candidate_ids = {row["orthogroup_id"] for row in rows}
    if expected_candidate_count is not None and len(candidate_ids) != expected_candidate_count:
        raise ValueError(f"Expected {expected_candidate_count} unique candidates, found {len(candidate_ids)}")
    expected_rows = len(candidate_ids) * len(CONTRASTS)
    if len(rows) != expected_rows:
        raise ValueError(f"Expected {expected_rows} bubble rows, found {len(rows)}")
    for row in rows:
        state = row["whole_transcriptome_state"]
        if state not in {"significant", "not_significant", "unavailable"}:
            raise ValueError(f"Invalid significance state for {row['orthogroup_id']}: {state}")
        validated_fdr(
            row["whole_transcriptome_padj"], row["orthogroup_id"], row["contrast_id"], state
        )
        candidate_value = str(row["candidate_set_simes_padj"] or "").strip()
        candidate_fdr = finite_float(candidate_value)
        if candidate_value and (candidate_fdr is None or not 0 <= candidate_fdr <= 1):
            raise ValueError(
                f"Candidate-set FDR must be finite for {row['orthogroup_id']}: {row['contrast_id']}"
            )
        if state != "unavailable" and finite_float(row["log2_fold_change"]) is None:
            raise ValueError(f"Available expression state requires finite log2 fold change: {row['orthogroup_id']}")


def validate_alluvial_tables(
    nodes: list[dict[str, str]],
    edges: list[dict[str, str]],
    candidate_by_id: dict[str, dict[str, str]],
) -> None:
    node_ids = [node["node_id"] for node in nodes]
    if len(node_ids) != len(set(node_ids)):
        raise ValueError("Alluvial node IDs must be unique")
    node_id_set = set(node_ids)
    incoming: dict[str, float] = defaultdict(float)
    outgoing: dict[str, float] = defaultdict(float)
    for edge in edges:
        if list(edge) != EDGE_FIELDS:
            raise ValueError("Alluvial edges must use the fixed schema")
        statistical_state = edge["statistical_state"]
        if edge["edge_class"] == "response_summary":
            if statistical_state not in {"significant", "not_significant", "unavailable"}:
                raise ValueError("Response edges require a whole-transcriptome statistical state")
        elif statistical_state != "not_applicable":
            raise ValueError("Non-response edges must use statistical_state=not_applicable")
        if edge["source_id"] not in node_id_set or edge["target_id"] not in node_id_set:
            raise ValueError("Alluvial edge references an unknown node")
        weight = finite_float(edge["weight"])
        if weight is None or weight <= 0:
            raise ValueError("Alluvial edge weights must be positive and finite")
        incoming[edge["target_id"]] += weight
        outgoing[edge["source_id"]] += weight
    for candidate_id in candidate_by_id:
        if not math.isclose(incoming[candidate_id], 1.0, rel_tol=0.0, abs_tol=1e-9):
            raise ValueError(f"Candidate incoming flow is not conserved: {candidate_id}")
        if not math.isclose(outgoing[candidate_id], 1.0, rel_tol=0.0, abs_tol=1e-9):
            raise ValueError(f"Candidate outgoing flow is not conserved: {candidate_id}")
    for node in nodes:
        if node["node_layer"] == "host_module" and not math.isclose(
            incoming[node["node_id"]], outgoing[node["node_id"]], rel_tol=0.0, abs_tol=1e-9
        ):
            raise ValueError(f"Host-module flow is not conserved: {node['node_id']}")


def validate_effector_context(
    effector_context: list[dict[str, str]],
    mechanism_edges: list[dict[str, str]],
) -> None:
    classes = [row.get("effector_class", "") for row in effector_context]
    if not classes or "" in classes or len(classes) != len(set(classes)):
        raise ValueError("Effector context classes must be non-empty and unique")
    functional_sources = {
        edge.get("source_id", "")
        for edge in mechanism_edges
        if edge.get("edge_type") == "functional_prior"
    }
    if set(classes) != functional_sources:
        raise ValueError("Effector context must match functional-prior source classes")
    for row in effector_context:
        boundary = row.get("claim_boundary", "").lower()
        if "direct_interaction" in boundary or "validated_binding" in boundary:
            raise ValueError(f"Unsupported effector-context boundary: {row.get('claim_boundary', '')}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--member-audit", type=Path, required=True)
    parser.add_argument("--mechanism-edges", type=Path, required=True)
    parser.add_argument("--effector-context", type=Path, required=True)
    parser.add_argument("--sectors", type=Path, required=True)
    parser.add_argument("--bubble-output", type=Path, required=True)
    parser.add_argument("--nodes-output", type=Path, required=True)
    parser.add_argument("--edges-output", type=Path, required=True)
    args = parser.parse_args()

    candidates = attach_whole_transcriptome_statistics(
        read_tsv(args.candidates), read_tsv(args.member_audit)
    )
    mechanism_edges = read_tsv(args.mechanism_edges)
    effector_context = read_tsv(args.effector_context)
    sectors = read_tsv(args.sectors)
    bubble_rows = build_bubble_rows(candidates)
    validate_bubble_rows(bubble_rows, expected_candidate_count=12)
    validate_effector_context(effector_context, mechanism_edges)
    nodes, edges = build_alluvial_tables(candidates, mechanism_edges, sectors)
    write_tsv(args.bubble_output, BUBBLE_FIELDS, bubble_rows)
    write_tsv(args.nodes_output, NODE_FIELDS, nodes)
    write_tsv(args.edges_output, EDGE_FIELDS, edges)
    print(
        f"Wrote {len(bubble_rows)} bubble rows, {len(nodes)} nodes, and {len(edges)} normalized edges",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
