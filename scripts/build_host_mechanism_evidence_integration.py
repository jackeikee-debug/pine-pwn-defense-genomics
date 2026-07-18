#!/usr/bin/env python3
"""Integrate host expression, network, and effector-anchor mechanism evidence."""

from __future__ import annotations

import argparse
import csv
import logging
from collections import defaultdict
from pathlib import Path


apoplastic_cell_wall_MODULES = {
    "hydraulic_xylem",
    "phenylpropanoid_lignin",
    "ros_detoxification",
    "wound_periderm",
}
PRIORITY_ORDER = {
    "mechanism_leading": 0,
    "mechanism_supported": 1,
    "mechanism_context": 2,
}
CLAIM = (
    "Integrated candidate mechanism evidence only; does not establish direct effector targeting, "
    "pine-specific physical interaction, causal resistance, coevolution, or a geographic effect"
)
PRIORITY_FIELDS = [
    "orthogroup_id", "candidate_origins", "evidence_tier", "axis_roles", "module_ids",
    "mechanism_axes", "pmas_simes_padj", "significant_member_count",
    "max_abs_interaction_lfc", "nonproxy_concordant_species_count",
    "nonproxy_concordant_species", "proxy_concordant_species_count",
    "proxy_concordant_species", "nonproxy_discordant_species_count",
    "nonproxy_discordant_species", "proxy_discordant_species_count",
    "proxy_discordant_species", "mixed_species_count", "mixed_species",
    "unavailable_evidence_count", "mapped_symbols", "mapped_symbol_count",
    "max_partner_count", "max_weighted_degree", "max_betweenness_centrality",
    "experimental_supported_symbol_count", "database_supported_symbol_count",
    "experimental_or_database_supported", "anchored_effector_ids",
    "anchored_effector_classes", "effector_anchor_routes", "anchored_effector_count",
    "structurally_supported_effector_count", "best_effector_evidence_tier",
    "best_effector_ptm", "best_effector_mean_plddt", "effector_target_network_priority_class",
    "apoplastic_cell_wall_priority_class", "overall_priority_class", "claim_ceiling",
]
AUDIT_FIELDS = [
    "orthogroup_id", "evidence_layer", "evidence_source", "species", "contrast",
    "mapping_evidence", "can_support_leading", "evidence_classification",
    "effector_id", "effector_evidence_tier", "source_file", "claim_ceiling",
]
CONVERGENCE_FIELDS = [
    "mechanism_axis", "candidate_orthogroup_count", "leading_orthogroup_count",
    "supported_orthogroup_count", "context_orthogroup_count",
    "qualifying_orthogroup_count", "qualifying_module_count", "qualifying_modules",
    "qualifying_nonproxy_species",
    "nonproxy_concordant_orthogroup_count", "backed_network_orthogroup_count",
    "structural_anchor_orthogroup_count", "multi_node_support",
]


def _split(value: str) -> list[str]:
    return [item.strip() for item in value.split(";") if item.strip()]


def _number(value: object, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _text(value: object) -> str:
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, float):
        return f"{value:.12g}"
    return str(value)


def _write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def candidate_universe(
    shortlist_rows: list[dict[str, str]],
    axis_rows: list[dict[str, str]],
) -> dict[str, dict[str, str]]:
    """Return Tier A/B shortlist candidates unioned with configured mechanism seeds."""
    records: dict[str, dict[str, object]] = {}
    for row in shortlist_rows:
        if (
            row.get("evidence_tier") not in {"Tier A", "Tier B"}
            or row.get("candidate_level", "orthogroup") != "orthogroup"
        ):
            continue
        orthogroup = row["orthogroup_id"]
        records[orthogroup] = {
            **row,
            "candidate_origins_set": {"interaction_shortlist"},
            "axis_roles_set": set(),
            "axis_set": set(_split(row.get("mechanism_axes", ""))),
            "module_set": set(_split(row.get("module_ids", ""))),
        }

    role_columns = {
        "primary_candidate_orthogroups": "primary",
        "supporting_candidate_orthogroups": "supporting",
        "deprioritized_candidate_orthogroups": "deprioritized",
    }
    for axis in axis_rows:
        axis_id = axis.get("mechanism_axis", "")
        for column, role in role_columns.items():
            for orthogroup in _split(axis.get(column, "")):
                record = records.setdefault(
                    orthogroup,
                    {
                        "orthogroup_id": orthogroup,
                        "evidence_tier": "",
                        "candidate_origins_set": set(),
                        "axis_roles_set": set(),
                        "axis_set": set(),
                        "module_set": set(),
                    },
                )
                record["candidate_origins_set"].add("mechanism_axis")
                record["axis_roles_set"].add(f"{axis_id}:{role}")
                if axis_id:
                    record["axis_set"].add(axis_id)
                record["module_set"].update(_split(axis.get("host_modules", "")))

    output: dict[str, dict[str, str]] = {}
    origin_order = {"interaction_shortlist": 0, "mechanism_axis": 1}
    for orthogroup, record in records.items():
        origins = sorted(record.pop("candidate_origins_set"), key=lambda item: origin_order[item])
        roles = sorted(record.pop("axis_roles_set"))
        axes = sorted(record.pop("axis_set"))
        modules = sorted(record.pop("module_set"))
        output[orthogroup] = {
            **{key: str(value) for key, value in record.items()},
            "candidate_origins": ";".join(origins),
            "axis_roles": ";".join(roles),
            "mechanism_axes": ";".join(axes),
            "module_ids": ";".join(modules),
        }
    return output


def summarize_expression(
    candidate_ids: set[str],
    validation_rows: list[dict[str, str]],
) -> dict[str, dict[str, object]]:
    """Separate nonproxy, proxy, discordant, mixed, and unavailable evidence."""
    buckets: dict[str, dict[str, set[str] | int]] = {
        orthogroup: {
            "nonproxy_concordant": set(),
            "proxy_concordant": set(),
            "nonproxy_discordant": set(),
            "proxy_discordant": set(),
            "mixed": set(),
            "unavailable_evidence_count": 0,
        }
        for orthogroup in candidate_ids
    }
    for row in validation_rows:
        orthogroup = row.get("orthogroup_id", "")
        if orthogroup not in buckets:
            continue
        species = row.get("species", "")
        classification = row.get("evidence_classification", "").lower()
        nonproxy = row.get("can_support_tier_a", "").lower() == "yes"
        prefix = "nonproxy" if nonproxy else "proxy"
        if classification == "directionally_concordant" and species:
            buckets[orthogroup][f"{prefix}_concordant"].add(species)
        elif classification == "directionally_discordant" and species:
            buckets[orthogroup][f"{prefix}_discordant"].add(species)
        elif classification.startswith("mixed") and species:
            buckets[orthogroup]["mixed"].add(species)
        elif classification == "unavailable":
            buckets[orthogroup]["unavailable_evidence_count"] += 1

    output: dict[str, dict[str, object]] = {}
    for orthogroup, bucket in buckets.items():
        row: dict[str, object] = {}
        for key in (
            "nonproxy_concordant",
            "proxy_concordant",
            "nonproxy_discordant",
            "proxy_discordant",
            "mixed",
        ):
            species = sorted(bucket[key])
            row[f"{key}_species_count"] = len(species)
            row[f"{key}_species"] = ";".join(species)
        row["unavailable_evidence_count"] = bucket["unavailable_evidence_count"]
        output[orthogroup] = row
    return output


def summarize_network(
    candidate_ids: set[str],
    centrality_rows: list[dict[str, str]],
    interaction_rows: list[dict[str, str]],
) -> dict[str, dict[str, object]]:
    """Aggregate interolog centrality and evidence channels by orthogroup."""
    incident: dict[str, list[dict[str, str]]] = defaultdict(list)
    for edge in interaction_rows:
        for endpoint in (edge.get("preferredName_A", ""), edge.get("preferredName_B", "")):
            if endpoint:
                incident[endpoint].append(edge)

    output = {
        orthogroup: {
            "mapped_symbols": "",
            "mapped_symbol_count": 0,
            "max_partner_count": 0,
            "max_weighted_degree": 0.0,
            "max_betweenness_centrality": 0.0,
            "experimental_supported_symbol_count": 0,
            "database_supported_symbol_count": 0,
            "experimental_or_database_supported": False,
        }
        for orthogroup in candidate_ids
    }
    symbols: dict[str, set[str]] = defaultdict(set)
    experimental: dict[str, set[str]] = defaultdict(set)
    database: dict[str, set[str]] = defaultdict(set)
    for row in centrality_rows:
        orthogroup = row.get("orthogroup_id", "")
        if orthogroup not in output or row.get("mapping_status") != "mapped":
            continue
        symbol = row.get("query_symbol", "")
        if symbol:
            symbols[orthogroup].add(symbol)
        output[orthogroup]["max_partner_count"] = max(
            int(output[orthogroup]["max_partner_count"]),
            int(_number(row.get("partner_count"))),
        )
        output[orthogroup]["max_weighted_degree"] = max(
            float(output[orthogroup]["max_weighted_degree"]),
            _number(row.get("weighted_degree_score_sum")),
        )
        output[orthogroup]["max_betweenness_centrality"] = max(
            float(output[orthogroup]["max_betweenness_centrality"]),
            _number(row.get("betweenness_centrality")),
        )
        for edge in incident.get(symbol, []):
            if _number(edge.get("escore")) > 0:
                experimental[orthogroup].add(symbol)
            if _number(edge.get("dscore")) > 0:
                database[orthogroup].add(symbol)

    for orthogroup in candidate_ids:
        output[orthogroup]["mapped_symbols"] = ";".join(sorted(symbols[orthogroup]))
        output[orthogroup]["mapped_symbol_count"] = len(symbols[orthogroup])
        output[orthogroup]["experimental_supported_symbol_count"] = len(experimental[orthogroup])
        output[orthogroup]["database_supported_symbol_count"] = len(database[orthogroup])
        output[orthogroup]["experimental_or_database_supported"] = bool(
            experimental[orthogroup] or database[orthogroup]
        )
    return output


def summarize_effector_anchors(
    candidate_ids: set[str],
    effector_target_network_rows: list[dict[str, str]],
    apoplastic_cell_wall_rows: list[dict[str, str]],
    effector_rows: list[dict[str, str]],
) -> dict[str, dict[str, object]]:
    """Summarize effector-anchor quality without changing host classification."""
    quality = {row.get("effector_id", ""): row for row in effector_rows}
    anchors: dict[str, dict[str, set[str]]] = {
        orthogroup: {"effectors": set(), "classes": set(), "routes": set()}
        for orthogroup in candidate_ids
    }
    for route, rows in (
        ("effector_target_network", effector_target_network_rows),
        ("hydraulic_xylem_collapse", apoplastic_cell_wall_rows),
    ):
        for row in rows:
            orthogroup = row.get("host_orthogroup_id", "")
            if orthogroup not in anchors:
                continue
            effector = row.get("effector_id", "")
            if effector:
                anchors[orthogroup]["effectors"].add(effector)
            if row.get("effector_class"):
                anchors[orthogroup]["classes"].add(row["effector_class"])
            elif row.get("protease_family"):
                anchors[orthogroup]["classes"].add(row["protease_family"])
            anchors[orthogroup]["routes"].add(route)

    tier_rank = {"A": 0, "B": 1, "C": 2}
    output: dict[str, dict[str, object]] = {}
    for orthogroup, values in anchors.items():
        effector_rows_for_host = [quality.get(effector, {}) for effector in values["effectors"]]
        tiers = [row.get("effector_evidence_tier", "") for row in effector_rows_for_host]
        tiers = [tier for tier in tiers if tier]
        structured = [
            row for row in effector_rows_for_host
            if row.get("best_cif_path") or row.get("ptm") or row.get("mean_ca_plddt")
        ]
        output[orthogroup] = {
            "anchored_effector_ids": ";".join(sorted(values["effectors"])),
            "anchored_effector_classes": ";".join(sorted(values["classes"])),
            "effector_anchor_routes": ";".join(sorted(values["routes"])),
            "anchored_effector_count": len(values["effectors"]),
            "structurally_supported_effector_count": len(structured),
            "best_effector_evidence_tier": min(tiers, key=lambda tier: tier_rank.get(tier, 99)) if tiers else "",
            "best_effector_ptm": max((_number(row.get("ptm")) for row in structured), default=0.0),
            "best_effector_mean_plddt": max(
                (_number(row.get("mean_ca_plddt")) for row in structured), default=0.0
            ),
        }
    return output


def classify_route(
    route: str,
    source: dict[str, object],
    expression: dict[str, object],
    network: dict[str, object],
) -> str:
    """Apply the predeclared route-specific host priority gate."""
    shortlist = source.get("evidence_tier") in {"Tier A", "Tier B"}
    concordant = int(expression.get("nonproxy_concordant_species_count", 0))
    backed_network = bool(network.get("experimental_or_database_supported", False))
    mapped = int(network.get("mapped_symbol_count", 0)) > 0
    if not shortlist or concordant < 1:
        return "mechanism_context"
    if route == "effector_target_network":
        return "mechanism_leading" if mapped and backed_network else "mechanism_supported"
    if route == "hydraulic_xylem_collapse":
        has_module = bool(set(_split(str(source.get("module_ids", "")))) & apoplastic_cell_wall_MODULES)
        if not has_module:
            return "mechanism_context"
        return "mechanism_leading" if backed_network or concordant >= 2 else "mechanism_supported"
    return "mechanism_context"


def summarize_route_convergence(
    candidate_rows: list[dict[str, str]],
) -> dict[str, dict[str, str]]:
    """Summarize independent orthogroup and module support for each mechanism route."""
    output: dict[str, dict[str, str]] = {}
    for route, class_column in (
        ("effector_target_network", "effector_target_network_priority_class"),
        ("hydraulic_xylem_collapse", "apoplastic_cell_wall_priority_class"),
    ):
        route_rows = [row for row in candidate_rows if route in _split(row.get("mechanism_axes", ""))]
        qualifying = [
            row for row in route_rows
            if row.get(class_column) in {"mechanism_leading", "mechanism_supported"}
        ]
        modules = {
            module
            for row in qualifying
            for module in _split(row.get("module_ids", ""))
            if route != "hydraulic_xylem_collapse" or module in apoplastic_cell_wall_MODULES
        }
        orthogroups = {row.get("orthogroup_id", "") for row in qualifying if row.get("orthogroup_id")}
        qualifying_species = {
            species
            for row in qualifying
            for species in _split(row.get("nonproxy_concordant_species", ""))
        }
        output[route] = {
            "mechanism_axis": route,
            "candidate_orthogroup_count": str(len(route_rows)),
            "leading_orthogroup_count": str(sum(row.get(class_column) == "mechanism_leading" for row in route_rows)),
            "supported_orthogroup_count": str(sum(row.get(class_column) == "mechanism_supported" for row in route_rows)),
            "context_orthogroup_count": str(sum(row.get(class_column) == "mechanism_context" for row in route_rows)),
            "qualifying_orthogroup_count": str(len(orthogroups)),
            "qualifying_module_count": str(len(modules)),
            "qualifying_modules": ";".join(sorted(modules)),
            "qualifying_nonproxy_species": ";".join(sorted(qualifying_species)),
            "nonproxy_concordant_orthogroup_count": str(sum(int(row.get("nonproxy_concordant_species_count", "0")) > 0 for row in route_rows)),
            "backed_network_orthogroup_count": str(sum(row.get("experimental_or_database_supported") == "true" for row in route_rows)),
            "structural_anchor_orthogroup_count": str(sum(int(row.get("structurally_supported_effector_count", "0")) > 0 for row in route_rows)),
            "multi_node_support": str(len(orthogroups) >= 2 and len(modules) >= 2).lower(),
        }
    return output


def _expression_audit(
    candidate_ids: set[str], rows: list[dict[str, str]]
) -> list[dict[str, str]]:
    return [
        {
            "orthogroup_id": row.get("orthogroup_id", ""),
            "evidence_layer": "infection_expression",
            "evidence_source": f"{row.get('species', '')}:{row.get('contrast', '')}",
            "species": row.get("species", ""),
            "contrast": row.get("contrast", ""),
            "mapping_evidence": row.get("mapping_evidence", ""),
            "can_support_leading": row.get("can_support_tier_a", ""),
            "evidence_classification": row.get("evidence_classification", ""),
            "effector_id": "",
            "effector_evidence_tier": "",
            "source_file": row.get("source_file", ""),
            "claim_ceiling": CLAIM,
        }
        for row in rows
        if row.get("orthogroup_id") in candidate_ids
    ]


def _network_audit(
    candidate_ids: set[str],
    rows: list[dict[str, str]],
    network: dict[str, dict[str, object]],
    source_file: str,
) -> list[dict[str, str]]:
    output = []
    for row in rows:
        orthogroup = row.get("orthogroup_id", "")
        if orthogroup not in candidate_ids:
            continue
        backed = bool(network[orthogroup]["experimental_or_database_supported"])
        output.append({
            "orthogroup_id": orthogroup,
            "evidence_layer": "interolog_network",
            "evidence_source": row.get("query_symbol", ""),
            "species": "Arabidopsis thaliana",
            "contrast": "",
            "mapping_evidence": "arabidopsis_interolog_proxy",
            "can_support_leading": str(backed).lower(),
            "evidence_classification": (
                "experimental_or_database_channel" if backed else row.get("ppi_support_status", "prediction_only")
            ),
            "effector_id": "",
            "effector_evidence_tier": "",
            "source_file": source_file,
            "claim_ceiling": CLAIM,
        })
    return output


def _anchor_audit(
    candidate_ids: set[str],
    route_rows: list[tuple[str, list[dict[str, str]], str]],
    effector_quality: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    output = []
    for route, rows, source_file in route_rows:
        for row in rows:
            orthogroup = row.get("host_orthogroup_id", "")
            if orthogroup not in candidate_ids:
                continue
            effector = row.get("effector_id", "")
            quality = effector_quality.get(effector, {})
            output.append({
                "orthogroup_id": orthogroup,
                "evidence_layer": "effector_anchor",
                "evidence_source": route,
                "species": "Bursaphelenchus xylophilus",
                "contrast": "",
                "mapping_evidence": "class_or_compartment_compatibility_hypothesis",
                "can_support_leading": "no",
                "evidence_classification": row.get("interaction_confidence", "hypothesis_only"),
                "effector_id": effector,
                "effector_evidence_tier": quality.get("effector_evidence_tier", ""),
                "source_file": source_file,
                "claim_ceiling": CLAIM,
            })
    return output


def build_host_mechanism_integration(
    shortlist_path: Path,
    axes_path: Path,
    validation_path: Path,
    centrality_path: Path,
    interactions_path: Path,
    effector_target_network_path: Path,
    apoplastic_cell_wall_path: Path,
    effectors_path: Path,
    priority_output: Path,
    audit_output: Path,
    convergence_output: Path,
    report_output: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]], str]:
    """Build host priorities, source audit, route convergence, and report."""
    shortlist_rows = _read_tsv(shortlist_path)
    eligible_shortlist = [
        row for row in shortlist_rows
        if row.get("evidence_tier") in {"Tier A", "Tier B"}
        and row.get("candidate_level", "orthogroup") == "orthogroup"
    ]
    duplicate_shortlist = [
        orthogroup for orthogroup in {row.get("orthogroup_id", "") for row in eligible_shortlist}
        if orthogroup and sum(row.get("orthogroup_id") == orthogroup for row in eligible_shortlist) > 1
    ]
    if duplicate_shortlist:
        raise ValueError(f"Duplicate shortlist orthogroups: {','.join(sorted(duplicate_shortlist))}")
    axis_rows = _read_tsv(axes_path)
    validation_rows = _read_tsv(validation_path)
    centrality_rows = _read_tsv(centrality_path)
    interaction_rows = _read_tsv(interactions_path)
    effector_target_network_rows = _read_tsv(effector_target_network_path)
    apoplastic_cell_wall_rows = _read_tsv(apoplastic_cell_wall_path)
    effector_rows = _read_tsv(effectors_path)

    universe = candidate_universe(shortlist_rows, axis_rows)
    candidate_ids = set(universe)
    expression = summarize_expression(candidate_ids, validation_rows)
    network = summarize_network(candidate_ids, centrality_rows, interaction_rows)
    anchors = summarize_effector_anchors(candidate_ids, effector_target_network_rows, apoplastic_cell_wall_rows, effector_rows)

    rows: list[dict[str, str]] = []
    for orthogroup in candidate_ids:
        source = universe[orthogroup]
        axes = set(_split(source.get("mechanism_axes", "")))
        effector_target_network_class = (
            classify_route("effector_target_network", source, expression[orthogroup], network[orthogroup])
            if "effector_target_network" in axes else "not_applicable"
        )
        apoplastic_cell_wall_class = (
            classify_route("hydraulic_xylem_collapse", source, expression[orthogroup], network[orthogroup])
            if "hydraulic_xylem_collapse" in axes else "not_applicable"
        )
        applicable = [value for value in (effector_target_network_class, apoplastic_cell_wall_class) if value in PRIORITY_ORDER]
        overall = min(applicable, key=lambda value: PRIORITY_ORDER[value]) if applicable else "mechanism_context"
        combined = {
            "orthogroup_id": orthogroup,
            **source,
            **expression[orthogroup],
            **network[orthogroup],
            **anchors[orthogroup],
            "effector_target_network_priority_class": effector_target_network_class,
            "apoplastic_cell_wall_priority_class": apoplastic_cell_wall_class,
            "overall_priority_class": overall,
            "claim_ceiling": CLAIM,
        }
        rows.append({field: _text(combined.get(field, "")) for field in PRIORITY_FIELDS})

    rows.sort(key=lambda row: (
        PRIORITY_ORDER.get(row["overall_priority_class"], 9),
        -int(row["nonproxy_concordant_species_count"]),
        row["experimental_or_database_supported"] != "true",
        -_number(row["max_betweenness_centrality"]),
        row["orthogroup_id"],
    ))

    effector_quality = {row.get("effector_id", ""): row for row in effector_rows}
    audit = _expression_audit(candidate_ids, validation_rows)
    audit.extend(_network_audit(candidate_ids, centrality_rows, network, str(centrality_path)))
    audit.extend(_anchor_audit(
        candidate_ids,
        [
            ("effector_target_network", effector_target_network_rows, str(effector_target_network_path)),
            ("hydraulic_xylem_collapse", apoplastic_cell_wall_rows, str(apoplastic_cell_wall_path)),
        ],
        effector_quality,
    ))
    audit.sort(key=lambda row: (row["orthogroup_id"], row["evidence_layer"], row["evidence_source"]))

    convergence_map = summarize_route_convergence(rows)
    convergence = [convergence_map[route] for route in ("effector_target_network", "hydraulic_xylem_collapse")]
    apoplastic_cell_wall_multinode = convergence_map["hydraulic_xylem_collapse"]["multi_node_support"] == "true"
    apoplastic_cell_wall_species = convergence_map["hydraulic_xylem_collapse"]["qualifying_nonproxy_species"]
    decision = "apoplastic_cell_wall_multinode_module_support" if apoplastic_cell_wall_multinode else "apoplastic_cell_wall_single_node_or_single_module_support"

    _write_tsv(priority_output, PRIORITY_FIELDS, rows)
    _write_tsv(audit_output, AUDIT_FIELDS, audit)
    _write_tsv(convergence_output, CONVERGENCE_FIELDS, convergence)
    report_output.parent.mkdir(parents=True, exist_ok=True)
    class_counts = {
        priority: sum(row["overall_priority_class"] == priority for row in rows)
        for priority in PRIORITY_ORDER
    }
    leading = [row["orthogroup_id"] for row in rows if row["overall_priority_class"] == "mechanism_leading"]
    consistently_concordant = [
        row["orthogroup_id"]
        for row in rows
        if row["overall_priority_class"] == "mechanism_leading"
        and int(row["nonproxy_concordant_species_count"]) >= 2
        and int(row["nonproxy_discordant_species_count"]) == 0
        and int(row["mixed_species_count"]) == 0
    ]
    direction_heterogeneous = [
        row["orthogroup_id"]
        for row in rows
        if row["overall_priority_class"] == "mechanism_leading"
        and (
            int(row["nonproxy_discordant_species_count"]) > 0
            or int(row["mixed_species_count"]) > 0
        )
    ]
    report = (
        "# Integrated host mechanism candidates\n\n"
        f"- Candidate orthogroups: {len(rows)}\n"
        f"- Mechanism leading: {class_counts['mechanism_leading']}\n"
        f"- Mechanism supported: {class_counts['mechanism_supported']}\n"
        f"- Mechanism context: {class_counts['mechanism_context']}\n"
        f"- Leading orthogroups: {';'.join(leading) if leading else 'none'}\n"
        f"- Consistently concordant multispecies leading orthogroups: {';'.join(consistently_concordant) if consistently_concordant else 'none'}\n"
        f"- Direction-heterogeneous leading orthogroups: {';'.join(direction_heterogeneous) if direction_heterogeneous else 'none'}\n"
        f"- Route 2 convergence decision: `{decision}`\n\n"
        f"- Qualifying Route 2 nonproxy species: {apoplastic_cell_wall_species if apoplastic_cell_wall_species else 'none'}\n\n"
        "Host priority is gated by infection expression, mapping quality, and route-specific independent "
        "support. Effector monomer structures describe anchor quality and do not promote host candidates.\n\n"
        f"Claim ceiling: {CLAIM}.\n"
    )
    report_output.write_text(report, encoding="utf-8")
    logging.info("Integrated %d host candidates and %d audit rows", len(rows), len(audit))
    return rows, audit, convergence, decision


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shortlist", type=Path, required=True)
    parser.add_argument("--axes", type=Path, required=True)
    parser.add_argument("--validation", type=Path, required=True)
    parser.add_argument("--centrality", type=Path, required=True)
    parser.add_argument("--interactions", type=Path, required=True)
    parser.add_argument("--effector_target_network-hypotheses", type=Path, required=True)
    parser.add_argument("--apoplastic_cell_wall-hypotheses", type=Path, required=True)
    parser.add_argument("--effectors", type=Path, required=True)
    parser.add_argument("--priority-output", type=Path, required=True)
    parser.add_argument("--audit-output", type=Path, required=True)
    parser.add_argument("--convergence-output", type=Path, required=True)
    parser.add_argument("--report-output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    build_host_mechanism_integration(
        shortlist_path=args.shortlist,
        axes_path=args.axes,
        validation_path=args.validation,
        centrality_path=args.centrality,
        interactions_path=args.interactions,
        effector_target_network_path=args.effector_target_network_hypotheses,
        apoplastic_cell_wall_path=args.apoplastic_cell_wall_hypotheses,
        effectors_path=args.effectors,
        priority_output=args.priority_output,
        audit_output=args.audit_output,
        convergence_output=args.convergence_output,
        report_output=args.report_output,
    )


if __name__ == "__main__":
    main()
