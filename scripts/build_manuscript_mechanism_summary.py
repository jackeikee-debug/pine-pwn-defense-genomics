#!/usr/bin/env python3
"""Build strict Route 2 and manuscript-facing mechanism evidence tables."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


LOGGER = logging.getLogger(__name__)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def as_int(value: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def split_values(value: str) -> set[str]:
    return {item.strip() for item in value.split(";") if item.strip()}


def direction_class(priority: dict[str, str]) -> str:
    concordant = as_int(priority.get("nonproxy_concordant_species_count", "0"))
    discordant = as_int(priority.get("nonproxy_discordant_species_count", "0"))
    if concordant and discordant:
        return "direction_heterogeneous"
    if concordant:
        return "concordant_only"
    if discordant:
        return "discordant_only"
    return "no_nonproxy_directional_support"


def sequence_support_class(sequence_row: dict[str, str]) -> str:
    tier = sequence_row.get("sequence_resolved_tier", "")
    if tier == "Tier 1_sequence_supported":
        return "tier1_candidate_substrate"
    if tier:
        return "localization_or_sequence_caution"
    return "not_sequence_resolved"


def build_apoplastic_cell_wall_node_audit(
    hypotheses_path: Path,
    priorities_path: Path,
    sequence_path: Path,
) -> list[dict[str, str]]:
    priorities = {row["orthogroup_id"]: row for row in read_tsv(priorities_path)}
    sequence_rows = {
        row["host_orthogroup_id"]: row for row in read_tsv(sequence_path)
    }
    nodes: dict[str, dict[str, str]] = {}
    for hypothesis in read_tsv(hypotheses_path):
        orthogroup_id = hypothesis["host_orthogroup_id"]
        if orthogroup_id in nodes:
            raise ValueError(f"Duplicate Route 2 orthogroup node: {orthogroup_id}")
        priority = priorities.get(orthogroup_id, {})
        sequence = sequence_rows.get(orthogroup_id, {})
        direction = direction_class(priority)
        overall_priority = priority.get("overall_priority_class", "not_in_integrated_priority")
        manuscript_qualifying = (
            overall_priority in {"mechanism_leading", "mechanism_supported"}
            and as_int(priority.get("nonproxy_concordant_species_count", "0")) > 0
        )
        nodes[orthogroup_id] = {
            "orthogroup_id": orthogroup_id,
            "effector_id": hypothesis.get("effector_id", ""),
            "host_symbols": hypothesis.get("host_symbols", ""),
            "host_family": hypothesis.get("host_family", ""),
            "module_ids": hypothesis.get("module_ids", ""),
            "target_priority_tier": hypothesis.get("target_priority_tier", ""),
            "overall_priority_class": overall_priority,
            "apoplastic_cell_wall_priority_class": priority.get("apoplastic_cell_wall_priority_class", ""),
            "nonproxy_concordant_species_count": priority.get(
                "nonproxy_concordant_species_count", "0"
            ),
            "nonproxy_concordant_species": priority.get("nonproxy_concordant_species", ""),
            "nonproxy_discordant_species_count": priority.get(
                "nonproxy_discordant_species_count", "0"
            ),
            "nonproxy_discordant_species": priority.get("nonproxy_discordant_species", ""),
            "nonproxy_direction_class": direction,
            "network_supported": priority.get("experimental_or_database_supported", "false"),
            "sequence_resolved_tier": sequence.get("sequence_resolved_tier", ""),
            "sequence_support_class": sequence_support_class(sequence),
            "manuscript_module_qualifying": "yes" if manuscript_qualifying else "no",
            "direct_cleavage_evidence": "none",
            "claim_ceiling": hypothesis.get("claim_ceiling", ""),
        }
    return [nodes[key] for key in sorted(nodes)]


def summarize_apoplastic_cell_wall_nodes(rows: list[dict[str, str]]) -> dict[str, str]:
    qualifying = [row for row in rows if row["manuscript_module_qualifying"] == "yes"]
    families = {row["host_family"] for row in qualifying if row["host_family"]}
    modules: set[str] = set()
    species: set[str] = set()
    for row in qualifying:
        modules.update(split_values(row["module_ids"]))
        species.update(split_values(row["nonproxy_concordant_species"]))
    strict = (
        len(qualifying) >= 3
        and len(families) >= 2
        and len(modules) >= 2
        and len(species) >= 2
    )
    peroxidases = [row for row in rows if row["host_family"] == "PER"]
    peroxidase_classes = {row["nonproxy_direction_class"] for row in peroxidases}
    if "concordant_only" in peroxidase_classes:
        peroxidase_status = "independent_concordant_support"
    elif "direction_heterogeneous" in peroxidase_classes:
        peroxidase_status = "direction_heterogeneous_support"
    elif "discordant_only" in peroxidase_classes:
        peroxidase_status = "discordant_only_evidence"
    else:
        peroxidase_status = "no_nonproxy_directional_support"
    return {
        "route": "hydraulic_xylem_collapse",
        "audited_independent_orthogroup_count": str(len(rows)),
        "qualifying_independent_orthogroup_count": str(len(qualifying)),
        "qualifying_host_family_count": str(len(families)),
        "qualifying_host_families": ";".join(sorted(families)),
        "qualifying_module_count": str(len(modules)),
        "qualifying_modules": ";".join(sorted(modules)),
        "qualifying_nonproxy_species_count": str(len(species)),
        "qualifying_nonproxy_species": ";".join(sorted(species)),
        "strict_module_multinode_support": "yes" if strict else "no",
        "tier1_candidate_substrate_count": str(
            sum(row["sequence_support_class"] == "tier1_candidate_substrate" for row in rows)
        ),
        "tier1_candidate_substrate_orthogroups": ";".join(
            row["orthogroup_id"]
            for row in rows
            if row["sequence_support_class"] == "tier1_candidate_substrate"
        ),
        "peroxidase_orthogroup_count": str(len(peroxidases)),
        "peroxidase_branch_status": peroxidase_status,
        "direct_cleavage_evidence": "none",
        "claim_ceiling": (
            "Multinode module evidence only; does not establish M8 binding, cleavage, "
            "causal hydraulic failure, or resistance."
        ),
    }


def og0005853_case_study_gate(audit_path: Path, expression_path: Path) -> bool:
    audit_rows = [
        row for row in read_tsv(audit_path) if row.get("orthogroup_id") == "OG0005853"
    ]
    if len(audit_rows) != 1:
        return False
    audit = audit_rows[0]
    complete = all(
        audit.get(field) == "yes"
        for field in (
            "annotation_coherent",
            "tree_membership_complete",
            "sequence_membership_complete",
        )
    )
    up_species = {
        row.get("species_id", "")
        for row in read_tsv(expression_path)
        if row.get("orthogroup_id") == "OG0005853"
        and row.get("de_interpretation") == "significant_up"
    }
    return complete and {"pden", "pstrobus"} <= up_species


def build_manuscript_candidates(
    priorities_path: Path,
    apoplastic_cell_wall_node_rows: list[dict[str, str]],
    og0005853_audit_path: Path,
    og0005853_expression_path: Path,
) -> list[dict[str, str]]:
    apoplastic_cell_wall = {row["orthogroup_id"]: row for row in apoplastic_cell_wall_node_rows}
    case_study_passes = og0005853_case_study_gate(
        og0005853_audit_path, og0005853_expression_path
    )
    output: list[dict[str, str]] = []
    for row in read_tsv(priorities_path):
        overall = row.get("overall_priority_class", "")
        if overall not in {"mechanism_leading", "mechanism_supported"}:
            continue
        orthogroup_id = row["orthogroup_id"]
        node = apoplastic_cell_wall.get(orthogroup_id, {})
        is_case_study = orthogroup_id == "OG0005853" and case_study_passes
        output.append(
            {
                "orthogroup_id": orthogroup_id,
                "manuscript_role": "cross_species_case_study" if is_case_study else overall,
                "overall_priority_class": overall,
                "mechanism_axes": row.get("mechanism_axes", ""),
                "module_ids": row.get("module_ids", ""),
                "nonproxy_concordant_species_count": row.get(
                    "nonproxy_concordant_species_count", "0"
                ),
                "nonproxy_concordant_species": row.get("nonproxy_concordant_species", ""),
                "nonproxy_discordant_species_count": row.get(
                    "nonproxy_discordant_species_count", "0"
                ),
                "nonproxy_discordant_species": row.get("nonproxy_discordant_species", ""),
                "mapped_symbols": row.get("mapped_symbols", ""),
                "experimental_or_database_supported": row.get(
                    "experimental_or_database_supported", "false"
                ),
                "anchored_effector_ids": row.get("anchored_effector_ids", ""),
                "structurally_supported_effector_count": row.get(
                    "structurally_supported_effector_count", "0"
                ),
                "effector_target_network_priority_class": row.get("effector_target_network_priority_class", ""),
                "apoplastic_cell_wall_priority_class": row.get("apoplastic_cell_wall_priority_class", ""),
                "apoplastic_cell_wall_host_family": node.get("host_family", ""),
                "apoplastic_cell_wall_nonproxy_direction_class": node.get(
                    "nonproxy_direction_class", "not_apoplastic_cell_wall_candidate"
                ),
                "apoplastic_cell_wall_sequence_support_class": node.get(
                    "sequence_support_class", "not_apoplastic_cell_wall_candidate"
                ),
                "cross_species_case_study": "yes" if is_case_study else "no",
                "claim_ceiling": row.get("claim_ceiling", ""),
            }
        )
    return output


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        raise ValueError(f"Cannot write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_apoplastic_cell_wall_report(
    path: Path,
    node_rows: list[dict[str, str]],
    summary: dict[str, str],
) -> None:
    heterogeneous_per = [
        row["orthogroup_id"]
        for row in node_rows
        if row["host_family"] == "PER"
        and row["nonproxy_direction_class"] == "direction_heterogeneous"
    ]
    lines = [
        "# Route 2 independent-node audit",
        "",
        f"- Audited independent orthogroups: {summary['audited_independent_orthogroup_count']}",
        f"- Qualifying independent orthogroups: {summary['qualifying_independent_orthogroup_count']}",
        f"- Qualifying host families: {summary['qualifying_host_family_count']} "
        f"({summary['qualifying_host_families']})",
        f"- Qualifying nonproxy species: {summary['qualifying_nonproxy_species']}",
        f"- Strict module-level multinode support: `{summary['strict_module_multinode_support']}`",
        f"- Tier 1 candidate substrates: {summary['tier1_candidate_substrate_count']} "
        f"({summary['tier1_candidate_substrate_orthogroups']})",
        f"- Peroxidase branch: `{summary['peroxidase_branch_status']}`",
        f"- Direction-heterogeneous PER orthogroups: {';'.join(heterogeneous_per) or 'none'}",
        f"- Direct cleavage evidence: `{summary['direct_cleavage_evidence']}`",
        "",
        "## Interpretation",
        "",
        "The current evidence supports Route 2 at the defense-module level through multiple "
        "independent cell-wall and redox orthogroups. Tier 1 sequence-resolved proteins are "
        "candidate substrates for follow-up, not validated substrates. The peroxidase branch "
        "is retained as direction-heterogeneous auxiliary evidence when concordant-only PER "
        "support is absent.",
        "",
        "This audit does not establish M8 binding, proteolysis, causal hydraulic collapse, "
        "or resistance.",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def write_candidate_report(
    path: Path,
    candidate_rows: list[dict[str, str]],
    apoplastic_cell_wall_summary: dict[str, str],
) -> None:
    leading = sum(row["overall_priority_class"] == "mechanism_leading" for row in candidate_rows)
    supported = sum(
        row["overall_priority_class"] == "mechanism_supported" for row in candidate_rows
    )
    case_studies = [
        row["orthogroup_id"]
        for row in candidate_rows
        if row["cross_species_case_study"] == "yes"
    ]
    lines = [
        "# Manuscript mechanism candidate summary",
        "",
        f"- Manuscript candidates: {len(candidate_rows)}",
        f"- Mechanism leading: {leading}",
        f"- Mechanism supported: {supported}",
        f"- Cross-species case studies: {len(case_studies)}",
        f"- Case-study orthogroups: {';'.join(case_studies) or 'none'}",
        f"- Route 2 strict multinode support: "
        f"`{apoplastic_cell_wall_summary['strict_module_multinode_support']}`",
        f"- Route 2 peroxidase branch: `{apoplastic_cell_wall_summary['peroxidase_branch_status']}`",
        "",
        "The table is a manuscript-prioritization layer. It does not establish direct "
        "effector targeting, pine physical interactions, causal resistance, coevolution, "
        "or geographic adaptation.",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apoplastic_cell_wall-hypotheses", type=Path, required=True)
    parser.add_argument("--host-priorities", type=Path, required=True)
    parser.add_argument("--sequence-shortlist", type=Path, required=True)
    parser.add_argument("--og0005853-audit", type=Path, required=True)
    parser.add_argument("--og0005853-expression", type=Path, required=True)
    parser.add_argument("--apoplastic_cell_wall-audit-output", type=Path, required=True)
    parser.add_argument("--apoplastic_cell_wall-summary-output", type=Path, required=True)
    parser.add_argument("--apoplastic_cell_wall-report-output", type=Path, required=True)
    parser.add_argument("--candidate-output", type=Path, required=True)
    parser.add_argument("--candidate-report-output", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    node_rows = build_apoplastic_cell_wall_node_audit(
        args.apoplastic_cell_wall_hypotheses, args.host_priorities, args.sequence_shortlist
    )
    apoplastic_cell_wall_summary = summarize_apoplastic_cell_wall_nodes(node_rows)
    candidates = build_manuscript_candidates(
        args.host_priorities,
        node_rows,
        args.og0005853_audit,
        args.og0005853_expression,
    )
    write_tsv(args.apoplastic_cell_wall_audit_output, node_rows)
    write_tsv(args.apoplastic_cell_wall_summary_output, [apoplastic_cell_wall_summary])
    write_apoplastic_cell_wall_report(args.apoplastic_cell_wall_report_output, node_rows, apoplastic_cell_wall_summary)
    write_tsv(args.candidate_output, candidates)
    write_candidate_report(args.candidate_report_output, candidates, apoplastic_cell_wall_summary)
    LOGGER.info(
        "Wrote %d Route 2 nodes and %d manuscript candidates",
        len(node_rows),
        len(candidates),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
