#!/usr/bin/env python3
"""Build traceable edge and sector tables for the manuscript mechanism chord diagram."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path


EDGE_FIELDS = ["source_id", "target_id", "edge_type", "weight", "evidence_boundary"]
SECTOR_FIELDS = [
    "sector_id",
    "sector_type",
    "label",
    "group_order",
    "display_order",
    "count",
    "priority_class",
    "expression_support",
    "network_support",
    "structure_support",
    "sequence_support",
    "direct_interaction_evidence",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def split_values(value: str) -> list[str]:
    return [item.strip() for item in value.split(";") if item.strip()]


def as_int(value: str) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def is_true(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def expression_support(row: dict[str, str]) -> str:
    count = as_int(row.get("nonproxy_concordant_species_count", "0"))
    if count >= 2:
        return "cross_species_concordant"
    if count == 1:
        return "single_species_concordant"
    return "none"


def build_chord_inputs(
    candidates: list[dict[str, str]],
    links: list[dict[str, str]],
    effectors: list[dict[str, str]],
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    candidate_ids = {row.get("orthogroup_id", "") for row in candidates}
    if "" in candidate_ids or len(candidate_ids) != len(candidates):
        raise ValueError("Candidate orthogroup IDs must be non-empty and unique")

    effector_counts = Counter(row.get("effector_class", "") for row in effectors)
    effector_counts.pop("", None)
    declared_classes = set(effector_counts)
    linked_modules = {row.get("host_module", "") for row in links if row.get("host_module", "")}
    candidate_modules = {
        module
        for candidate in candidates
        for module in split_values(candidate.get("module_ids", ""))
    }
    all_modules = linked_modules | candidate_modules

    edges: list[dict[str, str]] = []
    edge_keys: set[tuple[str, str, str]] = set()
    for link in links:
        effector_class = link.get("effector_class", "")
        host_module = link.get("host_module", "")
        if effector_class not in declared_classes:
            raise ValueError(f"undeclared effector class in functional prior: {effector_class}")
        if not host_module:
            raise ValueError("Functional-prior host module must be non-empty")
        key = (effector_class, host_module, "functional_prior")
        if key in edge_keys:
            raise ValueError(f"Duplicate functional-prior edge: {effector_class} -> {host_module}")
        edge_keys.add(key)
        edges.append(
            {
                "source_id": effector_class,
                "target_id": host_module,
                "edge_type": "functional_prior",
                "weight": str(max(1, as_int(link.get("host_candidate_count", "1")))),
                "evidence_boundary": "functional_hypothesis_not_interaction",
            }
        )

    for candidate in candidates:
        orthogroup_id = candidate["orthogroup_id"]
        modules = split_values(candidate.get("module_ids", ""))
        if not modules:
            raise ValueError(f"Candidate has no declared module: {orthogroup_id}")
        for module in modules:
            if module not in all_modules:
                raise ValueError(f"undeclared host module for candidate {orthogroup_id}: {module}")
            key = (module, orthogroup_id, "candidate_membership")
            if key in edge_keys:
                continue
            edge_keys.add(key)
            edges.append(
                {
                    "source_id": module,
                    "target_id": orthogroup_id,
                    "edge_type": "candidate_membership",
                    "weight": "1",
                    "evidence_boundary": "module_annotation_not_interaction",
                }
            )

    sectors: list[dict[str, str]] = []
    display_order = 1
    for effector_class in sorted(declared_classes):
        sectors.append(
            base_sector(
                effector_class,
                "effector_class",
                effector_class.replace("_", " "),
                1,
                display_order,
                effector_counts[effector_class],
            )
        )
        display_order += 1
    for module in sorted(all_modules):
        sectors.append(
            base_sector(module, "host_module", module.replace("_", " "), 2, display_order, 0)
        )
        display_order += 1
    for candidate in candidates:
        symbol = split_values(candidate.get("mapped_symbols", ""))
        label = candidate["orthogroup_id"]
        if symbol:
            label = f"{label} | {symbol[0]}"
        sector = base_sector(label=candidate["orthogroup_id"], sector_id=candidate["orthogroup_id"], sector_type="candidate", group_order=3, display_order=display_order, count=1)
        sector.update(
            {
                "label": label,
                "priority_class": candidate.get("overall_priority_class", ""),
                "expression_support": expression_support(candidate),
                "network_support": (
                    "experimental_or_database_supported"
                    if is_true(candidate.get("experimental_or_database_supported", ""))
                    else "none"
                ),
                "structure_support": (
                    "structurally_supported_anchor"
                    if as_int(candidate.get("structurally_supported_effector_count", "0")) > 0
                    else "none"
                ),
                "sequence_support": candidate.get("apoplastic_cell_wall_sequence_support_class", "") or "none",
                "direct_interaction_evidence": "none",
            }
        )
        sectors.append(sector)
        display_order += 1

    return edges, sectors


def base_sector(
    sector_id: str,
    sector_type: str,
    label: str,
    group_order: int,
    display_order: int,
    count: int,
) -> dict[str, str]:
    return {
        "sector_id": sector_id,
        "sector_type": sector_type,
        "label": label,
        "group_order": str(group_order),
        "display_order": str(display_order),
        "count": str(count),
        "priority_class": "",
        "expression_support": "",
        "network_support": "",
        "structure_support": "",
        "sequence_support": "",
        "direct_interaction_evidence": "none",
    }


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--links", type=Path, required=True)
    parser.add_argument("--effectors", type=Path, required=True)
    parser.add_argument("--edges-output", type=Path, required=True)
    parser.add_argument("--sectors-output", type=Path, required=True)
    args = parser.parse_args()

    edges, sectors = build_chord_inputs(
        read_tsv(args.candidates),
        read_tsv(args.links),
        read_tsv(args.effectors),
    )
    write_tsv(args.edges_output, EDGE_FIELDS, edges)
    write_tsv(args.sectors_output, SECTOR_FIELDS, sectors)
    print(
        f"Wrote {len(edges)} traceable chord edges and {len(sectors)} sectors",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
