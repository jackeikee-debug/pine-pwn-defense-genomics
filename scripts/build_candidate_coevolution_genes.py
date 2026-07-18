#!/usr/bin/env python3
"""Expand effector-host module links into candidate linked host orthogroups."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


FIELDS = [
    "orthogroup_id",
    "module_id",
    "effector_class",
    "candidate_signal",
    "interpretation",
    "link_type",
    "evidence_basis",
    "link_confidence",
    "effector_count",
    "effector_ids",
    "effector_confidences",
    "effector_evidence",
    "host_candidate_count",
    "host_bias_direction",
    "host_priority_rank",
    "host_priority_tier",
    "host_priority_score",
    "host_shortlist_rank",
    "host_shortlist_group_rank",
    "host_gene_count",
    "susceptible_gene_count",
    "tolerant_gene_count",
    "outgroup_gene_count",
    "species_count",
    "matched_keywords",
    "annotation_sources",
    "evidence_types",
    "mean_bitscore",
    "best_bitscore",
    "pilot_core_retention",
    "core_jaccard",
    "best_expansion_orthogroup",
    "notes",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_candidate_coevolution_genes(
    links_path: Path,
    host_shortlist_path: Path,
    effectors_path: Path,
    output_path: Path,
) -> list[dict[str, str]]:
    effectors_by_id = {row["protein_id"]: row for row in read_tsv(effectors_path)}
    host_rows = read_tsv(host_shortlist_path)
    host_by_module: dict[str, list[dict[str, str]]] = {}
    for row in host_rows:
        host_by_module.setdefault(row["module_id"], []).append(row)

    rows: list[dict[str, str]] = []
    for link in read_tsv(links_path):
        module_id = link["host_module"]
        linked_orthogroups = set(split_values(link.get("host_orthogroups", "")))
        effector_ids = split_values(link.get("effector_ids", ""))
        matching_hosts = [
            row
            for row in host_by_module.get(module_id, [])
            if not linked_orthogroups or row["orthogroup_id"] in linked_orthogroups
        ]
        matching_hosts.sort(key=lambda row: int_or_max(row.get("shortlist_rank", "")))

        for host in matching_hosts:
            rows.append(make_candidate_row(link, host, effector_ids, effectors_by_id))

    write_tsv(output_path, FIELDS, rows)
    return rows


def make_candidate_row(
    link: dict[str, str],
    host: dict[str, str],
    effector_ids: list[str],
    effectors_by_id: dict[str, dict[str, str]],
) -> dict[str, str]:
    module_id = link["host_module"]
    effector_class = link["effector_class"]
    confidence_values = {
        effectors_by_id[protein_id]["confidence"]
        for protein_id in effector_ids
        if protein_id in effectors_by_id and effectors_by_id[protein_id].get("confidence")
    }
    evidence_values = {
        effectors_by_id[protein_id]["evidence"]
        for protein_id in effector_ids
        if protein_id in effectors_by_id and effectors_by_id[protein_id].get("evidence")
    }
    if not evidence_values:
        evidence_values = set(split_values(link.get("effector_evidence", "")))

    return {
        "orthogroup_id": host["orthogroup_id"],
        "module_id": module_id,
        "effector_class": effector_class,
        "candidate_signal": "candidate_functional_link",
        "interpretation": build_interpretation(link, host),
        "link_type": link.get("link_type", ""),
        "evidence_basis": link.get("evidence_basis", ""),
        "link_confidence": link.get("link_confidence", "candidate"),
        "effector_count": link.get("effector_count", str(len(effector_ids))),
        "effector_ids": ";".join(effector_ids),
        "effector_confidences": join_values(confidence_values),
        "effector_evidence": join_values(evidence_values),
        "host_candidate_count": link.get("host_candidate_count", ""),
        "host_bias_direction": host.get("copy_bias_direction", ""),
        "host_priority_rank": host.get("priority_rank", ""),
        "host_priority_tier": host.get("priority_tier", ""),
        "host_priority_score": host.get("priority_score", ""),
        "host_shortlist_rank": host.get("shortlist_rank", ""),
        "host_shortlist_group_rank": host.get("shortlist_group_rank", ""),
        "host_gene_count": host.get("gene_count", ""),
        "susceptible_gene_count": host.get("susceptible_gene_count", ""),
        "tolerant_gene_count": host.get("tolerant_gene_count", ""),
        "outgroup_gene_count": host.get("outgroup_gene_count", ""),
        "species_count": host.get("species_count", ""),
        "matched_keywords": host.get("matched_keywords", ""),
        "annotation_sources": host.get("annotation_sources", ""),
        "evidence_types": host.get("evidence_types", ""),
        "mean_bitscore": host.get("mean_bitscore", ""),
        "best_bitscore": host.get("best_bitscore", ""),
        "pilot_core_retention": host.get("pilot_core_retention", ""),
        "core_jaccard": host.get("core_jaccard", ""),
        "best_expansion_orthogroup": host.get("best_expansion_orthogroup", ""),
        "notes": link.get("notes", ""),
    }


def build_interpretation(link: dict[str, str], host: dict[str, str]) -> str:
    return (
        f"Candidate linked module between {link['effector_class']} nematode effectors and "
        f"{host.get('copy_bias_direction', 'directional')} host {link['host_module']} "
        f"orthogroup {host['orthogroup_id']}; not direct interaction or coevolution evidence."
    )


def split_values(value: str) -> list[str]:
    if not value.strip():
        return []
    for separator in [";", ","]:
        if separator in value:
            return [part.strip() for part in value.split(separator) if part.strip()]
    return [value.strip()]


def join_values(values: set[str]) -> str:
    return ";".join(sorted(values))


def int_or_max(value: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return 10**9


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--links", type=Path, required=True)
    parser.add_argument("--host-shortlist", type=Path, required=True)
    parser.add_argument("--effectors", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = build_candidate_coevolution_genes(
        links_path=args.links,
        host_shortlist_path=args.host_shortlist,
        effectors_path=args.effectors,
        output_path=args.output,
    )
    print(f"Candidate linked host orthogroups written: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
