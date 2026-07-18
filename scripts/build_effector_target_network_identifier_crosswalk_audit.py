#!/usr/bin/env python3
"""Audit identifier namespace compatibility for route 1 expression supplements."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


SPECIES = ["pden", "pmas", "ptab"]

FIELDS = [
    "species_id",
    "candidate_gene_count",
    "candidate_id_namespace",
    "supplement_source_count",
    "supplement_sheet_count",
    "expression_id_namespace",
    "exact_gene_match_count",
    "exact_gene_match_orthogroups",
    "available_bridge_type",
    "crosswalk_status",
    "recommended_next_step",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_effector_target_network_identifier_crosswalk_audit(
    effector_target_network_path: Path,
    hits_path: Path,
    inventory_path: Path,
    output_path: Path,
    markdown_path: Path,
) -> list[dict[str, str]]:
    effector_target_network_rows = read_tsv(effector_target_network_path)
    hit_rows = read_tsv(hits_path)
    inventory_rows = read_tsv(inventory_path)

    rows = [
        build_species_audit_row(species, effector_target_network_rows, hit_rows, inventory_rows)
        for species in SPECIES
    ]
    write_tsv(output_path, rows)
    write_markdown(markdown_path, rows)
    return rows


def build_species_audit_row(
    species_id: str,
    effector_target_network_rows: list[dict[str, str]],
    hit_rows: list[dict[str, str]],
    inventory_rows: list[dict[str, str]],
) -> dict[str, str]:
    gene_ids = candidate_gene_ids(effector_target_network_rows, species_id)
    candidate_namespace = namespace_summary(gene_ids)
    species_inventory = rows_for_species(inventory_rows, species_id)
    expression_namespace = expression_namespace_summary(species_inventory)
    exact_hits, exact_ogs = exact_match_summary(rows_for_species(hit_rows, species_id))
    bridge_type, status, next_step = classify_crosswalk(
        species_id=species_id,
        candidate_namespace=candidate_namespace,
        expression_namespace=expression_namespace,
        exact_hits=exact_hits,
        supplement_source_count=len({row.get("source_id", "") for row in species_inventory}),
    )
    return {
        "species_id": species_id,
        "candidate_gene_count": str(len(gene_ids)),
        "candidate_id_namespace": candidate_namespace,
        "supplement_source_count": str(len({row.get("source_id", "") for row in species_inventory})),
        "supplement_sheet_count": str(len(species_inventory)),
        "expression_id_namespace": expression_namespace,
        "exact_gene_match_count": str(exact_hits),
        "exact_gene_match_orthogroups": ";".join(exact_ogs),
        "available_bridge_type": bridge_type,
        "crosswalk_status": status,
        "recommended_next_step": next_step,
    }


def candidate_gene_ids(rows: list[dict[str, str]], species_id: str) -> list[str]:
    field = f"{species_id}_gene_ids"
    values: set[str] = set()
    for row in rows:
        for value in split_multi(row.get(field, "")):
            values.add(strip_species_prefix(value))
    return sorted(values)


def rows_for_species(rows: list[dict[str, str]], species_id: str) -> list[dict[str, str]]:
    prefix = f"{species_id}_"
    return [row for row in rows if row.get("source_id", "").startswith(prefix)]


def exact_match_summary(rows: list[dict[str, str]]) -> tuple[int, list[str]]:
    exact_hits = 0
    orthogroups: set[str] = set()
    for row in rows:
        count = parse_int(row.get("exact_gene_match_count", "0"))
        exact_hits += count
        if count > 0 and row.get("orthogroup_id"):
            orthogroups.add(row["orthogroup_id"])
    return exact_hits, sorted(orthogroups)


def expression_namespace_summary(rows: list[dict[str, str]]) -> str:
    namespaces: set[str] = set()
    for row in rows:
        text = " ".join(
            [
                row.get("headers", ""),
                row.get("example_rows", ""),
                row.get("sheet_name", ""),
            ]
        )
        lowered = text.lower()
        if "trinity" in lowered:
            namespaces.add("TRINITY")
        if "arabidopsis" in lowered or "|at" in lowered:
            namespaces.add("Arabidopsis")
        if "plaza" in lowered:
            namespaces.add("PLAZA_Gymno")
        if "pita_" in lowered:
            namespaces.add("PITA")
    if not namespaces:
        return "none_detected"
    return ";".join(sorted(namespaces))


def namespace_summary(gene_ids: list[str]) -> str:
    namespaces = {classify_candidate_namespace(gene_id) for gene_id in gene_ids}
    namespaces.discard("")
    if not namespaces:
        return "none_detected"
    return ";".join(sorted(namespaces))


def classify_candidate_namespace(gene_id: str) -> str:
    if gene_id.startswith("Pd"):
        return "Pd"
    if gene_id.startswith("Pt"):
        return "Pt"
    if gene_id.startswith("gmmutg"):
        return "gmmutg"
    if gene_id.startswith(("STRG", "MSTRG")):
        return "STRG"
    if gene_id.startswith("PITA_"):
        return "PITA"
    if gene_id.startswith("TRINITY_"):
        return "TRINITY"
    return "other"


def classify_crosswalk(
    species_id: str,
    candidate_namespace: str,
    expression_namespace: str,
    exact_hits: int,
    supplement_source_count: int,
) -> tuple[str, str, str]:
    if exact_hits > 0:
        return (
            "exact_candidate_gene_matches",
            "exact_candidate_gene_hits_found",
            "Review exact matches against differential expression context before upgrading candidate-specific support.",
        )
    if supplement_source_count == 0:
        return (
            "none",
            "no_effector_target_network_supplement_source",
            "Do not infer candidate-specific expression support for this species until a route 1 compatible source is added.",
        )
    if species_id == "pden" and "TRINITY" in expression_namespace and "Pd" in candidate_namespace:
        return (
            "transcript_annotation_without_current_proteome_ids",
            "no_exact_crosswalk_trinity_to_pd",
            "Recover transcript FASTA or annotation mapping and run sequence-level alignment from Trinity transcripts to current Pd proteins.",
        )
    if species_id == "pmas" and "PITA" in expression_namespace and (
        "gmmutg" in candidate_namespace or "STRG" in candidate_namespace
    ):
        return (
            "expression_tables_without_current_proteome_ids",
            "no_exact_crosswalk_pita_to_gmmutg_or_strg",
            "Locate author ID mapping between PITA and gmmutg/STRG identifiers, or align PITA transcripts/proteins to the current pmas proteome.",
        )
    return (
        "unresolved_identifier_bridge",
        "no_exact_crosswalk_detected",
        "Treat the supplement evidence as family-level only until a reproducible gene ID bridge is added.",
    )


def split_multi(value: str) -> list[str]:
    return [item.strip() for item in value.replace(",", ";").split(";") if item.strip()]


def strip_species_prefix(value: str) -> str:
    return value.split("|", 1)[1] if "|" in value else value


def parse_int(value: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Route 1 Identifier Crosswalk Audit",
        "",
        "This report audits whether route 1 candidate protein identifiers can be connected directly to downloaded public expression supplements.",
        "",
        "| species | candidate namespace | expression namespace | exact hits | status |",
        "|---|---|---|---:|---|",
    ]
    for row in rows:
        lines.append(
            "| {species_id} | {candidate_id_namespace} | {expression_id_namespace} | {exact_gene_match_count} | {crosswalk_status} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "Interpretation: rows without exact candidate-gene matches remain family-level expression support, not candidate-specific DEG evidence.",
            "",
            "Recommended next steps:",
        ]
    )
    for row in rows:
        lines.append(f"- {row['species_id']}: {row['recommended_next_step']}")
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--effector_target_network", type=Path, required=True)
    parser.add_argument("--hits", type=Path, required=True)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = build_effector_target_network_identifier_crosswalk_audit(
        effector_target_network_path=args.effector_target_network,
        hits_path=args.hits,
        inventory_path=args.inventory,
        output_path=args.output,
        markdown_path=args.markdown,
    )
    unresolved = sum(row["exact_gene_match_count"] == "0" for row in rows)
    print(f"Route 1 identifier crosswalk rows: {len(rows)}; unresolved rows: {unresolved}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
