#!/usr/bin/env python3
"""Build route 1 expression support table for effector-target network candidates."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


FIELDS = [
    "orthogroup_id",
    "effector_target_network_role",
    "deprioritized_in_effector_target_network",
    "matched_keyword",
    "module_ids",
    "manuscript_section",
    "regional_direction",
    "evidence_tier",
    "pden_gene_ids",
    "pmas_gene_ids",
    "ptab_gene_ids",
    "effector_target_network_host_gene_count",
    "candidate_gene_annotations",
    "literature_hit_sources",
    "literature_hit_count",
    "literature_example_contexts",
    "supplement_downloaded_count",
    "expression_support_level",
    "claim_ceiling",
    "safe_wording",
    "avoid_wording",
]

effector_target_network_AXIS = "effector_target_network"
FOCAL_EXPRESSION_SPECIES = ["pden", "pmas", "ptab"]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_effector_target_network_expression_support(
    mechanism_axis_path: Path,
    candidates_path: Path,
    literature_path: Path,
    stable_genes_path: Path,
    annotations_path: Path,
    supplement_manifest_path: Path,
    output_path: Path,
    markdown_path: Path,
) -> list[dict[str, str]]:
    axis_row = find_effector_target_network_axis(read_tsv(mechanism_axis_path))
    effector_target_network_ids = effector_target_network_candidate_ids(axis_row)
    candidates = {row["orthogroup_id"]: row for row in read_tsv(candidates_path)}
    genes_by_og_species = index_genes(read_tsv(stable_genes_path))
    annotations = {row["gene_id"]: row.get("annotation_text", "") for row in read_tsv(annotations_path)}
    literature_hits = index_literature_hits(read_tsv(literature_path))
    supplement_count = count_downloaded_supplements(read_tsv(supplement_manifest_path))

    rows = []
    for orthogroup_id in effector_target_network_ids:
        candidate = candidates.get(orthogroup_id, {})
        row = build_row(
            orthogroup_id=orthogroup_id,
            axis_row=axis_row,
            candidate=candidate,
            genes_by_species=genes_by_og_species.get(orthogroup_id, {}),
            annotations=annotations,
            literature_hits=literature_hits.get(orthogroup_id, []),
            supplement_count=supplement_count,
        )
        rows.append(row)
    write_tsv(output_path, rows)
    write_markdown(markdown_path, rows)
    return rows


def find_effector_target_network_axis(rows: list[dict[str, str]]) -> dict[str, str]:
    for row in rows:
        if row.get("mechanism_axis") == effector_target_network_AXIS:
            return row
    raise ValueError(f"Missing mechanism axis: {effector_target_network_AXIS}")


def effector_target_network_candidate_ids(axis_row: dict[str, str]) -> list[str]:
    ids = []
    for field in [
        "primary_candidate_orthogroups",
        "supporting_candidate_orthogroups",
        "deprioritized_candidate_orthogroups",
    ]:
        ids.extend(split_multi(axis_row.get(field, "")))
    return sorted(set(ids))


def build_row(
    orthogroup_id: str,
    axis_row: dict[str, str],
    candidate: dict[str, str],
    genes_by_species: dict[str, list[str]],
    annotations: dict[str, str],
    literature_hits: list[dict[str, str]],
    supplement_count: int,
) -> dict[str, str]:
    pden_genes = genes_by_species.get("pden", [])
    pmas_genes = genes_by_species.get("pmas", [])
    ptab_genes = genes_by_species.get("ptab", [])
    effector_target_network_genes = pden_genes + pmas_genes + ptab_genes
    return {
        "orthogroup_id": orthogroup_id,
        "effector_target_network_role": effector_target_network_role(axis_row, orthogroup_id),
        "deprioritized_in_effector_target_network": yes_no(orthogroup_id in split_multi(axis_row.get("deprioritized_candidate_orthogroups", ""))),
        "matched_keyword": candidate.get("matched_keywords", ""),
        "module_ids": candidate.get("module_ids", ""),
        "manuscript_section": candidate.get("manuscript_section", ""),
        "regional_direction": candidate.get("regional_direction", ""),
        "evidence_tier": candidate.get("evidence_tier", ""),
        "pden_gene_ids": join(pden_genes),
        "pmas_gene_ids": join(pmas_genes),
        "ptab_gene_ids": join(ptab_genes),
        "effector_target_network_host_gene_count": str(len(effector_target_network_genes)),
        "candidate_gene_annotations": summarize_annotations(effector_target_network_genes, annotations),
        "literature_hit_sources": join_unique(row.get("source_id", "") for row in literature_hits),
        "literature_hit_count": str(len(literature_hits)),
        "literature_example_contexts": join_unique(row.get("example_context", "") for row in literature_hits),
        "supplement_downloaded_count": str(supplement_count),
        "expression_support_level": support_level(literature_hits),
        "claim_ceiling": (
            "family-level expression literature support; exact candidate gene/transcript differential "
            "expression remains unresolved"
        ),
        "safe_wording": candidate.get("safe_wording", ""),
        "avoid_wording": candidate.get("avoid_wording", ""),
    }


def effector_target_network_role(axis_row: dict[str, str], orthogroup_id: str) -> str:
    if orthogroup_id in split_multi(axis_row.get("primary_candidate_orthogroups", "")):
        return "primary"
    if orthogroup_id in split_multi(axis_row.get("supporting_candidate_orthogroups", "")):
        return "supporting"
    return "deprioritized"


def index_genes(rows: list[dict[str, str]]) -> dict[str, dict[str, list[str]]]:
    index: dict[str, dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))
    for row in rows:
        index[row["orthogroup_id"]][row["species_id"]].append(row["gene_id"])
    return {
        og: {species: sorted(genes) for species, genes in species_map.items()}
        for og, species_map in index.items()
    }


def index_literature_hits(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    hits: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row.get("mechanism_axis") == effector_target_network_AXIS and row.get("support_status") == "literature_keyword_hit":
            hits[row["orthogroup_id"]].append(row)
    return hits


def count_downloaded_supplements(rows: list[dict[str, str]]) -> int:
    return sum(row.get("status") == "ok" for row in rows)


def summarize_annotations(gene_ids: list[str], annotations: dict[str, str]) -> str:
    summaries = [f"{gene_id}:{annotations[gene_id]}" for gene_id in gene_ids if gene_id in annotations]
    return join_unique(summaries)


def support_level(literature_hits: list[dict[str, str]]) -> str:
    if literature_hits:
        return "family_level_literature_support"
    return "no_literature_keyword_hit"


def split_multi(value: str) -> list[str]:
    return [item for item in value.replace(",", ";").split(";") if item]


def join(values: list[str]) -> str:
    return ";".join(values)


def join_unique(values) -> str:
    return ";".join(sorted({value for value in values if value}))


def yes_no(value: bool) -> str:
    return "yes" if value else "no"


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    supported = [row for row in rows if row["expression_support_level"] == "family_level_literature_support"]
    supplement_count = rows[0]["supplement_downloaded_count"] if rows else "0"
    lines = [
        "# Route 1 Expression Support",
        "",
        "This report summarizes expression support for the effector-target network vulnerability route.",
        "Current support is family-level literature evidence unless a future parser links exact genes or transcripts.",
        "",
        f"- Route 1 candidate orthogroups: {len(rows)}",
        f"- Family-level literature-supported orthogroups: {len(supported)}",
        f"- Downloaded supplementary files available for parsing: {supplement_count}",
        "",
    ]
    for row in rows:
        lines.append(
            f"- {row['orthogroup_id']} ({row['effector_target_network_role']}, {row['matched_keyword']}): "
            f"{row['expression_support_level']}; sources={row['literature_hit_sources'] or 'none'}"
        )
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mechanism-axis", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--literature", type=Path, required=True)
    parser.add_argument("--stable-genes", type=Path, required=True)
    parser.add_argument("--annotations", type=Path, required=True)
    parser.add_argument("--supplement-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = build_effector_target_network_expression_support(
        mechanism_axis_path=args.mechanism_axis,
        candidates_path=args.candidates,
        literature_path=args.literature,
        stable_genes_path=args.stable_genes,
        annotations_path=args.annotations,
        supplement_manifest_path=args.supplement_manifest,
        output_path=args.output,
        markdown_path=args.markdown,
    )
    supported = sum(row["expression_support_level"] == "family_level_literature_support" for row in rows)
    print(f"Route 1 expression support rows written: {len(rows)}; literature-supported: {supported}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
