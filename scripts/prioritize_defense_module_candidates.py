#!/usr/bin/env python3
"""Prioritize defense-module orthogroup candidates for follow-up."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


FIELDS = [
    "priority_rank",
    "priority_tier",
    "priority_score",
    "module_id",
    "orthogroup_id",
    "copy_bias_direction",
    "species_count",
    "gene_count",
    "susceptible_species_count",
    "tolerant_species_count",
    "outgroup_species_count",
    "susceptible_gene_count",
    "tolerant_gene_count",
    "outgroup_gene_count",
    "max_species_gene_count",
    "min_species_gene_count",
    "copy_number_range",
    "mean_pident",
    "mean_bitscore",
    "best_bitscore",
    "matched_keywords",
    "annotation_sources",
    "evidence_types",
    "best_expansion_orthogroup",
    "pilot_core_retention",
    "core_jaccard",
]


def split_values(value: str) -> list[str]:
    if not value.strip():
        return []
    for separator in [";", ","]:
        if separator in value:
            return [part.strip() for part in value.split(separator) if part.strip()]
    return [value.strip()]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_phenotype_groups(path: Path) -> dict[str, str]:
    groups: dict[str, str] = {}
    for row in read_tsv(path):
        response = row.get("pwd_response", "").lower()
        score = row.get("response_score", "")
        if "outgroup" in response:
            group = "outgroup"
        elif score == "0" or response.startswith("susceptible"):
            group = "susceptible"
        elif score == "2" or "tolerant" in response or "resistant" in response:
            group = "tolerant"
        else:
            group = "other"
        groups[row["species_id"]] = group
    return groups


def read_diamond_scores(path: Path) -> dict[str, dict[str, float]]:
    scores: dict[str, dict[str, float]] = {}
    for row in read_tsv(path):
        scores[row["gene_id"]] = {
            "pident": float(row.get("pident", "0") or 0),
            "bitscore": float(row.get("bitscore", "0") or 0),
        }
    return scores


def read_stable_orthogroups(path: Path) -> dict[str, dict[str, str]]:
    return {row["orthogroup_id"]: row for row in read_tsv(path)}


def prioritize_candidates(
    matrix_path: Path,
    gene_hits_path: Path,
    diamond_hits_path: Path,
    stable_orthogroups_path: Path,
    phenotype_path: Path,
    output_path: Path,
) -> list[dict[str, str]]:
    phenotype_groups = read_phenotype_groups(phenotype_path)
    diamond_scores = read_diamond_scores(diamond_hits_path)
    stable_rows = read_stable_orthogroups(stable_orthogroups_path)

    buckets: dict[tuple[str, str], dict[str, object]] = defaultdict(new_bucket)
    for row in read_tsv(matrix_path):
        key = (row["orthogroup_id"], row["module_id"])
        bucket = buckets[key]
        species_id = row["species_id"]
        gene_count = int(row.get("gene_count", "0") or 0)
        bucket["species_ids"].add(species_id)
        bucket["species_gene_counts"][species_id] += gene_count
        bucket["gene_count"] += gene_count
        bucket["matched_keywords"].update(split_values(row.get("matched_keywords", "")))
        bucket["annotation_sources"].update(split_values(row.get("annotation_sources", "")))
        bucket["evidence_types"].update(split_values(row.get("evidence_type", "")))

    for row in read_tsv(gene_hits_path):
        key = (row["orthogroup_id"], row["module_id"])
        bucket = buckets[key]
        gene_id = row["gene_id"]
        bucket["gene_ids"].add(gene_id)
        score = diamond_scores.get(gene_id)
        if score:
            bucket["pidents"].append(score["pident"])
            bucket["bitscores"].append(score["bitscore"])

    rows: list[dict[str, str]] = []
    for (orthogroup_id, module_id), bucket in buckets.items():
        species_gene_counts = bucket["species_gene_counts"]
        group_counts = count_by_group(species_gene_counts, phenotype_groups)
        pidents = bucket["pidents"]
        bitscores = bucket["bitscores"]
        mean_pident = sum(pidents) / len(pidents) if pidents else 0.0
        mean_bitscore = sum(bitscores) / len(bitscores) if bitscores else 0.0
        best_bitscore = max(bitscores) if bitscores else 0.0
        species_counts = list(species_gene_counts.values())
        min_count = min(species_counts) if species_counts else 0
        max_count = max(species_counts) if species_counts else 0
        copy_range = max_count - min_count
        stable = stable_rows.get(orthogroup_id, {})
        score = priority_score(
            species_count=len(bucket["species_ids"]),
            susceptible_species_count=group_counts["susceptible_species_count"],
            tolerant_species_count=group_counts["tolerant_species_count"],
            copy_number_range=copy_range,
            mean_pident=mean_pident,
            mean_bitscore=mean_bitscore,
            pilot_core_retention=float(stable.get("pilot_core_retention", "0") or 0),
            core_jaccard=float(stable.get("core_jaccard", "0") or 0),
        )
        rows.append(
            {
                "priority_rank": "0",
                "priority_tier": priority_tier(score),
                "priority_score": str(score),
                "module_id": module_id,
                "orthogroup_id": orthogroup_id,
                "copy_bias_direction": copy_bias_direction(group_counts),
                "species_count": str(len(bucket["species_ids"])),
                "gene_count": str(bucket["gene_count"]),
                "susceptible_species_count": str(group_counts["susceptible_species_count"]),
                "tolerant_species_count": str(group_counts["tolerant_species_count"]),
                "outgroup_species_count": str(group_counts["outgroup_species_count"]),
                "susceptible_gene_count": str(group_counts["susceptible_gene_count"]),
                "tolerant_gene_count": str(group_counts["tolerant_gene_count"]),
                "outgroup_gene_count": str(group_counts["outgroup_gene_count"]),
                "max_species_gene_count": str(max_count),
                "min_species_gene_count": str(min_count),
                "copy_number_range": str(copy_range),
                "mean_pident": f"{mean_pident:.2f}",
                "mean_bitscore": f"{mean_bitscore:.2f}",
                "best_bitscore": f"{best_bitscore:.2f}",
                "matched_keywords": join_set(bucket["matched_keywords"]),
                "annotation_sources": join_set(bucket["annotation_sources"]),
                "evidence_types": join_set(bucket["evidence_types"]),
                "best_expansion_orthogroup": stable.get("best_expansion_orthogroup", "NA"),
                "pilot_core_retention": stable.get("pilot_core_retention", "NA"),
                "core_jaccard": stable.get("core_jaccard", "NA"),
            }
        )

    rows.sort(
        key=lambda row: (
            -int(row["priority_score"]),
            row["module_id"],
            row["orthogroup_id"],
        )
    )
    for index, row in enumerate(rows, start=1):
        row["priority_rank"] = str(index)
    write_tsv(output_path, FIELDS, rows)
    return rows


def new_bucket() -> dict[str, object]:
    return {
        "species_ids": set(),
        "species_gene_counts": defaultdict(int),
        "gene_count": 0,
        "gene_ids": set(),
        "pidents": [],
        "bitscores": [],
        "matched_keywords": set(),
        "annotation_sources": set(),
        "evidence_types": set(),
    }


def count_by_group(species_gene_counts: dict[str, int], phenotype_groups: dict[str, str]) -> dict[str, int]:
    counts = {
        "susceptible_species_count": 0,
        "tolerant_species_count": 0,
        "outgroup_species_count": 0,
        "susceptible_gene_count": 0,
        "tolerant_gene_count": 0,
        "outgroup_gene_count": 0,
    }
    for species_id, gene_count in species_gene_counts.items():
        group = phenotype_groups.get(species_id, "other")
        if group == "susceptible":
            counts["susceptible_species_count"] += 1
            counts["susceptible_gene_count"] += gene_count
        elif group == "tolerant":
            counts["tolerant_species_count"] += 1
            counts["tolerant_gene_count"] += gene_count
        elif group == "outgroup":
            counts["outgroup_species_count"] += 1
            counts["outgroup_gene_count"] += gene_count
    return counts


def copy_bias_direction(counts: dict[str, int]) -> str:
    susceptible_species = counts["susceptible_species_count"]
    tolerant_species = counts["tolerant_species_count"]
    outgroup_species = counts["outgroup_species_count"]
    if susceptible_species and not tolerant_species:
        return "susceptible_only"
    if tolerant_species and not susceptible_species:
        return "tolerant_only"
    if outgroup_species and not susceptible_species and not tolerant_species:
        return "outgroup_only"
    susceptible_mean = counts["susceptible_gene_count"] / susceptible_species if susceptible_species else 0
    tolerant_mean = counts["tolerant_gene_count"] / tolerant_species if tolerant_species else 0
    if susceptible_mean >= tolerant_mean * 1.5 and susceptible_mean > tolerant_mean:
        return "susceptible_enriched"
    if tolerant_mean >= susceptible_mean * 1.5 and tolerant_mean > susceptible_mean:
        return "tolerant_enriched"
    return "balanced_or_mixed"


def priority_score(
    species_count: int,
    susceptible_species_count: int,
    tolerant_species_count: int,
    copy_number_range: int,
    mean_pident: float,
    mean_bitscore: float,
    pilot_core_retention: float,
    core_jaccard: float,
) -> int:
    score = 0
    score += 2 if species_count >= 4 else 1 if species_count >= 2 else 0
    score += 2 if susceptible_species_count and tolerant_species_count else 0
    score += 2 if copy_number_range >= 3 else 1 if copy_number_range >= 1 else 0
    score += 2 if mean_pident >= 50 else 1 if mean_pident >= 35 else 0
    score += 2 if mean_bitscore >= 150 else 1 if mean_bitscore >= 80 else 0
    score += 1 if pilot_core_retention >= 0.95 else 0
    score += 1 if core_jaccard >= 0.80 else 0
    return score


def priority_tier(score: int) -> str:
    if score >= 9:
        return "high"
    if score >= 6:
        return "medium"
    return "low"


def join_set(value: object) -> str:
    return ";".join(sorted(value))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--gene-hits", type=Path, required=True)
    parser.add_argument("--diamond-hits", type=Path, required=True)
    parser.add_argument("--stable-orthogroups", type=Path, required=True)
    parser.add_argument("--phenotype", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = prioritize_candidates(
        matrix_path=args.matrix,
        gene_hits_path=args.gene_hits,
        diamond_hits_path=args.diamond_hits,
        stable_orthogroups_path=args.stable_orthogroups,
        phenotype_path=args.phenotype,
        output_path=args.output,
    )
    print(f"Defense module candidate priority rows written: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
