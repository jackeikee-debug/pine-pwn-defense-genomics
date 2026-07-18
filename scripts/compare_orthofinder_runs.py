#!/usr/bin/env python3
"""Compare core orthogroup stability between two OrthoFinder runs."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


STABILITY_FIELDS = [
    "pilot_orthogroup",
    "best_expansion_orthogroup",
    "pilot_core_genes",
    "best_shared_core_genes",
    "expansion_core_genes",
    "pilot_core_retention",
    "expansion_core_purity",
    "core_jaccard",
    "status",
]


SUMMARY_FIELDS = [
    "pilot_orthogroups",
    "stable_orthogroups",
    "split_or_partial_orthogroups",
    "unmapped_orthogroups",
    "stable_fraction",
    "mean_core_retention",
    "mean_core_jaccard",
]


def split_genes(value: str) -> list[str]:
    if not value.strip():
        return []
    return [gene.strip() for gene in value.split(",") if gene.strip()]


def read_orthogroups(path: Path, species: list[str]) -> dict[str, set[str]]:
    groups: dict[str, set[str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            genes: set[str] = set()
            for species_id in species:
                genes.update(split_genes(row.get(species_id, "")))
            groups[row["Orthogroup"]] = genes
    return groups


def compare_runs(
    pilot_path: Path,
    expansion_path: Path,
    core_species: list[str],
    output_path: Path,
    summary_path: Path,
    stable_retention_threshold: float = 0.95,
    stable_jaccard_threshold: float = 0.80,
) -> tuple[list[dict[str, str]], dict[str, int | str]]:
    pilot_groups = read_orthogroups(pilot_path, core_species)
    expansion_groups = read_orthogroups(expansion_path, core_species)

    gene_to_expansion_og: dict[str, str] = {}
    for orthogroup, genes in expansion_groups.items():
        for gene in genes:
            gene_to_expansion_og[gene] = orthogroup

    rows: list[dict[str, str]] = []
    retention_values: list[float] = []
    jaccard_values: list[float] = []
    status_counts: Counter[str] = Counter()

    for pilot_og, pilot_genes in pilot_groups.items():
        overlap_counts = Counter(gene_to_expansion_og[gene] for gene in pilot_genes if gene in gene_to_expansion_og)
        if not pilot_genes or not overlap_counts:
            best_expansion_og = "NA"
            shared = 0
            expansion_genes: set[str] = set()
        else:
            best_expansion_og, shared = overlap_counts.most_common(1)[0]
            expansion_genes = expansion_groups[best_expansion_og]

        union_count = len(pilot_genes | expansion_genes)
        retention = shared / len(pilot_genes) if pilot_genes else 0.0
        purity = shared / len(expansion_genes) if expansion_genes else 0.0
        jaccard = shared / union_count if union_count else 0.0

        if best_expansion_og == "NA":
            status = "unmapped"
        elif retention >= stable_retention_threshold and jaccard >= stable_jaccard_threshold:
            status = "stable"
        else:
            status = "split_or_partial"

        status_counts[status] += 1
        retention_values.append(retention)
        jaccard_values.append(jaccard)
        rows.append(
            {
                "pilot_orthogroup": pilot_og,
                "best_expansion_orthogroup": best_expansion_og,
                "pilot_core_genes": str(len(pilot_genes)),
                "best_shared_core_genes": str(shared),
                "expansion_core_genes": str(len(expansion_genes)),
                "pilot_core_retention": f"{retention:.3f}",
                "expansion_core_purity": f"{purity:.3f}",
                "core_jaccard": f"{jaccard:.3f}",
                "status": status,
            }
        )

    summary_row: dict[str, int | str] = {
        "pilot_orthogroups": len(rows),
        "stable_orthogroups": status_counts["stable"],
        "split_or_partial_orthogroups": status_counts["split_or_partial"],
        "unmapped_orthogroups": status_counts["unmapped"],
        "stable_fraction": f"{status_counts['stable'] / len(rows):.3f}" if rows else "0.000",
        "mean_core_retention": f"{sum(retention_values) / len(retention_values):.3f}" if retention_values else "0.000",
        "mean_core_jaccard": f"{sum(jaccard_values) / len(jaccard_values):.3f}" if jaccard_values else "0.000",
    }

    write_tsv(output_path, STABILITY_FIELDS, rows)
    write_tsv(summary_path, SUMMARY_FIELDS, [summary_row])
    return rows, summary_row


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot-orthogroups", type=Path, required=True)
    parser.add_argument("--expansion-orthogroups", type=Path, required=True)
    parser.add_argument("--core-species", nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    compare_runs(
        pilot_path=args.pilot_orthogroups,
        expansion_path=args.expansion_orthogroups,
        core_species=args.core_species,
        output_path=args.output,
        summary_path=args.summary,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
