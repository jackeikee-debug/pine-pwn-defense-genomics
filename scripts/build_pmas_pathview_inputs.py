#!/usr/bin/env python3
"""Prepare conservative TAIR-locus fold changes for P. massoniana Pathview overlays."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path


HEADERLESS_FIELDS = ["feature_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
OUTPUT_FIELDS = [
    "tair_locus",
    "contrast_id",
    "pmas_gene_id",
    "log2_fold_change",
    "padj",
    "collapse_count",
]


def safe_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def read_deseq2(path: Path, contrast_id: str) -> list[dict[str, object]]:
    with path.open(encoding="utf-8", newline="") as handle:
        first_line = handle.readline()
        handle.seek(0)
        has_header = first_line.split("\t", 1)[0] == "feature_id"
        reader = csv.DictReader(
            handle,
            delimiter="\t",
            fieldnames=None if has_header else HEADERLESS_FIELDS,
        )
        rows: list[dict[str, object]] = []
        for row in reader:
            if contrast_id == "genotype_by_inoculum":
                lfc = row.get("interaction_log2_fold_change", "")
                padj = row.get("interaction_padj", "")
            else:
                lfc = row.get("log2FoldChange", "")
                padj = row.get("padj", "")
            rows.append(
                {
                    "feature_id": row.get("feature_id", ""),
                    "contrast_id": contrast_id,
                    "log2_fold_change": safe_float(str(lfc)),
                    "padj": safe_float(str(padj)),
                }
            )
    return rows


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def format_float(value: float) -> str:
    return "" if not math.isfinite(value) else f"{value:.12g}"


def aggregate_to_tair(
    de_rows: list[dict[str, object]],
    crosswalk: list[dict[str, str]],
) -> list[dict[str, str]]:
    tair_by_gene = {
        row["pmas_gene_id"]: row["tair_locus"].upper()
        for row in crosswalk
        if row.get("pmas_gene_id") and row.get("tair_locus")
    }
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in de_rows:
        tair = tair_by_gene.get(str(row["feature_id"]), "")
        lfc = float(row["log2_fold_change"])
        if tair and math.isfinite(lfc):
            grouped[(tair, str(row["contrast_id"]))].append(row)

    output: list[dict[str, str]] = []
    for (tair, contrast), rows in sorted(grouped.items()):
        representative = max(
            rows,
            key=lambda row: (abs(float(row["log2_fold_change"])), str(row["feature_id"])),
        )
        output.append(
            {
                "tair_locus": tair,
                "contrast_id": contrast,
                "pmas_gene_id": str(representative["feature_id"]),
                "log2_fold_change": format_float(float(representative["log2_fold_change"])),
                "padj": format_float(float(representative["padj"])),
                "collapse_count": str(len(rows)),
            }
        )
    return output


def assess_coverage(mapped_loci: set[str], minimum: int) -> str:
    return "render" if len(mapped_loci) >= minimum else "insufficient_coverage"


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resistant", type=Path, required=True)
    parser.add_argument("--susceptible", type=Path, required=True)
    parser.add_argument("--interaction", type=Path, required=True)
    parser.add_argument("--crosswalk", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    de_rows = []
    de_rows.extend(read_deseq2(args.resistant, "resistant_pwn_vs_water"))
    de_rows.extend(read_deseq2(args.susceptible, "susceptible_pwn_vs_water"))
    de_rows.extend(read_deseq2(args.interaction, "genotype_by_inoculum"))
    output = aggregate_to_tair(de_rows, read_tsv(args.crosswalk))
    write_tsv(args.output, output)
    contrast_counts: dict[str, int] = defaultdict(int)
    for row in output:
        contrast_counts[row["contrast_id"]] += 1
    print(
        "Wrote ortholog-projected Pathview values: "
        + "; ".join(f"{key}={value}" for key, value in sorted(contrast_counts.items())),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
