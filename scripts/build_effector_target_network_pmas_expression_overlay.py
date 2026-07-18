#!/usr/bin/env python3
"""Overlay exact P. massoniana Route 1 orthogroup members with 1 dpi DESeq2."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


DE_COLUMNS = ["feature_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
FIELDS = [
    "pmas_feature_id", "orthogroup_id", "effector_target_network_role", "host_module", "host_network_role",
    "resistant_log2_fold_change", "resistant_padj", "resistant_fdr_status",
    "susceptible_log2_fold_change", "susceptible_padj", "susceptible_fdr_status", "claim_ceiling",
]
CLAIM = "Exact P. massoniana orthogroup membership plus 1 dpi expression is infection-response evidence only; not direct effector targeting, PPI, or causal resistance evidence"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_de(path: Path) -> dict[str, dict[str, str]]:
    rows = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for values in csv.reader(handle, delimiter="\t"):
            if len(values) != 7:
                raise ValueError(f"Expected 7 DESeq2 columns in {path}, found {len(values)}")
            rows[values[0]] = dict(zip(DE_COLUMNS, values))
    return rows


def fdr_status(row: dict[str, str] | None) -> str:
    if not row:
        return "not_tested_or_filtered"
    try:
        return "fdr_lt_0.05" if float(row["padj"]) < 0.05 else "fdr_ge_0.05"
    except (TypeError, ValueError):
        return "not_tested_or_filtered"


def build_overlay(seeds_path: Path, orthogroups_path: Path, resistant_path: Path, susceptible_path: Path, output_path: Path) -> list[dict[str, str]]:
    seeds = {row["orthogroup_id"]: row for row in read_tsv(seeds_path)}
    resistant, susceptible = read_de(resistant_path), read_de(susceptible_path)
    rows = []
    for group in read_tsv(orthogroups_path):
        og = group["Orthogroup"]
        if og not in seeds:
            continue
        for member in filter(None, (value.strip() for value in group.get("pmas", "").split(","))):
            feature_id = member.removeprefix("pmas|")
            seed, r, s = seeds[og], resistant.get(feature_id), susceptible.get(feature_id)
            rows.append({
                "pmas_feature_id": feature_id, "orthogroup_id": og,
                "effector_target_network_role": seed.get("effector_target_network_role", ""), "host_module": seed.get("host_module", ""), "host_network_role": seed.get("host_network_role", ""),
                "resistant_log2_fold_change": r.get("log2FoldChange", "") if r else "", "resistant_padj": r.get("padj", "") if r else "", "resistant_fdr_status": fdr_status(r),
                "susceptible_log2_fold_change": s.get("log2FoldChange", "") if s else "", "susceptible_padj": s.get("padj", "") if s else "", "susceptible_fdr_status": fdr_status(s),
                "claim_ceiling": CLAIM,
            })
    rows.sort(key=lambda row: (row["orthogroup_id"], row["pmas_feature_id"]))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=Path, required=True)
    parser.add_argument("--orthogroups", type=Path, required=True)
    parser.add_argument("--resistant", type=Path, required=True)
    parser.add_argument("--susceptible", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = build_overlay(args.seeds, args.orthogroups, args.resistant, args.susceptible, args.output)
    print(f"P. massoniana Route 1 expression overlay: {len(rows)} exact members")


if __name__ == "__main__":
    main()
