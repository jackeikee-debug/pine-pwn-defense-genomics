#!/usr/bin/env python3
"""Freeze counts and set intersections for the PWN secretome audit figure."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


SOURCE_IDS = (
    "shinya2013_secretome",
    "cardoso2016_bxpe",
    "silva2021_stimulus_secretome",
)
STAGE_FIELDS = ("stage_id", "stage_label", "protein_count", "set_relation", "interpretation")
COMBINATION_FIELDS = (
    "combination_id",
    "shinya",
    "cardoso",
    "silva",
    "published_support",
    "protein_count",
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def yes(row: dict[str, str], field: str) -> bool:
    return row.get(field, "").strip().lower() == "yes"


def build_tables(
    integrated: list[dict[str, str]],
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    if len(integrated) != 2164:
        raise ValueError(f"Expected 2,164 standard secreted-soluble candidates, found {len(integrated)}")

    functional_count = sum(yes(row, "domain_support") for row in integrated)
    published_count = sum(yes(row, "published_secretome_support") for row in integrated)
    if functional_count != 1265 or published_count != 401:
        raise ValueError(
            f"Frozen audit counts changed: functional={functional_count}, published={published_count}"
        )

    stages = [
        {
            "stage_id": "total_proteome",
            "stage_label": "BXYJv5 proteins",
            "protein_count": "15860",
            "set_relation": "starting_set",
            "interpretation": "all annotated proteins",
        },
        {
            "stage_id": "signal_peptide",
            "stage_label": "DeepSig-supported signal peptide",
            "protein_count": "2741",
            "set_relation": "nested_filter",
            "interpretation": "signal-peptide prediction",
        },
        {
            "stage_id": "secreted_soluble",
            "stage_label": "Secreted-soluble candidates",
            "protein_count": "2164",
            "set_relation": "nested_filter",
            "interpretation": "no mature-protein transmembrane helix",
        },
        {
            "stage_id": "functional_domain",
            "stage_label": "Functional-domain support",
            "protein_count": str(functional_count),
            "set_relation": "independent_support_within_2164",
            "interpretation": "qualifying InterPro signatures",
        },
        {
            "stage_id": "published_secretome",
            "stage_label": "Published-secretome support",
            "protein_count": str(published_count),
            "set_relation": "independent_support_within_2164",
            "interpretation": "conservative mapping to three proteomic studies",
        },
    ]

    combinations: Counter[tuple[bool, bool, bool]] = Counter()
    for row in integrated:
        sources = {
            item.strip()
            for item in row.get("published_secretome_sources", "").split(";")
            if item.strip() and item.strip() != "NA"
        }
        combinations[tuple(source in sources for source in SOURCE_IDS)] += 1

    combination_rows = []
    for membership, protein_count in sorted(
        combinations.items(), key=lambda item: (-item[1], item[0])
    ):
        labels = [name for name, present in zip(("Shinya", "Cardoso", "Silva"), membership) if present]
        combination_rows.append(
            {
                "combination_id": " + ".join(labels) if labels else "Current candidates only",
                "shinya": "yes" if membership[0] else "no",
                "cardoso": "yes" if membership[1] else "no",
                "silva": "yes" if membership[2] else "no",
                "published_support": "yes" if any(membership) else "no",
                "protein_count": str(protein_count),
            }
        )
    return stages, combination_rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--integrated", type=Path, required=True)
    parser.add_argument("--stages-output", type=Path, required=True)
    parser.add_argument("--combinations-output", type=Path, required=True)
    args = parser.parse_args()

    stages, combinations = build_tables(read_tsv(args.integrated))
    write_tsv(args.stages_output, STAGE_FIELDS, stages)
    write_tsv(args.combinations_output, COMBINATION_FIELDS, combinations)


if __name__ == "__main__":
    main()
