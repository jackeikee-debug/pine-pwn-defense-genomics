#!/usr/bin/env python3
"""Audit whether core host candidates depend on flagged proteome assemblies."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


PANEL_SPECIES = ("pden", "pmas", "ptab", "plam", "plon", "ptae")
ALL_SPECIES = ("pabi",) + PANEL_SPECIES
SCENARIOS = (
    ("baseline", "All seven proteomes", ()),
    ("exclude_pabi", "Exclude low-completeness P. abies", ("pabi",)),
    ("exclude_ptae", "Exclude low-completeness P. taeda", ("ptae",)),
    ("exclude_plon", "Exclude high-duplication P. longaeva", ("plon",)),
    (
        "quality_restricted",
        "Exclude P. abies, P. taeda, and P. longaeva",
        ("pabi", "ptae", "plon"),
    ),
)
FIELDS = (
    "scenario_id",
    "scenario_label",
    "excluded_species",
    "orthogroup_id",
    "remaining_member_count",
    "remaining_species_count",
    "remaining_panel_species_count",
    "panel_species_available",
    "panel_prevalence_threshold",
    "panel_prevalence_pass",
    "candidate_retained",
    "multi_species_supported",
    "single_excluded_proteome_dependency",
    "selection_path",
    "interpretation",
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def count(row: dict[str, str], species_id: str) -> int:
    value = row.get(f"{species_id}_count", "0").strip()
    return int(float(value or 0))


def proportional_prevalence_threshold(panel_species_available: int) -> int:
    # The original conservative screen required five of six pine proteomes.
    return max(1, round(panel_species_available * 5 / 6))


def build_sensitivity_rows(
    stable_rows: list[dict[str, str]],
    candidate_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    stable_by_id = {row["orthogroup_id"]: row for row in stable_rows}
    candidate_ids = [row["orthogroup_id"] for row in candidate_rows]
    missing = sorted(set(candidate_ids) - set(stable_by_id))
    if missing:
        raise ValueError(f"Candidates missing from stable orthogroups: {', '.join(missing)}")

    output: list[dict[str, str]] = []
    for scenario_id, scenario_label, excluded in SCENARIOS:
        excluded_set = set(excluded)
        retained_species = [sid for sid in ALL_SPECIES if sid not in excluded_set]
        retained_panel = [sid for sid in PANEL_SPECIES if sid not in excluded_set]
        threshold = proportional_prevalence_threshold(len(retained_panel))
        for candidate_id in candidate_ids:
            stable = stable_by_id[candidate_id]
            remaining_counts = {sid: count(stable, sid) for sid in retained_species}
            panel_counts = {sid: count(stable, sid) for sid in retained_panel}
            remaining_species_count = sum(value > 0 for value in remaining_counts.values())
            remaining_panel_species_count = sum(value > 0 for value in panel_counts.values())
            remaining_member_count = sum(remaining_counts.values())
            candidate_retained = remaining_member_count > 0
            prevalence_pass = remaining_panel_species_count >= threshold
            single_dependency = (
                len(excluded) == 1
                and count(stable, excluded[0]) > 0
                and remaining_member_count == 0
            )
            selection_path = (
                "conservative_prevalence_screen"
                if prevalence_pass
                else "expression_mechanism_rescue"
            )
            output.append(
                {
                    "scenario_id": scenario_id,
                    "scenario_label": scenario_label,
                    "excluded_species": ";".join(excluded) or "none",
                    "orthogroup_id": candidate_id,
                    "remaining_member_count": str(remaining_member_count),
                    "remaining_species_count": str(remaining_species_count),
                    "remaining_panel_species_count": str(remaining_panel_species_count),
                    "panel_species_available": str(len(retained_panel)),
                    "panel_prevalence_threshold": str(threshold),
                    "panel_prevalence_pass": "yes" if prevalence_pass else "no",
                    "candidate_retained": "yes" if candidate_retained else "no",
                    "multi_species_supported": "yes" if remaining_species_count >= 2 else "no",
                    "single_excluded_proteome_dependency": "yes" if single_dependency else "no",
                    "selection_path": selection_path,
                    "interpretation": (
                        "candidate retains orthogroup members after exclusion; exclusion does not remove the family"
                        if candidate_retained
                        else "candidate family is lost after exclusion"
                    ),
                }
            )
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stable-orthogroups", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows = build_sensitivity_rows(read_tsv(args.stable_orthogroups), read_tsv(args.candidates))
    write_tsv(args.output, rows)
    single = [row for row in rows if row["scenario_id"].startswith("exclude_")]
    retained = sum(row["candidate_retained"] == "yes" for row in single)
    print(f"Wrote {len(rows)} sensitivity records; {retained}/{len(single)} single-exclusion records retained", file=sys.stderr)


if __name__ == "__main__":
    main()
