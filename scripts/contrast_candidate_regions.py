#!/usr/bin/env python3
"""Contrast candidate host orthogroup copy number between East Asia and North America."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


FIELDS = [
    "orthogroup_id",
    "module_id",
    "effector_class",
    "regional_copy_direction",
    "regional_presence_pattern",
    "contrast_basis",
    "interpretation",
    "east_asia_total_copy",
    "north_america_total_copy",
    "outgroup_total_copy",
    "east_asia_species_count",
    "north_america_species_count",
    "east_asia_species_present",
    "north_america_species_present",
    "east_asia_mean_copy",
    "north_america_mean_copy",
    "log2_east_asia_vs_north_america_mean_copy_ratio",
    "east_asia_species_counts",
    "north_america_species_counts",
    "outgroup_species_counts",
    "host_bias_direction",
    "host_priority_score",
    "host_priority_rank",
    "host_shortlist_rank",
    "matched_keywords",
    "mean_bitscore",
    "best_bitscore",
    "effector_count",
    "effector_ids",
    "effector_confidences",
    "effector_evidence",
    "candidate_signal",
    "link_type",
    "evidence_basis",
    "link_confidence",
    "notes",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def contrast_candidate_regions(
    candidates_path: Path,
    stable_orthogroups_path: Path,
    species_metadata_path: Path,
    output_path: Path,
    ratio_threshold: float = 1.5,
    pseudocount: float = 0.5,
) -> list[dict[str, str]]:
    species_by_region = collect_species_by_region(species_metadata_path)
    stable_rows = read_tsv(stable_orthogroups_path)
    validate_species_count_columns(stable_rows, species_by_region, stable_orthogroups_path)
    stable_by_og = {row["orthogroup_id"]: row for row in stable_rows}
    log2_threshold = math.log2(ratio_threshold)

    rows: list[dict[str, str]] = []
    for candidate in read_tsv(candidates_path):
        stable = stable_by_og.get(candidate["orthogroup_id"], {})
        contrast = build_contrast_row(
            candidate=candidate,
            stable=stable,
            species_by_region=species_by_region,
            log2_threshold=log2_threshold,
            pseudocount=pseudocount,
        )
        rows.append(contrast)

    rows.sort(
        key=lambda row: (
            row["module_id"],
            row["effector_class"],
            row["regional_copy_direction"],
            -abs(to_float(row["log2_east_asia_vs_north_america_mean_copy_ratio"])),
            to_int(row["host_shortlist_rank"]),
            row["orthogroup_id"],
        )
    )
    write_tsv(output_path, FIELDS, rows)
    return rows


def collect_species_by_region(path: Path) -> dict[str, list[str]]:
    grouped: dict[str, list[str]] = {"East_Asia": [], "North_America": [], "Outgroup": []}
    for row in read_tsv(path):
        species_id = row["species_id"]
        if row.get("protein_fasta", "NA") == "NA":
            continue
        region = row.get("region", "")
        if region in grouped:
            grouped[region].append(species_id)
    for region in grouped:
        grouped[region].sort()
    return grouped


def validate_species_count_columns(
    stable_rows: list[dict[str, str]],
    species_by_region: dict[str, list[str]],
    stable_orthogroups_path: Path,
) -> None:
    if not stable_rows:
        return
    available_fields = set(stable_rows[0])
    required_species = [
        species_id
        for region in ("East_Asia", "North_America", "Outgroup")
        for species_id in species_by_region[region]
    ]
    missing = [
        f"{species_id}_count"
        for species_id in required_species
        if f"{species_id}_count" not in available_fields
    ]
    if missing:
        missing_text = ", ".join(missing)
        raise ValueError(
            f"{stable_orthogroups_path} is missing required species count columns: "
            f"{missing_text}. Rebuild the orthogroup tables or remove species without "
            "orthogroup counts from the focal metadata before regional contrast."
        )


def build_contrast_row(
    candidate: dict[str, str],
    stable: dict[str, str],
    species_by_region: dict[str, list[str]],
    log2_threshold: float,
    pseudocount: float,
) -> dict[str, str]:
    east_counts = counts_for_species(stable, species_by_region["East_Asia"])
    north_counts = counts_for_species(stable, species_by_region["North_America"])
    outgroup_counts = counts_for_species(stable, species_by_region["Outgroup"])

    east_total = sum(east_counts.values())
    north_total = sum(north_counts.values())
    outgroup_total = sum(outgroup_counts.values())
    east_mean = mean(list(east_counts.values()))
    north_mean = mean(list(north_counts.values()))
    log2_ratio = math.log2((east_mean + pseudocount) / (north_mean + pseudocount))
    direction = classify_direction(log2_ratio, log2_threshold)
    pattern = classify_presence_pattern(east_total, north_total)

    return {
        "orthogroup_id": candidate["orthogroup_id"],
        "module_id": candidate["module_id"],
        "effector_class": candidate["effector_class"],
        "regional_copy_direction": direction,
        "regional_presence_pattern": pattern,
        "contrast_basis": "mean_copy_per_available_species;not_statistical_significance",
        "interpretation": build_interpretation(candidate, direction, pattern, log2_ratio),
        "east_asia_total_copy": str(east_total),
        "north_america_total_copy": str(north_total),
        "outgroup_total_copy": str(outgroup_total),
        "east_asia_species_count": str(len(east_counts)),
        "north_america_species_count": str(len(north_counts)),
        "east_asia_species_present": str(sum(1 for count in east_counts.values() if count > 0)),
        "north_america_species_present": str(sum(1 for count in north_counts.values() if count > 0)),
        "east_asia_mean_copy": format_float(east_mean),
        "north_america_mean_copy": format_float(north_mean),
        "log2_east_asia_vs_north_america_mean_copy_ratio": format_float(log2_ratio),
        "east_asia_species_counts": format_species_counts(east_counts),
        "north_america_species_counts": format_species_counts(north_counts),
        "outgroup_species_counts": format_species_counts(outgroup_counts),
        "host_bias_direction": candidate.get("host_bias_direction", ""),
        "host_priority_score": candidate.get("host_priority_score", ""),
        "host_priority_rank": candidate.get("host_priority_rank", ""),
        "host_shortlist_rank": candidate.get("host_shortlist_rank", ""),
        "matched_keywords": candidate.get("matched_keywords", ""),
        "mean_bitscore": candidate.get("mean_bitscore", ""),
        "best_bitscore": candidate.get("best_bitscore", ""),
        "effector_count": candidate.get("effector_count", ""),
        "effector_ids": candidate.get("effector_ids", ""),
        "effector_confidences": candidate.get("effector_confidences", ""),
        "effector_evidence": candidate.get("effector_evidence", ""),
        "candidate_signal": candidate.get("candidate_signal", ""),
        "link_type": candidate.get("link_type", ""),
        "evidence_basis": candidate.get("evidence_basis", ""),
        "link_confidence": candidate.get("link_confidence", ""),
        "notes": candidate.get("notes", ""),
    }


def counts_for_species(stable: dict[str, str], species_ids: list[str]) -> dict[str, int]:
    return {species_id: to_int(stable.get(f"{species_id}_count", "0")) for species_id in species_ids}


def classify_direction(log2_ratio: float, log2_threshold: float) -> str:
    if log2_ratio >= log2_threshold:
        return "East_Asia_enriched"
    if log2_ratio <= -log2_threshold:
        return "North_America_enriched"
    return "balanced"


def classify_presence_pattern(east_total: int, north_total: int) -> str:
    if east_total > 0 and north_total > 0:
        return "shared"
    if east_total > 0:
        return "East_Asia_only"
    if north_total > 0:
        return "North_America_only"
    return "absent_in_focus_regions"


def build_interpretation(
    candidate: dict[str, str],
    direction: str,
    pattern: str,
    log2_ratio: float,
) -> str:
    return (
        f"{direction} copy-number contrast for candidate {candidate['module_id']} orthogroup "
        f"{candidate['orthogroup_id']} linked to {candidate['effector_class']} effectors "
        f"({pattern}; log2 mean-copy ratio={log2_ratio:.2f}); candidate contrast only."
    )


def format_species_counts(values: dict[str, int]) -> str:
    return ";".join(f"{species_id}:{count}" for species_id, count in sorted(values.items()))


def to_int(value: str) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def to_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def mean(values: list[int]) -> float:
    return sum(values) / len(values) if values else 0.0


def format_float(value: float) -> str:
    return f"{value:.2f}"


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--stable-orthogroups", type=Path, required=True)
    parser.add_argument("--species-metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ratio-threshold", type=float, default=1.5)
    parser.add_argument("--pseudocount", type=float, default=0.5)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = contrast_candidate_regions(
        candidates_path=args.candidates,
        stable_orthogroups_path=args.stable_orthogroups,
        species_metadata_path=args.species_metadata,
        output_path=args.output,
        ratio_threshold=args.ratio_threshold,
        pseudocount=args.pseudocount,
    )
    print(f"East Asia vs North America candidate contrast rows written: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
