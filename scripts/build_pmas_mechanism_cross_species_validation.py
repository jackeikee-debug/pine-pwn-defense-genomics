#!/usr/bin/env python3
"""Build targeted cross-species expression validation for Pmas candidate OGs."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


FIELDS = [
    "orthogroup_id", "mechanism_axes", "selection_reasons", "species", "contrast",
    "mapping_evidence", "can_support_tier_a", "source_feature_count",
    "cross_species_log2_fold_change", "cross_species_padj", "evidence_classification",
    "pmas_direction_status", "source_file", "claim_ceiling",
]

PDEN_CONTRASTS = [
    ("pathogen_associated", "primary_log2_fold_change", "primary_padj"),
    ("bxyl_vs_water", "bxyl_vs_water_log2_fold_change", "bxyl_vs_water_padj"),
    ("bthai_vs_water", "bthai_vs_water_log2_fold_change", "bthai_vs_water_padj"),
]
PSTROBUS_CONTRASTS = [
    ("two_week_vs_zero", "two_week_log2_fold_change", "two_week_padj"),
    ("four_week_vs_zero", "four_week_log2_fold_change", "four_week_padj"),
]
CLAIM = "Expression supports infection response only; it is not direct effector targeting, causal resistance, or regional evidence"


def _read(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _number(value: str) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _format(value: float | None) -> str:
    return "" if value is None else f"{value:.12g}"


def _classification(direction: str, effect: float | None, padj: float | None) -> str:
    if padj is None or padj >= 0.05 or effect is None:
        return "not_significant"
    if direction not in {"positive", "negative"} or effect == 0:
        return "direction_not_comparable"
    same = (direction == "positive" and effect > 0) or (direction == "negative" and effect < 0)
    return "directionally_concordant" if same else "directionally_discordant"


def _best(rows: list[dict[str, str]], effect_col: str, padj_col: str) -> tuple[float | None, float | None]:
    candidates = [(_number(row.get(effect_col, "")), _number(row.get(padj_col, ""))) for row in rows]
    finite = [(effect, padj) for effect, padj in candidates if effect is not None and padj is not None]
    if finite:
        return min(finite, key=lambda pair: (pair[1], -abs(pair[0])))
    effects = [effect for effect, _ in candidates if effect is not None]
    return (max(effects, key=abs), None) if effects else (None, None)


def build_cross_species_validation(
    features_path: Path,
    ranked_path: Path,
    orthogroups_path: Path,
    pden_path: Path,
    pstrobus_path: Path,
    output_path: Path,
) -> list[dict[str, str]]:
    feature_to_og = {row["feature_id"]: row.get("orthogroup_id", "") for row in _read(features_path)}
    og_rows = {row["orthogroup_id"]: row for row in _read(orthogroups_path)}
    reasons: dict[str, set[str]] = defaultdict(set)
    for row in _read(ranked_path):
        if row.get("set_level") != "mechanism_axis" or row.get("competition_status") != "significant":
            continue
        for feature in row.get("leading_edge_features", "").split(";"):
            orthogroup = feature_to_og.get(feature, "")
            if orthogroup:
                reasons[orthogroup].add(f"significant_gsea_leading_edge:{row['set_id']}")
    for orthogroup, row in og_rows.items():
        simes = _number(row.get("simes_padj", ""))
        if simes is not None and simes < 0.05 and int(row.get("significant_member_count", "0") or 0) > 0:
            reasons[orthogroup].add("simes_fdr_and_foreground")

    pden_by_og: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in _read(pden_path):
        for orthogroup in row.get("orthogroup_ids", "").split(";"):
            if orthogroup:
                pden_by_og[orthogroup].append(row)
    pstrobus_by_og: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in _read(pstrobus_path):
        if row.get("orthogroup_id"):
            pstrobus_by_og[row["orthogroup_id"]].append(row)

    output: list[dict[str, str]] = []
    for orthogroup in sorted(reasons):
        og = og_rows[orthogroup]
        common = {
            "orthogroup_id": orthogroup,
            "mechanism_axes": og.get("mechanism_axes", ""),
            "selection_reasons": ";".join(sorted(reasons[orthogroup])),
            "pmas_direction_status": og.get("direction_status", "not_estimable"),
            "claim_ceiling": CLAIM,
        }
        for species, source, source_rows, contrasts, mapping, tier_a in [
            ("Pinus densiflora", Path(pden_path).name, pden_by_og.get(orthogroup, []), PDEN_CONTRASTS, "exact_or_curated_orthogroup", "yes"),
            ("Pinus strobus", Path(pstrobus_path).name, pstrobus_by_og.get(orthogroup, []), PSTROBUS_CONTRASTS, "homology_proxy", "no"),
        ]:
            for contrast, effect_col, padj_col in contrasts:
                effect, padj = _best(source_rows, effect_col, padj_col)
                available = bool(source_rows)
                output.append({
                    **common,
                    "species": species,
                    "contrast": contrast,
                    "mapping_evidence": mapping if available else "unavailable",
                    "can_support_tier_a": tier_a if available else "no",
                    "source_feature_count": str(len(source_rows)),
                    "cross_species_log2_fold_change": _format(effect),
                    "cross_species_padj": _format(padj),
                    "evidence_classification": _classification(common["pmas_direction_status"], effect, padj) if available else "unavailable",
                    "source_file": source,
                })

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with Path(output_path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--ranked-enrichment", type=Path, required=True)
    parser.add_argument("--orthogroups", type=Path, required=True)
    parser.add_argument("--pden", type=Path, required=True)
    parser.add_argument("--pstrobus", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = build_cross_species_validation(args.features, args.ranked_enrichment, args.orthogroups, args.pden, args.pstrobus, args.output)
    print(f"Wrote {len(rows)} cross-species validation rows to {args.output}")


if __name__ == "__main__":
    main()
