#!/usr/bin/env python3
"""Create the P. massoniana mechanism competition shortlist and report."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


FIELDS = [
    "candidate_id", "candidate_level", "orthogroup_id", "feature_ids", "module_ids",
    "mechanism_axes", "evidence_tier", "selection_basis", "pmas_simes_padj",
    "significant_member_count", "max_abs_interaction_lfc", "cross_species_support",
    "network_followup", "claim_ceiling",
]
CLAIM = "Candidate mechanism only; does not establish direct effector targeting, causal resistance, coevolution, or a regional difference"


def _read(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _significant(row: dict[str, str]) -> bool:
    return row.get("competition_status") == "significant"


def _route_statuses(
    ranked: list[dict[str, str]], ora: list[dict[str, str]], axes: list[dict[str, str]]
) -> dict[str, str]:
    ranked_by_id = {row["set_id"]: row for row in ranked if row.get("set_level") == "mechanism_axis"}
    ora_by_id = {row["set_id"]: row for row in ora if row.get("set_level") == "mechanism_axis"}
    statuses = {}
    for axis in axes:
        axis_id = axis["mechanism_axis"]
        ranked_row, ora_row = ranked_by_id.get(axis_id, {}), ora_by_id.get(axis_id, {})
        if not axis.get("host_modules", "").strip() or ranked_row.get("competition_status") == "not_testable":
            statuses[axis_id] = "not_testable"
        elif _significant(ranked_row) and _significant(ora_row):
            statuses[axis_id] = "supported"
        elif _significant(ranked_row) or _significant(ora_row):
            statuses[axis_id] = "suggestive"
        else:
            statuses[axis_id] = "not_supported"
    return statuses


def build_report(
    features_path: Path,
    ranked_path: Path,
    ora_path: Path,
    orthogroups_path: Path,
    validation_path: Path,
    axes_path: Path,
    shortlist_path: Path,
    report_path: Path,
) -> tuple[list[dict[str, str]], dict[str, str]]:
    features, ranked, ora = _read(features_path), _read(ranked_path), _read(ora_path)
    orthogroups, validation, axes = _read(orthogroups_path), _read(validation_path), _read(axes_path)
    statuses = _route_statuses(ranked, ora, axes)

    exact_concordant: set[str] = set()
    for row in validation:
        if row.get("can_support_tier_a") == "yes" and row.get("evidence_classification") == "directionally_concordant":
            exact_concordant.add(row["orthogroup_id"])

    candidates: list[dict[str, str]] = []
    covered_features: set[str] = set()
    for row in orthogroups:
        significant_count = int(row.get("significant_member_count", "0") or 0)
        if significant_count == 0 or not row.get("mechanism_axes", "").strip():
            continue
        orthogroup = row["orthogroup_id"]
        tier = "Tier A" if orthogroup in exact_concordant else "Tier B"
        significant_features = row.get("significant_feature_ids", "")
        covered_features.update(part for part in significant_features.split(";") if part)
        candidates.append({
            "candidate_id": orthogroup,
            "candidate_level": "orthogroup",
            "orthogroup_id": orthogroup,
            "feature_ids": significant_features,
            "module_ids": row.get("module_ids", ""),
            "mechanism_axes": row.get("mechanism_axes", ""),
            "evidence_tier": tier,
            "selection_basis": "significant_pmas_orthogroup_with_nonproxy_cross_species_support" if tier == "Tier A" else "significant_pmas_orthogroup_with_mechanism_annotation",
            "pmas_simes_padj": row.get("simes_padj", ""),
            "significant_member_count": str(significant_count),
            "max_abs_interaction_lfc": row.get("max_abs_interaction_lfc", ""),
            "cross_species_support": "nonproxy_directionally_concordant" if tier == "Tier A" else "incomplete_or_proxy_only",
            "network_followup": "yes",
            "claim_ceiling": CLAIM,
        })

    for row in features:
        if row.get("foreground_850") != "yes" or row["feature_id"] in covered_features:
            continue
        if row.get("orthogroup_id", "").strip() and row.get("mechanism_axes", "").strip():
            continue
        candidates.append({
            "candidate_id": row["feature_id"],
            "candidate_level": "feature",
            "orthogroup_id": row.get("orthogroup_id", ""),
            "feature_ids": row["feature_id"],
            "module_ids": row.get("module_ids", ""),
            "mechanism_axes": row.get("mechanism_axes", ""),
            "evidence_tier": "Tier C",
            "selection_basis": "significant_pmas_feature_with_incomplete_orthogroup_or_mechanism_mapping",
            "pmas_simes_padj": "",
            "significant_member_count": "1",
            "max_abs_interaction_lfc": row.get("interaction_log2_fold_change", ""),
            "cross_species_support": "unavailable",
            "network_followup": "no",
            "claim_ceiling": CLAIM,
        })

    tier_order = {"Tier A": 0, "Tier B": 1, "Tier C": 2}
    candidates.sort(key=lambda row: (tier_order[row["evidence_tier"]], row["candidate_id"]))
    Path(shortlist_path).parent.mkdir(parents=True, exist_ok=True)
    with Path(shortlist_path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(candidates)

    ranked_by_id = {row["set_id"]: row for row in ranked if row.get("set_level") == "mechanism_axis"}
    ora_by_id = {row["set_id"]: row for row in ora if row.get("set_level") == "mechanism_axis"}
    mapped = sum(bool(row.get("orthogroup_id", "").strip()) for row in features)
    finite_padj = sum(row.get("finite_padj") == "yes" for row in features)
    foreground = sum(row.get("foreground_850") == "yes" for row in features)
    tiers = Counter(row["evidence_tier"] for row in candidates)
    validation_classes = Counter(row.get("evidence_classification", "") for row in validation)
    lines = [
        "# P. massoniana mechanism competition",
        "",
        "## Analysis universe",
        "",
        f"The analysis contains {len(features):,} interaction-tested features; {mapped:,} map to a stable orthogroup, {len(features)-mapped:,} are unmapped, {finite_padj:,} have finite adjusted P values, and {foreground:,} meet the strict foreground definition.",
        "",
        "## Route competition",
        "",
        "| Route | Outcome | Ranked NES | Ranked FDR | ORA odds ratio | ORA FDR |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for axis in axes:
        axis_id = axis["mechanism_axis"]
        rank, over = ranked_by_id.get(axis_id, {}), ora_by_id.get(axis_id, {})
        lines.append(f"| {axis.get('route_label', axis_id)}: {axis_id} | {statuses[axis_id]} | {rank.get('nes','NA') or 'NA'} | {rank.get('padj','NA') or 'NA'} | {over.get('odds_ratio','NA') or 'NA'} | {over.get('padj','NA') or 'NA'} |")
    lines += [
        "",
        "Route outcomes are categorical evidence decisions, not a composite score. A supported route requires significant ranked enrichment plus a concordant secondary layer; a route without configured direct host modules is not testable.",
        "",
        "## Candidate shortlist",
        "",
        f"The shortlist contains {len(candidates):,} candidates: {tiers['Tier A']:,} Tier A, {tiers['Tier B']:,} Tier B, and {tiers['Tier C']:,} Tier C. Only Tier A and Tier B proceed to predicted interolog/PPI follow-up.",
        "",
        "## Cross-species evidence",
        "",
        f"Validation rows: {len(validation):,}; directionally concordant: {validation_classes['directionally_concordant']:,}; directionally discordant: {validation_classes['directionally_discordant']:,}; not significant: {validation_classes['not_significant']:,}; unavailable: {validation_classes['unavailable']:,}.",
        "",
        "P. strobus mappings labeled as homology proxies cannot support Tier A. Unavailable evidence is retained as unavailable rather than treated as a negative result.",
        "P. densiflora Tier A support uses same-species targeted homology to the current proteome; it is stronger than a cross-species proxy but is not an exact transcript-to-gene identifier match.",
        "P. rigida x P. taeda expression was quantified against P. taeda CDS and is retained as a parental-reference proxy; it can support reproducibility but cannot independently create Tier A evidence.",
        "P. thunbergii expression was quantified against P. taeda CDS and is retained as a surrogate-reference proxy; it can support infection-response reproducibility but cannot independently create Tier A evidence.",
        "",
        "## Claim ceiling",
        "",
        "This analysis identifies statistically coordinated genotype-dependent infection-response candidates. It does not establish direct effector targeting, causal resistance, host-parasite coevolution, or a general East Asia versus North America difference.",
        "",
        "## Negative results and gaps",
        "",
        "Routes without ranked support are not rescued by candidate-level interpretation. Cross-species datasets with no matching candidate orthogroup define a validation gap, not evidence against that candidate.",
    ]
    Path(report_path).parent.mkdir(parents=True, exist_ok=True)
    Path(report_path).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return candidates, statuses


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ["features", "ranked", "ora", "orthogroups", "validation", "axes", "shortlist", "report"]:
        parser.add_argument(f"--{argument}", type=Path, required=True)
    args = parser.parse_args()
    rows, statuses = build_report(args.features, args.ranked, args.ora, args.orthogroups, args.validation, args.axes, args.shortlist, args.report)
    print(f"Wrote {len(rows)} candidates; route outcomes: {statuses}")


if __name__ == "__main__":
    main()
