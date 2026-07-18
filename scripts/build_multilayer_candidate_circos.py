#!/usr/bin/env python3
"""Build a validated 12-candidate table for the multilayer mechanism Circos."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path


EAST_ASIA = ("pden", "pmas", "ptab")
NORTH_AMERICA = ("plam", "plon", "ptae")
COPY_CONTEXT_SPECIES = EAST_ASIA + NORTH_AMERICA
CONTRASTS = (
    "resistant_pwn_vs_water",
    "susceptible_pwn_vs_water",
    "genotype_by_inoculum",
)
HEADERLESS_FIELDS = ["feature_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
FUNCTIONAL_SUPPORT_FIELDS = ("domain_support", "merops_support", "dbcan_support")


def safe_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def simes(values: list[float]) -> float:
    invalid = [value for value in values if math.isfinite(value) and not 0 <= value <= 1]
    if invalid:
        raise ValueError("P values must be within [0, 1]")
    finite = sorted(value for value in values if math.isfinite(value))
    if not finite:
        return math.nan
    n = len(finite)
    return min(1.0, min(n * value / rank for rank, value in enumerate(finite, 1)))


def benjamini_hochberg(values: list[float]) -> list[float]:
    invalid = [value for value in values if math.isfinite(value) and not 0 <= value <= 1]
    if invalid:
        raise ValueError("P values must be within [0, 1]")
    indexed = [(index, value) for index, value in enumerate(values) if math.isfinite(value)]
    ordered = sorted(indexed, key=lambda item: item[1])
    adjusted = [math.nan] * len(values)
    running = 1.0
    for rank_index in range(len(ordered) - 1, -1, -1):
        original_index, value = ordered[rank_index]
        rank = rank_index + 1
        running = min(running, value * len(values) / rank)
        adjusted[original_index] = min(1.0, running)
    return adjusted


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_deseq2(path: Path, interaction: bool = False) -> list[dict[str, object]]:
    with path.open(encoding="utf-8", newline="") as handle:
        first = handle.readline()
        handle.seek(0)
        has_header = first.split("\t", 1)[0] == "feature_id"
        reader = csv.DictReader(handle, delimiter="\t", fieldnames=None if has_header else HEADERLESS_FIELDS)
        rows: list[dict[str, object]] = []
        for row in reader:
            rows.append(
                {
                    "feature_id": row.get("feature_id", ""),
                    "log2_fold_change": safe_float(
                        row.get("interaction_log2_fold_change" if interaction else "log2FoldChange", "")
                    ),
                    "pvalue": safe_float(row.get("interaction_pvalue" if interaction else "pvalue", "")),
                }
            )
    return rows


def normalize_feature(identifier: str) -> str:
    return identifier.split("|", 1)[1] if "|" in identifier else identifier


def aggregate_contrast(
    candidate_ids: list[str],
    feature_rows: list[dict[str, object]],
    orthogroup_by_feature: dict[str, str],
) -> dict[str, dict[str, object]]:
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    candidate_set = set(candidate_ids)
    for row in feature_rows:
        feature = normalize_feature(str(row.get("feature_id", "")))
        orthogroup = orthogroup_by_feature.get(feature, "")
        if orthogroup in candidate_set:
            grouped[orthogroup].append({**row, "feature_id": feature})

    output: dict[str, dict[str, object]] = {}
    for orthogroup in candidate_ids:
        rows = grouped.get(orthogroup, [])
        finite_effects = [row for row in rows if math.isfinite(float(row["log2_fold_change"]))]
        finite_pvalues = [float(row["pvalue"]) for row in rows if math.isfinite(float(row["pvalue"]))]
        representative = None
        if finite_effects:
            representative = sorted(
                finite_effects,
                key=lambda row: (-abs(float(row["log2_fold_change"])), str(row["feature_id"])),
            )[0]
        output[orthogroup] = {
            "representative_feature": "" if representative is None else str(representative["feature_id"]),
            "representative_log2_fold_change": (
                math.nan if representative is None else float(representative["log2_fold_change"])
            ),
            "tested_member_count": len(rows),
            "finite_pvalue_member_count": len(finite_pvalues),
            "simes_pvalue": simes(finite_pvalues),
            "simes_padj": math.nan,
            "state": "unavailable",
        }
    return output


def assign_contrast_fdr(summaries: dict[str, dict[str, object]], candidate_ids: list[str]) -> None:
    adjusted = benjamini_hochberg([float(summaries[og]["simes_pvalue"]) for og in candidate_ids])
    for orthogroup, padj in zip(candidate_ids, adjusted):
        summary = summaries[orthogroup]
        summary["simes_padj"] = padj
        if not math.isfinite(padj):
            summary["state"] = "unavailable"
        elif padj < 0.05:
            summary["state"] = "significant"
        else:
            summary["state"] = "not_significant"


def regional_effect(rows: list[dict[str, str]], orthogroup_id: str) -> dict[str, object]:
    matched = [row for row in rows if row.get("orthogroup_id") == orthogroup_id]
    if not matched:
        return {"regional_log2_ratio": math.nan, "regional_direction": "unavailable"}
    values = {
        safe_float(row.get("log2_east_asia_vs_north_america_mean_copy_ratio", ""))
        for row in matched
    }
    finite = {value for value in values if math.isfinite(value)}
    if len(finite) != 1:
        raise ValueError(f"Conflicting regional effects for {orthogroup_id}: {sorted(finite)}")
    directions = {row.get("regional_copy_direction", "") for row in matched}
    if len(directions) != 1:
        raise ValueError(f"Conflicting regional directions for {orthogroup_id}: {sorted(directions)}")
    return {"regional_log2_ratio": finite.pop(), "regional_direction": directions.pop()}


def stable_membership(rows: list[dict[str, str]]) -> tuple[dict[str, str], dict[str, Counter[str]]]:
    orthogroup_by_feature: dict[str, str] = {}
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in rows:
        orthogroup = row.get("orthogroup_id", "")
        species = row.get("species_id", "")
        feature = normalize_feature(row.get("gene_id", ""))
        if species == "pmas" and feature:
            previous = orthogroup_by_feature.get(feature)
            if previous and previous != orthogroup:
                raise ValueError(f"Conflicting orthogroup assignment for {feature}: {previous}, {orthogroup}")
            orthogroup_by_feature[feature] = orthogroup
        if orthogroup and species:
            counts[orthogroup][species] += 1
    return orthogroup_by_feature, counts


def copy_effect(counts: Counter[str]) -> dict[str, object]:
    east_mean = sum(counts[species] for species in EAST_ASIA) / len(EAST_ASIA)
    north_mean = sum(counts[species] for species in NORTH_AMERICA) / len(NORTH_AMERICA)
    ratio = math.log2((east_mean + 0.5) / (north_mean + 0.5))
    threshold = math.log2(1.5)
    direction = "East_Asia_enriched" if ratio >= threshold else "North_America_enriched" if ratio <= -threshold else "balanced"
    return {
        "east_asia_mean_copy": east_mean,
        "north_america_mean_copy": north_mean,
        "regional_log2_ratio": ratio,
        "regional_direction": direction,
        "east_asia_species_counts": ";".join(f"{species}:{counts[species]}" for species in EAST_ASIA),
        "north_america_species_counts": ";".join(f"{species}:{counts[species]}" for species in NORTH_AMERICA),
    }


def center_species_copy_context(counts: Counter[str]) -> dict[str, object]:
    for species, count in counts.items():
        if count < 0:
            raise ValueError(f"Negative copy count for {species}: {count}")
        if species not in COPY_CONTEXT_SPECIES and count != 0:
            raise ValueError(f"Unknown species copy count for {species}: {count}")

    transformed = {
        species: math.log2(counts[species] + 0.5)
        for species in COPY_CONTEXT_SPECIES
    }
    mean_transformed = sum(transformed.values()) / len(COPY_CONTEXT_SPECIES)
    context: dict[str, object] = {
        f"copy_{species}": counts[species]
        for species in COPY_CONTEXT_SPECIES
    }
    context.update(
        {
            f"centered_log2_copy_{species}": transformed[species] - mean_transformed
            for species in COPY_CONTEXT_SPECIES
        }
    )
    context["copy_context_claim_boundary"] = (
        "descriptive within-orthogroup copy context; not formal gene-family expansion/contraction"
    )
    return context


def cross_species_state(
    candidate: dict[str, str],
    evidence_rows: list[dict[str, str]],
) -> str:
    concordant = int(candidate.get("nonproxy_concordant_species_count", "0") or 0)
    discordant = int(candidate.get("nonproxy_discordant_species_count", "0") or 0)
    if concordant and discordant:
        return "heterogeneous"
    if concordant:
        return "concordant"
    if discordant:
        return "discordant"
    orthogroup = candidate.get("orthogroup_id", "")
    tested = any(
        row.get("orthogroup_id") == orthogroup
        and int(row.get("finite_padj_gene_count", "0") or 0) > 0
        for row in evidence_rows
    )
    if tested:
        return "nonsignificant"
    return "unavailable"


def truthy(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes"}


def format_value(value: object) -> str:
    if isinstance(value, float):
        return "" if not math.isfinite(value) else f"{value:.12g}"
    return str(value)


def validate_interaction(
    summaries: dict[str, dict[str, object]],
    existing_rows: list[dict[str, str]],
    candidate_ids: list[str],
) -> None:
    existing = {row["orthogroup_id"]: row for row in existing_rows if row.get("orthogroup_id") in candidate_ids}
    for orthogroup in candidate_ids:
        if orthogroup not in existing:
            raise ValueError(f"Missing interaction summary for {orthogroup}")
        current = summaries[orthogroup]
        prior = existing[orthogroup]
        for field in ("tested_member_count", "finite_pvalue_member_count"):
            if int(current[field]) != int(prior[field]):
                raise ValueError(f"Interaction {field} disagreement for {orthogroup}")
        prior_simes = safe_float(prior.get("simes_pvalue", ""))
        current_simes = float(current["simes_pvalue"])
        if math.isfinite(prior_simes) != math.isfinite(current_simes) or (
            math.isfinite(prior_simes) and not math.isclose(prior_simes, current_simes, rel_tol=1e-8, abs_tol=1e-12)
        ):
            raise ValueError(f"Interaction Simes disagreement for {orthogroup}")
        prior_effect = safe_float(prior.get("max_abs_interaction_lfc", ""))
        current_effect = abs(float(current["representative_log2_fold_change"]))
        if math.isfinite(prior_effect) != math.isfinite(current_effect) or (
            math.isfinite(prior_effect)
            and not math.isclose(prior_effect, current_effect, rel_tol=1e-8, abs_tol=1e-12)
        ):
            raise ValueError(f"Interaction representative effect disagreement for {orthogroup}")


def build_effector_class_context(
    effectors: list[dict[str, str]],
    audit_rows: list[dict[str, str]],
    edge_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    allowed_classes = list(
        dict.fromkeys(
            row.get("source_id", "")
            for row in edge_rows
            if row.get("edge_type") == "functional_prior" and row.get("source_id", "")
        )
    )
    if not allowed_classes:
        raise ValueError("No effector classes with functional-prior links were found")

    audit_by_id: dict[str, dict[str, str]] = {}
    for row in audit_rows:
        protein_id = row.get("protein_id", "")
        if not protein_id:
            continue
        if protein_id in audit_by_id:
            raise ValueError(f"Duplicate secretome-audit protein identifier: {protein_id}")
        audit_by_id[protein_id] = row

    seen_effectors: set[str] = set()
    counts = {
        effector_class: {"standard": 0, "evaluated": 0, "supported": 0}
        for effector_class in allowed_classes
    }
    for row in effectors:
        protein_id = row.get("protein_id", "")
        if not protein_id:
            continue
        if protein_id in seen_effectors:
            raise ValueError(f"Duplicate effector protein identifier: {protein_id}")
        seen_effectors.add(protein_id)
        if protein_id not in audit_by_id:
            raise ValueError(f"Effector protein {protein_id} is missing from the secretome audit")
        effector_class = row.get("effector_class", "")
        if effector_class not in counts:
            continue
        audit = audit_by_id[protein_id]
        if audit.get("standard_secreted_soluble_status") != "pass":
            continue
        counts[effector_class]["standard"] += 1
        support_states = tuple(audit.get(field, "not_evaluated") for field in FUNCTIONAL_SUPPORT_FIELDS)
        if any(state != "not_evaluated" for state in support_states):
            counts[effector_class]["evaluated"] += 1
        if any(state == "yes" for state in support_states):
            counts[effector_class]["supported"] += 1

    output: list[dict[str, str]] = []
    for display_order, effector_class in enumerate(allowed_classes, 1):
        standard = counts[effector_class]["standard"]
        evaluated = counts[effector_class]["evaluated"]
        supported = counts[effector_class]["supported"]
        fraction = supported / evaluated if evaluated else math.nan
        state = "supported" if supported else "evaluated_none" if evaluated else "unavailable"
        output.append(
            {
                "display_order": str(display_order),
                "effector_class": effector_class,
                "standard_candidate_count": str(standard),
                "functional_evaluated_count": str(evaluated),
                "functional_supported_count": str(supported),
                "functional_supported_fraction": format_value(fraction),
                "functional_annotation_state": state,
                "claim_boundary": "class-level functional annotation context; not direct host targeting",
            }
        )
    return output


def build_candidate_layers(
    candidates: list[dict[str, str]],
    stable_rows: list[dict[str, str]],
    regional_rows: list[dict[str, str]],
    contrasts: dict[str, list[dict[str, object]]],
    interaction_summary: list[dict[str, str]],
    cross_species_evidence: list[dict[str, str]],
) -> list[dict[str, str]]:
    candidate_ids = [row.get("orthogroup_id", "") for row in candidates]
    if len(candidate_ids) != 12 or len(set(candidate_ids)) != 12 or any(not value for value in candidate_ids):
        raise ValueError("Expected exactly 12 unique candidate orthogroups")
    orthogroup_by_feature, all_counts = stable_membership(stable_rows)
    summaries = {
        contrast: aggregate_contrast(candidate_ids, rows, orthogroup_by_feature)
        for contrast, rows in contrasts.items()
    }
    if set(summaries) != set(CONTRASTS):
        raise ValueError(f"Expected contrasts {CONTRASTS}, observed {sorted(summaries)}")
    for contrast in CONTRASTS:
        assign_contrast_fdr(summaries[contrast], candidate_ids)
    validate_interaction(summaries["genotype_by_inoculum"], interaction_summary, candidate_ids)

    output: list[dict[str, str]] = []
    for display_order, candidate in enumerate(candidates, 1):
        orthogroup = candidate["orthogroup_id"]
        copy = copy_effect(all_counts[orthogroup])
        registered = regional_effect(regional_rows, orthogroup)
        if math.isfinite(float(registered["regional_log2_ratio"])) and not math.isclose(
            float(registered["regional_log2_ratio"]), float(copy["regional_log2_ratio"]), abs_tol=0.02
        ):
            raise ValueError(f"Regional contrast disagreement for {orthogroup}")
        copy_context_counts = Counter(
            {species: all_counts[orthogroup][species] for species in COPY_CONTEXT_SPECIES}
        )
        row: dict[str, object] = {
            "display_order": display_order,
            "orthogroup_id": orthogroup,
            "manuscript_role": candidate.get("manuscript_role", ""),
            "priority_class": candidate.get("overall_priority_class", ""),
            "mechanism_axes": candidate.get("mechanism_axes", ""),
            "module_ids": candidate.get("module_ids", ""),
            "mapped_symbols": candidate.get("mapped_symbols", ""),
            **copy,
            **center_species_copy_context(copy_context_counts),
            "cross_species_state": cross_species_state(candidate, cross_species_evidence),
            "network_support": "yes" if truthy(candidate.get("experimental_or_database_supported", "")) else "no",
            "structure_support": "yes" if int(candidate.get("structurally_supported_effector_count", "0") or 0) > 0 else "no",
            "sequence_support": "yes" if candidate.get("apoplastic_cell_wall_sequence_support_class") == "tier1_candidate_substrate" else "no",
            "claim_boundary": "candidate synthesis only; does not establish direct effector targeting, causal resistance, or geographic adaptation",
        }
        for contrast in CONTRASTS:
            for key, value in summaries[contrast][orthogroup].items():
                row[f"{contrast}_{key}"] = value
        output.append({key: format_value(value) for key, value in row.items()})
    return output


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--stable-genes", type=Path, required=True)
    parser.add_argument("--regional-contrasts", type=Path, required=True)
    parser.add_argument("--resistant", type=Path, required=True)
    parser.add_argument("--susceptible", type=Path, required=True)
    parser.add_argument("--interaction", type=Path, required=True)
    parser.add_argument("--interaction-summary", type=Path, required=True)
    parser.add_argument("--cross-species-validation", type=Path, required=True)
    parser.add_argument("--effectors", type=Path, required=True)
    parser.add_argument("--secretome-audit", type=Path, required=True)
    parser.add_argument("--edges", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--effector-context-output", type=Path, required=True)
    args = parser.parse_args()
    rows = build_candidate_layers(
        read_tsv(args.candidates),
        read_tsv(args.stable_genes),
        read_tsv(args.regional_contrasts),
        {
            "resistant_pwn_vs_water": read_deseq2(args.resistant),
            "susceptible_pwn_vs_water": read_deseq2(args.susceptible),
            "genotype_by_inoculum": read_deseq2(args.interaction, interaction=True),
        },
        read_tsv(args.interaction_summary),
        read_tsv(args.cross_species_validation),
    )
    write_tsv(args.output, rows)
    effector_context = build_effector_class_context(
        read_tsv(args.effectors),
        read_tsv(args.secretome_audit),
        read_tsv(args.edges),
    )
    write_tsv(args.effector_context_output, effector_context)
    significant = {
        contrast: sum(row[f"{contrast}_state"] == "significant" for row in rows)
        for contrast in CONTRASTS
    }
    print(
        f"Wrote {len(rows)} multilayer candidates and {len(effector_context)} effector classes; "
        f"significant={significant}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
