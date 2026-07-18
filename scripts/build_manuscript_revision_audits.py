#!/usr/bin/env python3
"""Build traceable evidence-audit tables used by the revised manuscript."""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter
from pathlib import Path


OUTPUT_NAMES = {
    "expression_dataset_summary": "manuscript_expression_dataset_summary.tsv",
    "candidate_attrition": "manuscript_candidate_attrition.tsv",
    "priority_rule_matrix": "manuscript_priority_rule_matrix.tsv",
    "expression_member_audit": "manuscript_expression_member_audit.tsv",
    "biochemical_test_candidates": "manuscript_biochemical_test_candidates.tsv",
}


def read_tsv(path: Path, fieldnames: list[str] | None = None) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        if fieldnames is None:
            return list(csv.DictReader(handle, delimiter="\t"))
        return list(csv.DictReader(handle, delimiter="\t", fieldnames=fieldnames))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def safe_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def clean_number(value: object) -> str:
    number = safe_float(value)
    return "" if not math.isfinite(number) else f"{number:.12g}"


def normalize_feature(value: str) -> str:
    return value.split("|", 1)[1] if "|" in value else value


def count_samples(path: Path) -> int:
    return len(read_tsv(path))


def build_dataset_summary(root: Path) -> list[dict[str, object]]:
    return [
        {
            "species_id": "pmas",
            "species": "Pinus massoniana",
            "bioproject": "PRJNA382462",
            "sra_study": "SRP103562",
            "sample_count": count_samples(root / "config/effector_target_network_pmas_1dpi_primary_sample_sheet.tsv"),
            "tissue_and_time": "stem tissue; 1 day post inoculation",
            "experimental_groups": "resistant-PWN; resistant-water; susceptible-PWN; susceptible-water; n=3 each",
            "reference": "GigaDB 102688 P. massoniana CDS; 76912 unique targets after exact-duplicate removal",
            "quantification": "FastQC 0.12.1; Salmon 1.10.1; 39-51 nt single-end unstranded reads; sequence-bias correction; mean mapping 58.79%",
            "gene_level_method": "tximport 1.30.0; countsFromAbundance=lengthScaledTPM; identity tx2gene",
            "differential_expression": "DESeq2 1.42.0; Wald tests; Benjamini-Hochberg adjustment",
            "design_formula": "~ phenotype + inoculum + phenotype:inoculum",
            "contrasts": "resistant PWN vs water; susceptible PWN vs water; genotype-by-inoculum interaction",
            "evidence_scope": "controlled_infection_and_genotype_by_inoculum_response",
            "limitations": "A defensible balanced 1-dpi subset was used because several public run labels outside this subset were ambiguous.",
        },
        {
            "species_id": "pden",
            "species": "Pinus densiflora",
            "bioproject": "PRJNA496563",
            "sra_study": "SRP165817",
            "sample_count": count_samples(root / "config/effector_target_network_pden_deseq2_sample_sheet.tsv"),
            "tissue_and_time": "cambium; 4 weeks post inoculation",
            "experimental_groups": "B. xylophilus; B. thailandae; water; n=3 each",
            "reference": "published P. densiflora Trinity transcriptome; transcript-to-gene collapse",
            "quantification": "FastQC 0.12.1; Galaxy Salmon 1.10.1+galaxy4; 101 nt paired-end reads; sequence- and GC-bias correction; mean mapping 91.53%",
            "gene_level_method": "transcript-to-gene aggregation using the published Trinity identifier structure",
            "differential_expression": "Galaxy DESeq2 wrapper 2.11.40.8+galaxy1; Benjamini-Hochberg adjustment",
            "design_formula": "~ condition",
            "contrasts": "B. xylophilus vs B. thailandae; B. xylophilus vs water; B. thailandae vs water",
            "evidence_scope": "pathogenic_species_contrast_and_inoculation_response",
            "limitations": "No explicit count prefilter was applied; DESeq2 independent filtering remained enabled. B. xylophilus vs B. thailandae is a nematode-species contrast and is not by itself a pathogenicity test.",
        },
        {
            "species_id": "pstrobus",
            "species": "Pinus strobus",
            "bioproject": "PRJNA1127083",
            "sra_study": "SRP515438",
            "sample_count": count_samples(root / "config/effector_target_network_pstrobus_deseq2_sample_sheet.tsv"),
            "tissue_and_time": "post-inoculation series; 0, 2, and 4 weeks; n=3 each",
            "experimental_groups": "PWN-inoculated temporal series without a contemporaneous mock control",
            "reference": "P. strobus GIIE transcript reference",
            "quantification": "Salmon 1.10.1/1.10.2; source recorded per sample",
            "gene_level_method": "tximport; countsFromAbundance=lengthScaledTPM",
            "differential_expression": "DESeq2; timepoint design; Wald tests; Benjamini-Hochberg adjustment",
            "design_formula": "~ timepoint",
            "contrasts": "2 weeks vs 0 weeks; 4 weeks vs 0 weeks",
            "evidence_scope": "post_inoculation_temporal_response_without_mock_control",
            "limitations": "Changes cannot be separated from development, handling, or inoculation effects; orthogroup support is a consensus-homology projection.",
        },
    ]


def unique_orthogroups(path: Path) -> set[str]:
    return {row["orthogroup_id"] for row in read_tsv(path) if row.get("orthogroup_id")}


def linked_orthogroups(path: Path) -> set[str]:
    output: set[str] = set()
    for row in read_tsv(path):
        output.update(item for item in row.get("host_orthogroups", "").split(";") if item)
    return output


def build_attrition(root: Path) -> list[dict[str, object]]:
    stages = [
        ("stable_orthology", unique_orthogroups(root / "results/tables/stable_orthogroup_genes.tsv"),
         "Orthogroups retained in the stable six-pine OrthoFinder comparison."),
        ("defense_annotation", unique_orthogroups(root / "results/tables/orthogroup_defense_matrix.tsv"),
         "At least one keyword-supported defense-module annotation in the six-pine panel."),
        ("effector_informed_module_link", linked_orthogroups(root / "results/tables/effector_host_module_links.tsv"),
         "Belongs to a host module linked to a predicted secreted-protein functional class; this is not a physical target call."),
        ("integrated_mechanism_review", unique_orthogroups(root / "results/tables/mechanism_host_candidate_priority.tsv"),
         "Entered explicit expression, network, compartment, structure, and sequence evidence review."),
        ("core_manuscript_candidates", unique_orthogroups(root / "results/tables/manuscript_mechanism_candidate_summary.tsv"),
         "Overall priority was mechanism_leading or mechanism_supported; context-only candidates were excluded."),
    ]
    rows = []
    previous = None
    for order, (stage, values, rule) in enumerate(stages, 1):
        count = len(values)
        rows.append({
            "stage_order": order,
            "stage": stage,
            "unit": "orthogroup",
            "retained_count": count,
            "retention_from_previous": "" if previous is None or previous == 0 else f"{count / previous:.6f}",
            "selection_rule": rule,
            "claim_boundary": "Counts describe computational filtering, not validated resistance genes or direct effector targets.",
        })
        previous = count
    return rows


def build_priority_rules(root: Path) -> list[dict[str, object]]:
    priorities = read_tsv(root / "results/tables/mechanism_host_candidate_priority.tsv")
    counts = Counter(row["overall_priority_class"] for row in priorities)
    definitions = [
        ("mechanism_leading", "yes", "Convergent evidence from the focal P. massoniana analysis plus at least one independent evidence layer; conflicting contrasts remain visible."),
        ("mechanism_supported", "yes", "Focal expression evidence plus qualified cross-species, network, compartment, or sequence support, without meeting the leading tier."),
        ("mechanism_context", "no", "Biologically relevant context but insufficient independent directional or integrated evidence for the core set."),
    ]
    return [{
        "priority_class": priority,
        "orthogroup_count": counts[priority],
        "included_in_core_12": included,
        "operational_rule": rule,
        "direct_target_claim_allowed": "no",
        "causal_resistance_claim_allowed": "no",
    } for priority, included, rule in definitions]


def pmas_feature_map(root: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for row in read_tsv(root / "results/tables/stable_orthogroup_genes.tsv"):
        if row.get("species_id") == "pmas":
            mapping[normalize_feature(row.get("gene_id", ""))] = row["orthogroup_id"]
    return mapping


def read_pmas_result(path: Path, interaction: bool) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        has_header = handle.readline().split("\t", 1)[0] == "feature_id"
    headerless_fields = [
        "feature_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"
    ]
    rows = read_tsv(path, None if has_header else headerless_fields)
    lfc_key = "interaction_log2_fold_change" if interaction else "log2FoldChange"
    p_key = "interaction_pvalue" if interaction else "pvalue"
    padj_key = "interaction_padj" if interaction else "padj"
    return [{
        "feature_id": normalize_feature(row.get("feature_id", "")),
        "log2_fold_change": clean_number(row.get(lfc_key, "")),
        "pvalue": clean_number(row.get(p_key, "")),
        "padj": clean_number(row.get(padj_key, "")),
    } for row in rows]


def build_member_audit(root: Path) -> list[dict[str, object]]:
    core_rows = read_tsv(root / "results/tables/manuscript_multilayer_circos_candidates.tsv")
    core = {row["orthogroup_id"]: row for row in core_rows}
    feature_to_og = pmas_feature_map(root)
    contrasts = [
        ("resistant_pwn_vs_water", root / "results/tables/effector_target_network_pmas_GigaDB102688_DESeq2_resistant_pwn_vs_water.tsv", False),
        ("susceptible_pwn_vs_water", root / "results/tables/effector_target_network_pmas_GigaDB102688_DESeq2_susceptible_pwn_vs_water.tsv", False),
        ("genotype_by_inoculum", root / "results/tables/effector_target_network_pmas_GigaDB102688_DESeq2_genotype_by_inoculum.tsv", True),
    ]
    output: list[dict[str, object]] = []
    for contrast, path, interaction in contrasts:
        representative_key = f"{contrast}_representative_feature"
        for feature in read_pmas_result(path, interaction):
            orthogroup = feature_to_og.get(feature["feature_id"], "")
            if orthogroup not in core:
                continue
            summary = core[orthogroup]
            displayed = feature["feature_id"] == summary.get(representative_key, "")
            output.append({
                "orthogroup_id": orthogroup,
                "species": "Pinus massoniana",
                "contrast": contrast,
                "evidence_unit": "orthogroup_member",
                "feature_id": feature["feature_id"],
                "is_displayed_member": "yes" if displayed else "no",
                "log2_fold_change": feature["log2_fold_change"],
                "raw_pvalue": feature["pvalue"],
                "whole_transcriptome_padj": feature["padj"],
                "whole_transcriptome_test_status": (
                    "finite_adjusted_pvalue" if feature["padj"] else
                    "raw_pvalue_only_after_independent_filtering" if feature["pvalue"] else
                    "not_testable"
                ),
                "contributes_to_candidate_simes_test": "yes" if feature["pvalue"] else "no",
                "candidate_set_simes_pvalue": summary.get(f"{contrast}_simes_pvalue", ""),
                "candidate_set_simes_padj": summary.get(f"{contrast}_simes_padj", ""),
                "candidate_set_fdr_scope": "12 displayed orthogroups only",
                "directional_classification": "",
                "contrast_summary_status": "contrast_specific_opposite_directions" if orthogroup == "OG0020287" else "not_flagged",
                "source_file": path.name,
            })

    cross_species = read_tsv(root / "results/tables/pmas_mechanism_cross_species_validation.tsv")
    for row in cross_species:
        orthogroup = row.get("orthogroup_id", "")
        if orthogroup not in core:
            continue
        output.append({
            "orthogroup_id": orthogroup,
            "species": row.get("species", ""),
            "contrast": row.get("contrast", ""),
            "evidence_unit": "mapped_gene_set_summary",
            "feature_id": "",
            "is_displayed_member": "no",
            "log2_fold_change": row.get("cross_species_log2_fold_change", ""),
            "raw_pvalue": "",
            "whole_transcriptome_padj": row.get("cross_species_padj", ""),
            "whole_transcriptome_test_status": (
                "finite_adjusted_pvalue" if row.get("cross_species_padj", "") else "not_testable"
            ),
            "contributes_to_candidate_simes_test": "no",
            "candidate_set_simes_pvalue": "",
            "candidate_set_simes_padj": "",
            "candidate_set_fdr_scope": "not_applicable",
            "directional_classification": row.get("evidence_classification", ""),
            "contrast_summary_status": "contrast_specific_opposite_directions" if orthogroup == "OG0020287" else "not_flagged",
            "source_file": row.get("source_file", ""),
        })
    return sorted(output, key=lambda row: (str(row["orthogroup_id"]), str(row["species"]), str(row["contrast"]), str(row["feature_id"])))


def build_biochemical_candidates(root: Path) -> list[dict[str, object]]:
    core = unique_orthogroups(root / "results/tables/manuscript_mechanism_candidate_summary.tsv")
    rows = read_tsv(root / "results/tables/apoplastic_cell_wall_multinode_independence_audit.tsv")
    output = []
    for row in rows:
        if not row.get("sequence_resolved_tier", "").startswith("Tier 1"):
            continue
        included = row["orthogroup_id"] in core
        output.append({
            "orthogroup_id": row["orthogroup_id"],
            "effector_id": row.get("effector_id", ""),
            "host_family": row.get("host_family", ""),
            "host_symbols": row.get("host_symbols", ""),
            "sequence_status": "sequence_resolved_cleavage_test_candidate",
            "included_in_core_12": "yes" if included else "no",
            "exclusion_reason": "" if included else "sequence_resolved_but_lacks_independent_directional_expression_support",
            "nonproxy_direction_class": row.get("nonproxy_direction_class", ""),
            "direct_cleavage_evidence": row.get("direct_cleavage_evidence", "none"),
            "recommended_assay": "purified-protein or peptide cleavage assay followed by targeted binding/localization tests",
            "claim_boundary": "Sequence and compartment compatibility prioritize a test; they do not establish binding, cleavage, or hydraulic failure.",
        })
    return output


def build_revision_audits(root: Path, output_dir: Path) -> dict[str, Path]:
    builders = {
        "expression_dataset_summary": build_dataset_summary,
        "candidate_attrition": build_attrition,
        "priority_rule_matrix": build_priority_rules,
        "expression_member_audit": build_member_audit,
        "biochemical_test_candidates": build_biochemical_candidates,
    }
    outputs = {key: output_dir / name for key, name in OUTPUT_NAMES.items()}
    for key, builder in builders.items():
        rows = builder(root)
        write_tsv(outputs[key], rows)
        print(f"Wrote {len(rows)} rows to {outputs[key]}")
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--output-dir", type=Path, default=Path("results/tables"))
    args = parser.parse_args()
    build_revision_audits(args.root.resolve(), args.output_dir.resolve())


if __name__ == "__main__":
    main()
