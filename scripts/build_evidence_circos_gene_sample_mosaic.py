#!/usr/bin/env python3
"""Freeze exact gene-by-sample and gene-by-contrast data for the evidence Circos."""

from __future__ import annotations

import argparse
import csv
import logging
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


LOG = logging.getLogger(__name__)
DE_COLUMNS = ["gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
OUTPUT_COLUMNS = [
    "record_type", "dataset_id", "species", "sector_id", "orthogroup_id", "og_order",
    "gene_id", "gene_order", "gene_count", "sample_id", "sample_order", "sample_count",
    "condition", "replicate", "expression_value", "expression_scale", "log2_expression",
    "expression_zscore", "contrast_id", "contrast_order", "contrast_count", "log2fc",
    "pvalue", "fdr", "significant", "mapping_evidence", "mapping_tier", "is_proxy",
    "is_representative", "source_file", "claim_boundary",
]


DATASETS = {
    "pmas": {
        "species": "Pinus massoniana",
        "mapping_evidence": "exact candidate orthogroup member",
        "mapping_tier": "direct_orthogroup_member",
        "is_proxy": "no",
        "claim": "Controlled expression evidence for candidate orthogroup members; interaction causality and resistance causality are not established.",
    },
    "pden": {
        "species": "Pinus densiflora",
        "mapping_evidence": "same-species targeted homology",
        "mapping_tier": "same_species_targeted_homology",
        "is_proxy": "no",
        "claim": "Same-species targeted homology and expression support; exact orthology and resistance causality are not established.",
    },
    "pstrobus": {
        "species": "Pinus strobus",
        "mapping_evidence": "multispecies orthogroup consensus homology",
        "mapping_tier": "multispecies_consensus_homology",
        "is_proxy": "yes",
        "claim": "Cross-species consensus-homology expression proxy; exact orthology and resistance causality are not established.",
    },
    "pthun": {
        "species": "Pinus thunbergii",
        "mapping_evidence": "P. taeda surrogate-reference gene",
        "mapping_tier": "P_taeda_surrogate_reference_proxy",
        "is_proxy": "yes",
        "claim": "P. taeda surrogate-reference expression proxy for P. thunbergii; native gene identity and resistance causality are not established.",
    },
    "prig_xtae": {
        "species": "Pinus rigida x Pinus taeda",
        "mapping_evidence": "P. taeda parental-reference gene",
        "mapping_tier": "P_taeda_parental_reference_proxy",
        "is_proxy": "yes",
        "claim": "P. taeda parental-reference expression proxy for the hybrid; parental origin and resistance causality are not established.",
    },
}


def read_tsv(path: Path, **kwargs: object) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, **kwargs)


def finite_number(value: object) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return math.nan
    return number if math.isfinite(number) else math.nan


def read_de(path: Path) -> pd.DataFrame:
    first = path.open(encoding="utf-8").readline().rstrip("\n\r").split("\t")
    header = 0 if first and first[0] == "gene_id" else None
    frame = pd.read_csv(path, sep="\t", dtype=str, header=header, keep_default_na=False)
    if header is None:
        if frame.shape[1] != len(DE_COLUMNS):
            raise ValueError(f"Expected seven DESeq2 columns in {path}, found {frame.shape[1]}")
        frame.columns = DE_COLUMNS
    return frame.set_index("gene_id", drop=False)


def unique_gene_mappings(mapping_path: Path, tx2gene_path: Path) -> dict[str, set[str]]:
    mapping = read_tsv(mapping_path)
    mapping = mapping[mapping["mapping_status"].isin(["accepted", "accepted_consensus"])].copy()
    tx2gene = read_tsv(tx2gene_path).set_index("transcript_id")["gene_id"].to_dict()
    mapping["gene_id"] = mapping["transcript_id"].map(tx2gene)
    if mapping["gene_id"].eq("").any() or mapping["gene_id"].isna().any():
        missing = mapping.loc[mapping["gene_id"].isna() | mapping["gene_id"].eq(""), "transcript_id"].iloc[0]
        raise ValueError(f"Accepted transcript missing from tx2gene: {missing}")
    gene_og_count = mapping.groupby("gene_id")["orthogroup_id"].nunique()
    mapping = mapping[mapping["gene_id"].map(gene_og_count).eq(1)]
    return {
        orthogroup: set(group["gene_id"])
        for orthogroup, group in mapping.groupby("orthogroup_id", sort=False)
    }


def candidate_gene_sets(root: Path, candidate_ids: set[str]) -> dict[str, dict[str, set[str]]]:
    audit = read_tsv(root / "results/tables/manuscript_expression_member_audit.tsv")
    pmas = audit[
        audit["species"].eq("Pinus massoniana")
        & audit["evidence_unit"].eq("orthogroup_member")
        & audit["orthogroup_id"].isin(candidate_ids)
        & audit["feature_id"].ne("")
    ]
    pmas_genes = {
        orthogroup: set(group["feature_id"])
        for orthogroup, group in pmas.groupby("orthogroup_id", sort=False)
    }
    pden_genes = unique_gene_mappings(
        root / "results/tables/cross_species_expression_pden_transcript_mapping.tsv",
        root / "data/interim/pden_published_trinity_tx2gene.tsv",
    )
    pstrobus_genes = unique_gene_mappings(
        root / "results/tables/pstrobus_orthogroup_consensus_mapping.tsv",
        root / "results/tables/effector_target_network_pstrobus_GIIE_tx2gene.tsv",
    )
    stable = read_tsv(root / "results/tables/stable_orthogroup_genes.tsv")
    stable = stable[stable["species_id"].eq("ptae") & stable["orthogroup_id"].isin(candidate_ids)].copy()
    stable["gene_id"] = stable["gene_id"].str.split("|", n=1).str[-1]
    ptae_genes = {
        orthogroup: set(group["gene_id"])
        for orthogroup, group in stable.groupby("orthogroup_id", sort=False)
    }
    return {
        "pmas": {og: genes for og, genes in pmas_genes.items() if og in candidate_ids},
        "pden": {og: genes for og, genes in pden_genes.items() if og in candidate_ids},
        "pstrobus": {og: genes for og, genes in pstrobus_genes.items() if og in candidate_ids},
        "pthun": ptae_genes,
        "prig_xtae": {og: set(genes) for og, genes in ptae_genes.items()},
    }


def de_sources(root: Path) -> dict[str, list[tuple[str, str, Path]]]:
    return {
        "pden": [
            ("pathogen_associated", "pden_expr", root / "data/external/galaxy/effector_target_network_pden_full_deseq2/effector_target_network_bxyl_vs_bthai_result.tabular"),
            ("bxyl_vs_water", "pden_expr", root / "data/external/galaxy/effector_target_network_pden_full_deseq2/effector_target_network_bxyl_vs_water_result.tabular"),
            ("bthai_vs_water", "pden_expr", root / "data/external/galaxy/effector_target_network_pden_full_deseq2/effector_target_network_bthai_vs_water_result.tabular"),
        ],
        "pstrobus": [
            ("two_week_vs_zero", "pstrobus_expr", root / "results/tables/effector_target_network_pstrobus_GIIE_DESeq2_2w_vs_0w.tsv"),
            ("four_week_vs_zero", "pstrobus_expr", root / "results/tables/effector_target_network_pstrobus_GIIE_DESeq2_4w_vs_0w.tsv"),
        ],
        "pthun": [
            ("two_week_vs_zero", "pthun_expr", root / "results/tables/apoplastic_cell_wall_pthun_DESeq2_2w_vs_0w.tsv"),
            ("four_week_vs_zero", "pthun_expr", root / "results/tables/apoplastic_cell_wall_pthun_DESeq2_4w_vs_0w.tsv"),
        ],
        "prig_xtae": [
            ("two_week_vs_zero", "prigtae_expr", root / "results/tables/apoplastic_cell_wall_prig_xtae_DESeq2_2w_vs_0w.tsv"),
            ("four_week_vs_zero", "prigtae_expr", root / "results/tables/apoplastic_cell_wall_prig_xtae_DESeq2_4w_vs_0w.tsv"),
        ],
    }


def build_de_records(
    root: Path,
    genes_by_dataset: dict[str, dict[str, set[str]]],
    og_order: dict[str, int],
) -> tuple[list[dict[str, object]], dict[str, dict[str, list[str]]]]:
    records: list[dict[str, object]] = []
    ordering_metrics: dict[str, dict[str, dict[str, object]]] = defaultdict(lambda: defaultdict(lambda: {"significant": False, "max_abs_fc": -1.0}))
    sources = de_sources(root)

    audit = read_tsv(root / "results/tables/manuscript_expression_member_audit.tsv")
    audit = audit[
        audit["species"].eq("Pinus massoniana")
        & audit["evidence_unit"].eq("orthogroup_member")
        & audit["feature_id"].ne("")
    ].copy()
    pmas_contrasts = [
        ("resistant_pwn_vs_water", "pmas_r"),
        ("susceptible_pwn_vs_water", "pmas_s"),
        ("genotype_by_inoculum", "pmas_gxi"),
    ]
    sources["pmas"] = [(contrast, sector, Path("manuscript_expression_member_audit.tsv")) for contrast, sector in pmas_contrasts]

    for dataset_id, contrasts in sources.items():
        metadata = DATASETS[dataset_id]
        contrast_count = len(contrasts)
        for contrast_order, (contrast_id, sector_id, path) in enumerate(contrasts, start=1):
            if dataset_id == "pmas":
                subset = audit[audit["contrast"].eq(contrast_id)].set_index(["orthogroup_id", "feature_id"])
                de = None
            else:
                subset = None
                de = read_de(path)
            contrast_rows: list[dict[str, object]] = []
            for orthogroup_id in sorted(genes_by_dataset[dataset_id], key=lambda og: og_order.get(og, 10**9)):
                for gene_id in sorted(genes_by_dataset[dataset_id][orthogroup_id]):
                    if dataset_id == "pmas" and (orthogroup_id, gene_id) in subset.index:
                        hit = subset.loc[(orthogroup_id, gene_id)]
                        if isinstance(hit, pd.DataFrame):
                            hit = hit.iloc[0]
                        log2fc = finite_number(hit["log2_fold_change"])
                        pvalue = finite_number(hit["raw_pvalue"])
                        fdr = finite_number(hit["whole_transcriptome_padj"])
                    elif dataset_id != "pmas" and gene_id in de.index:
                        hit = de.loc[gene_id]
                        if isinstance(hit, pd.DataFrame):
                            hit = hit.iloc[0]
                        log2fc = finite_number(hit["log2FoldChange"])
                        pvalue = finite_number(hit["pvalue"])
                        fdr = finite_number(hit["padj"])
                    else:
                        log2fc = pvalue = fdr = math.nan
                    significant = "yes" if math.isfinite(fdr) and fdr < 0.05 else ("no" if math.isfinite(fdr) else "not_tested")
                    metric = ordering_metrics[dataset_id][gene_id]
                    metric["significant"] = bool(metric["significant"] or significant == "yes")
                    if math.isfinite(log2fc):
                        metric["max_abs_fc"] = max(float(metric["max_abs_fc"]), abs(log2fc))
                    contrast_rows.append({
                        "record_type": "contrast_de", "dataset_id": dataset_id,
                        "species": metadata["species"], "sector_id": sector_id,
                        "orthogroup_id": orthogroup_id, "og_order": og_order[orthogroup_id],
                        "gene_id": gene_id, "sample_id": "", "sample_order": "", "sample_count": "",
                        "condition": "", "replicate": "", "expression_value": "", "expression_scale": "",
                        "log2_expression": "", "expression_zscore": "", "contrast_id": contrast_id,
                        "contrast_order": contrast_order, "contrast_count": contrast_count,
                        "log2fc": log2fc, "pvalue": pvalue, "fdr": fdr, "significant": significant,
                        "mapping_evidence": metadata["mapping_evidence"], "mapping_tier": metadata["mapping_tier"],
                        "is_proxy": metadata["is_proxy"], "is_representative": "no",
                        "source_file": path.name, "claim_boundary": metadata["claim"],
                    })
            by_og: dict[str, list[dict[str, object]]] = defaultdict(list)
            for row in contrast_rows:
                by_og[str(row["orthogroup_id"])].append(row)
            for group in by_og.values():
                tested = [row for row in group if math.isfinite(float(row["fdr"]))]
                pool = tested or [row for row in group if math.isfinite(float(row["log2fc"]))]
                if pool:
                    representative = min(
                        pool,
                        key=lambda row: (
                            float(row["fdr"]) if math.isfinite(float(row["fdr"])) else math.inf,
                            -abs(float(row["log2fc"])) if math.isfinite(float(row["log2fc"])) else math.inf,
                            str(row["gene_id"]),
                        ),
                    )
                    representative["is_representative"] = "yes"
            records.extend(contrast_rows)

    gene_orders: dict[str, dict[str, list[str]]] = defaultdict(dict)
    for dataset_id, by_og in genes_by_dataset.items():
        for orthogroup_id, genes in by_og.items():
            gene_orders[dataset_id][orthogroup_id] = sorted(
                genes,
                key=lambda gene: (
                    not bool(ordering_metrics[dataset_id][gene]["significant"]),
                    -float(ordering_metrics[dataset_id][gene]["max_abs_fc"]),
                    gene,
                ),
            )
    for row in records:
        order = gene_orders[str(row["dataset_id"])][str(row["orthogroup_id"])]
        row["gene_order"] = order.index(str(row["gene_id"])) + 1
        row["gene_count"] = len(order)
    return records, gene_orders


def read_quant_gene_tpm(path: Path, tx2gene: dict[str, str] | None, selected: set[str]) -> dict[str, float]:
    frame = pd.read_csv(path, sep="\t", usecols=["Name", "TPM"])
    if tx2gene is None:
        frame["gene_id"] = frame["Name"]
    else:
        frame["gene_id"] = frame["Name"].map(tx2gene)
    frame = frame[frame["gene_id"].isin(selected)]
    return frame.groupby("gene_id", sort=False)["TPM"].sum().to_dict()


def read_pden_counts(path: Path) -> tuple[list[str], dict[str, dict[str, float]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        sample_ids = next(reader)
        values: dict[str, dict[str, float]] = {}
        for row in reader:
            if len(row) != len(sample_ids) + 1:
                raise ValueError(f"Malformed P. densiflora count row for {row[0] if row else 'unknown'}")
            values[row[0]] = {sample: finite_number(value) for sample, value in zip(sample_ids, row[1:])}
    return sample_ids, values


def sample_sheets(root: Path) -> dict[str, pd.DataFrame]:
    sheets = {
        "pmas": read_tsv(root / "config/effector_target_network_pmas_1dpi_primary_sample_sheet.tsv"),
        "pden": read_tsv(root / "config/effector_target_network_pden_deseq2_sample_sheet.tsv"),
        "pstrobus": read_tsv(root / "config/effector_target_network_pstrobus_deseq2_sample_sheet.tsv"),
        "pthun": read_tsv(root / "config/apoplastic_cell_wall_pthun_deseq2_sample_sheet.tsv"),
        "prig_xtae": read_tsv(root / "config/apoplastic_cell_wall_prig_xtae_deseq2_sample_sheet.tsv"),
    }
    for dataset_id, sheet in sheets.items():
        sheet["condition_plot"] = sheet["condition"] if "condition" in sheet else sheet["timepoint"]
        sheets[dataset_id] = sheet
    return sheets


def expression_matrices(
    root: Path,
    genes_by_dataset: dict[str, dict[str, set[str]]],
    sheets: dict[str, pd.DataFrame],
) -> tuple[dict[str, pd.DataFrame], dict[str, str], dict[str, dict[str, str]]]:
    matrices: dict[str, pd.DataFrame] = {}
    scales = {"pmas": "TPM", "pden": "tximport_count", "pstrobus": "TPM", "pthun": "TPM", "prig_xtae": "TPM"}
    source_by_sample: dict[str, dict[str, str]] = defaultdict(dict)
    selected = {dataset: set().union(*by_og.values()) for dataset, by_og in genes_by_dataset.items()}

    for dataset_id, tx2gene_path, quant_root in [
        ("pmas", root / "results/tables/effector_target_network_pmas_GigaDB102688_tx2gene.tsv", root / "data/external/galaxy/effector_target_network_pmas_1dpi_quant"),
        ("pthun", None, root / "data/external/galaxy/apoplastic_cell_wall_pthun_salmon"),
        ("prig_xtae", None, root / "data/external/galaxy/apoplastic_cell_wall_prig_xtae_salmon"),
    ]:
        tx2gene = None if tx2gene_path is None else read_tsv(tx2gene_path).set_index("transcript_id")["gene_id"].to_dict()
        sample_values: dict[str, dict[str, float]] = {}
        for _, sample in sheets[dataset_id].iterrows():
            folder = sample["sample_id"] if dataset_id == "pmas" else sample["run_accession"]
            path = quant_root / folder / "quant.sf"
            sample_values[sample["sample_id"]] = read_quant_gene_tpm(path, tx2gene, selected[dataset_id])
            source_by_sample[dataset_id][sample["sample_id"]] = str(path.relative_to(root))
        matrices[dataset_id] = pd.DataFrame(sample_values).reindex(sorted(selected[dataset_id])).fillna(0.0)

    pden_path = root / "data/external/galaxy/effector_target_network_pden_full_salmon/effector_target_network_pden_tximport_gene_level_values.tabular"
    pden_samples, pden_values = read_pden_counts(pden_path)
    matrices["pden"] = pd.DataFrame.from_dict(pden_values, orient="index").reindex(sorted(selected["pden"])).fillna(0.0)
    matrices["pden"] = matrices["pden"].reindex(columns=pden_samples)
    source_by_sample["pden"] = {sample: str(pden_path.relative_to(root)) for sample in pden_samples}

    manifest = read_tsv(root / "results/tables/effector_target_network_pstrobus_quant_manifest.tsv")
    manifest_paths = {row["sample_id"]: root / Path(row["local_path"]) for _, row in manifest.iterrows()}
    pstrobus_tx2gene = read_tsv(root / "results/tables/effector_target_network_pstrobus_GIIE_tx2gene.tsv").set_index("transcript_id")["gene_id"].to_dict()
    pstrobus_values: dict[str, dict[str, float]] = {}
    for _, sample in sheets["pstrobus"].iterrows():
        path = manifest_paths[sample["sample_id"]]
        pstrobus_values[sample["sample_id"]] = read_quant_gene_tpm(path, pstrobus_tx2gene, selected["pstrobus"])
        source_by_sample["pstrobus"][sample["sample_id"]] = str(path.relative_to(root))
    matrices["pstrobus"] = pd.DataFrame(pstrobus_values).reindex(sorted(selected["pstrobus"])).fillna(0.0)
    return matrices, scales, source_by_sample


def build_sample_records(
    matrices: dict[str, pd.DataFrame],
    scales: dict[str, str],
    source_by_sample: dict[str, dict[str, str]],
    sheets: dict[str, pd.DataFrame],
    genes_by_dataset: dict[str, dict[str, set[str]]],
    gene_orders: dict[str, dict[str, list[str]]],
    og_order: dict[str, int],
) -> list[dict[str, object]]:
    sector_filters = {
        "pmas": [
            ("pmas_r", lambda sheet: sheet["phenotype"].eq("resistant")),
            ("pmas_s", lambda sheet: sheet["phenotype"].eq("susceptible")),
            ("pmas_gxi", lambda sheet: pd.Series(True, index=sheet.index)),
        ],
        "pden": [("pden_expr", lambda sheet: pd.Series(True, index=sheet.index))],
        "pstrobus": [("pstrobus_expr", lambda sheet: pd.Series(True, index=sheet.index))],
        "pthun": [("pthun_expr", lambda sheet: pd.Series(True, index=sheet.index))],
        "prig_xtae": [("prigtae_expr", lambda sheet: pd.Series(True, index=sheet.index))],
    }
    rows: list[dict[str, object]] = []
    for dataset_id, matrix in matrices.items():
        metadata = DATASETS[dataset_id]
        logged = np.log2(matrix.astype(float) + 1.0)
        means = logged.mean(axis=1)
        stds = logged.std(axis=1, ddof=0)
        zscores = logged.sub(means, axis=0).div(stds.replace(0, np.nan), axis=0).fillna(0.0)
        for sector_id, choose in sector_filters[dataset_id]:
            sector_samples = sheets[dataset_id][choose(sheets[dataset_id])].copy()
            sample_count = len(sector_samples)
            for sample_order, (_, sample) in enumerate(sector_samples.iterrows(), start=1):
                sample_id = sample["sample_id"]
                for orthogroup_id in sorted(genes_by_dataset[dataset_id], key=lambda og: og_order.get(og, 10**9)):
                    order = gene_orders[dataset_id][orthogroup_id]
                    for gene_order, gene_id in enumerate(order, start=1):
                        rows.append({
                            "record_type": "sample_expression", "dataset_id": dataset_id,
                            "species": metadata["species"], "sector_id": sector_id,
                            "orthogroup_id": orthogroup_id, "og_order": og_order[orthogroup_id],
                            "gene_id": gene_id, "gene_order": gene_order, "gene_count": len(order),
                            "sample_id": sample_id, "sample_order": sample_order, "sample_count": sample_count,
                            "condition": sample["condition_plot"], "replicate": sample["replicate"],
                            "expression_value": float(matrix.at[gene_id, sample_id]), "expression_scale": scales[dataset_id],
                            "log2_expression": float(logged.at[gene_id, sample_id]),
                            "expression_zscore": float(zscores.at[gene_id, sample_id]),
                            "contrast_id": "", "contrast_order": "", "contrast_count": "",
                            "log2fc": "", "pvalue": "", "fdr": "", "significant": "",
                            "mapping_evidence": metadata["mapping_evidence"], "mapping_tier": metadata["mapping_tier"],
                            "is_proxy": metadata["is_proxy"], "is_representative": "no",
                            "source_file": source_by_sample[dataset_id][sample_id],
                            "claim_boundary": metadata["claim"],
                        })
    return rows


def build_mosaic(root: Path, output: Path) -> pd.DataFrame:
    candidates = read_tsv(root / "results/tables/manuscript_multilayer_circos_candidates.tsv")
    candidates["display_order"] = pd.to_numeric(candidates["display_order"], errors="raise").astype(int)
    candidates = candidates.sort_values("display_order")
    og_order = dict(zip(candidates["orthogroup_id"], candidates["display_order"]))
    genes_by_dataset = candidate_gene_sets(root, set(og_order))
    for dataset_id, by_og in genes_by_dataset.items():
        LOG.info("%s: %d genes across %d candidate orthogroups", dataset_id, sum(map(len, by_og.values())), len(by_og))
    de_records, gene_orders = build_de_records(root, genes_by_dataset, og_order)
    sheets = sample_sheets(root)
    matrices, scales, source_by_sample = expression_matrices(root, genes_by_dataset, sheets)
    sample_records = build_sample_records(
        matrices, scales, source_by_sample, sheets, genes_by_dataset, gene_orders, og_order
    )
    result = pd.DataFrame(sample_records + de_records, columns=OUTPUT_COLUMNS)
    result = result.sort_values(
        ["sector_id", "record_type", "og_order", "gene_order", "sample_order", "contrast_order"],
        kind="stable",
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, sep="\t", index=False, na_rep="NA")
    LOG.info("Wrote %d gene-sample/contrast records to %s", len(result), output)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    build_mosaic(args.project_root.resolve(), args.output)


if __name__ == "__main__":
    main()
