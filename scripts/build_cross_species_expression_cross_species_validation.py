#!/usr/bin/env python3
"""Integrate targeted Tier B transcript mappings with five pine DE contrasts."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


DE_COLUMNS = ["gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
FIELDS = [
    "orthogroup_id", "mechanism_axes", "species", "contrast", "mapping_evidence",
    "can_support_tier_a", "mapped_transcript_count", "mapped_gene_count",
    "finite_padj_gene_count", "significant_gene_count", "cross_species_log2_fold_change",
    "cross_species_padj", "evidence_classification", "pmas_direction_status",
    "source_file", "claim_ceiling",
]
CLAIM = "Targeted homology and expression support are candidate evidence only; they do not establish exact orthology, direct effector targeting, causal resistance, coevolution, or a regional difference"
CONSENSUS_EVIDENCE = "multispecies_orthogroup_consensus_homology"


def _read(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _read_de(path: Path) -> dict[str, dict[str, str]]:
    rows = {}
    with Path(path).open(encoding="utf-8", newline="") as handle:
        for values in csv.reader(handle, delimiter="\t"):
            if values and values[0] == "gene_id":
                continue
            if len(values) != len(DE_COLUMNS):
                raise ValueError(f"Expected {len(DE_COLUMNS)} DESeq2 columns in {path}, found {len(values)}")
            rows[values[0]] = dict(zip(DE_COLUMNS, values))
    return rows


def _number(value: str) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _format(value: float | None) -> str:
    return "" if value is None else f"{value:.12g}"


def _gene_mappings(
    mapping_path: Path,
    tx2gene_path: Path,
) -> tuple[dict[str, set[str]], dict[str, int], dict[str, set[str]]]:
    tx2gene = {row["transcript_id"]: row["gene_id"] for row in _read(tx2gene_path)}
    gene_ogs: dict[str, set[str]] = defaultdict(set)
    transcript_counts: dict[str, int] = defaultdict(int)
    accepted = [
        row for row in _read(mapping_path)
        if row["mapping_status"] in {"accepted", "accepted_consensus"}
    ]
    for row in accepted:
        transcript = row["transcript_id"]
        if transcript not in tx2gene:
            raise ValueError(f"Accepted transcript missing from tx2gene: {transcript}")
        gene_ogs[tx2gene[transcript]].add(row["orthogroup_id"])
    result: dict[str, set[str]] = defaultdict(set)
    evidence_types: dict[str, set[str]] = defaultdict(set)
    for row in accepted:
        gene = tx2gene[row["transcript_id"]]
        if len(gene_ogs[gene]) == 1:
            orthogroup = next(iter(gene_ogs[gene]))
            result[orthogroup].add(gene)
            transcript_counts[orthogroup] += 1
            evidence_types[orthogroup].add(row.get("evidence_type", ""))
    return result, transcript_counts, evidence_types


def _summarize(genes: set[str], de: dict[str, dict[str, str]], pmas_direction: str) -> tuple[int, int, float | None, float | None, str]:
    tested = [de[gene] for gene in genes if gene in de]
    finite = [row for row in tested if _number(row["padj"]) is not None]
    significant = [row for row in finite if _number(row["padj"]) < 0.05]
    best_pool = finite or tested
    best = min(best_pool, key=lambda row: (_number(row["padj"]) if _number(row["padj"]) is not None else math.inf, -abs(_number(row["log2FoldChange"]) or 0))) if best_pool else None
    effect = _number(best["log2FoldChange"]) if best else None
    padj = _number(best["padj"]) if best else None
    if not significant:
        classification = "not_significant"
    else:
        signs = {1 if (_number(row["log2FoldChange"]) or 0) > 0 else -1 for row in significant if (_number(row["log2FoldChange"]) or 0) != 0}
        if len(signs) > 1:
            classification = "mixed_cross_species_response"
        elif pmas_direction not in {"positive", "negative"} or not signs:
            classification = "direction_not_comparable"
        else:
            expected = 1 if pmas_direction == "positive" else -1
            classification = "directionally_concordant" if signs == {expected} else "directionally_discordant"
    return len(finite), len(significant), effect, padj, classification


def build_validation(
    shortlist_path: Path,
    orthogroups_path: Path,
    pden_mapping_path: Path,
    pstrobus_mapping_path: Path,
    pden_tx2gene_path: Path,
    pstrobus_tx2gene_path: Path,
    pden_contrasts: dict[str, Path],
    pstrobus_contrasts: dict[str, Path],
    output_path: Path,
    stable_genes_path: Path | None = None,
    prig_xtae_contrasts: dict[str, Path] | None = None,
    pthun_contrasts: dict[str, Path] | None = None,
) -> list[dict[str, str]]:
    candidates = [row for row in _read(shortlist_path) if row["evidence_tier"] in {"Tier A", "Tier B"}]
    directions = {row["orthogroup_id"]: row["direction_status"] for row in _read(orthogroups_path)}
    pden_genes, pden_tx_counts, _pden_evidence = _gene_mappings(pden_mapping_path, pden_tx2gene_path)
    pstrobus_genes, pstrobus_tx_counts, pstrobus_evidence = _gene_mappings(pstrobus_mapping_path, pstrobus_tx2gene_path)
    datasets = [
        ("Pinus densiflora", "same_species_targeted_homology", "yes", pden_genes, pden_tx_counts, pden_contrasts, {}),
        ("Pinus strobus", "cross_species_homology_proxy", "no", pstrobus_genes, pstrobus_tx_counts, pstrobus_contrasts, pstrobus_evidence),
    ]
    if stable_genes_path is not None and prig_xtae_contrasts is not None:
        hybrid_genes: dict[str, set[str]] = defaultdict(set)
        for row in _read(stable_genes_path):
            if row["species_id"] == "ptae":
                hybrid_genes[row["orthogroup_id"]].add(row["gene_id"].split("|", 1)[-1])
        hybrid_counts = {orthogroup: len(genes) for orthogroup, genes in hybrid_genes.items()}
        datasets.append(("Pinus rigida x Pinus taeda", "P_taeda_parental_reference_proxy", "no", hybrid_genes, hybrid_counts, prig_xtae_contrasts, {}))
        if pthun_contrasts is not None:
            datasets.append(("Pinus thunbergii", "P_taeda_surrogate_reference_proxy", "no", hybrid_genes, hybrid_counts, pthun_contrasts, {}))
    elif pthun_contrasts is not None:
        raise ValueError("P. thunbergii contrasts require stable P. taeda orthogroup genes")
    all_contrasts = [pden_contrasts, pstrobus_contrasts] + ([prig_xtae_contrasts] if prig_xtae_contrasts else []) + ([pthun_contrasts] if pthun_contrasts else [])
    de_cache = {Path(path): _read_de(path) for contrasts in all_contrasts for path in contrasts.values()}
    rows: list[dict[str, str]] = []
    for candidate in sorted(candidates, key=lambda row: row["orthogroup_id"]):
        orthogroup = candidate["orthogroup_id"]
        for species, default_evidence, default_tier_a, genes_by_og, tx_counts, contrasts, evidence_by_og in datasets:
            genes = genes_by_og.get(orthogroup, set())
            consensus_only = evidence_by_og.get(orthogroup, set()) == {CONSENSUS_EVIDENCE}
            mapping_evidence = CONSENSUS_EVIDENCE if consensus_only else default_evidence
            tier_a = "yes" if consensus_only else default_tier_a
            for contrast, path in contrasts.items():
                if genes:
                    finite, significant, effect, padj, classification = _summarize(genes, de_cache[Path(path)], directions[orthogroup])
                    evidence = mapping_evidence
                else:
                    finite, significant, effect, padj, classification, evidence = 0, 0, None, None, "unavailable", "unavailable"
                rows.append({
                    "orthogroup_id": orthogroup,
                    "mechanism_axes": candidate.get("mechanism_axes", ""),
                    "species": species,
                    "contrast": contrast,
                    "mapping_evidence": evidence,
                    "can_support_tier_a": tier_a if genes else "no",
                    "mapped_transcript_count": str(tx_counts.get(orthogroup, 0)),
                    "mapped_gene_count": str(len(genes)),
                    "finite_padj_gene_count": str(finite),
                    "significant_gene_count": str(significant),
                    "cross_species_log2_fold_change": _format(effect),
                    "cross_species_padj": _format(padj),
                    "evidence_classification": classification,
                    "pmas_direction_status": directions[orthogroup],
                    "source_file": Path(path).name,
                    "claim_ceiling": CLAIM,
                })
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with Path(output_path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["shortlist", "orthogroups", "pden-mapping", "pstrobus-mapping", "pden-tx2gene", "pstrobus-tx2gene", "pden-primary", "pden-bxyl-water", "pden-bthai-water", "pstrobus-two-week", "pstrobus-four-week", "stable-genes", "prig-xtae-two-week", "prig-xtae-four-week", "output"]:
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--pthun-two-week", type=Path)
    parser.add_argument("--pthun-four-week", type=Path)
    args = parser.parse_args()
    if bool(args.pthun_two_week) != bool(args.pthun_four_week):
        parser.error("--pthun-two-week and --pthun-four-week must be supplied together")
    rows = build_validation(
        args.shortlist, args.orthogroups, args.pden_mapping, args.pstrobus_mapping,
        args.pden_tx2gene, args.pstrobus_tx2gene,
        {"pathogen_associated": args.pden_primary, "bxyl_vs_water": args.pden_bxyl_water, "bthai_vs_water": args.pden_bthai_water},
        {"two_week_vs_zero": args.pstrobus_two_week, "four_week_vs_zero": args.pstrobus_four_week}, args.output,
        stable_genes_path=args.stable_genes,
        prig_xtae_contrasts={"two_week_vs_zero": args.prig_xtae_two_week, "four_week_vs_zero": args.prig_xtae_four_week},
        pthun_contrasts=(
            {"two_week_vs_zero": args.pthun_two_week, "four_week_vs_zero": args.pthun_four_week}
            if args.pthun_two_week else None
        ),
    )
    print(f"Wrote {len(rows)} Tier B validation rows")


if __name__ == "__main__":
    main()
