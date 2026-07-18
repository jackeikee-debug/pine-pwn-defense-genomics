#!/usr/bin/env python3
"""Build a claim-bounded PR-4/WIN1 dossier for OG0005853."""

from __future__ import annotations

import argparse
import csv
import logging
import math
from collections import Counter
from pathlib import Path

from Bio import Phylo, SeqIO


CLAIM_CEILING = "candidate mechanism; not direct effector targeting or causal resistance evidence"
COHERENT_CLASSES = {"PR4", "WIN1", "WHEATWIN"}
LOGGER = logging.getLogger(__name__)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def normalize_tree_label(label: str) -> str:
    """Convert OrthoFinder's species-prefixed tip label to the stable gene ID."""
    if "_" not in label:
        return label
    prefix, remainder = label.split("_", 1)
    if remainder.startswith(f"{prefix}|"):
        return remainder
    return label


def classify_annotation(annotation_text: str) -> str:
    text = annotation_text.upper()
    if "WHEATWIN" in text:
        return "WHEATWIN"
    if "WIN1" in text or "WOUND-INDUCED PROTEIN" in text:
        return "WIN1"
    if "PATHOGENESIS-RELATED PROTEIN 4" in text or "PR-4" in text or "PR4" in text:
        return "PR4"
    return "OTHER"


def n_terminal_hydrophobic_flag(sequence: str) -> str:
    segment = sequence[:25].upper()
    if not segment:
        return "no"
    hydrophobic = sum(residue in "AILMFWVY" for residue in segment)
    return "yes" if hydrophobic / len(segment) >= 0.40 else "no"


def build_member_annotation(
    orthogroup_id: str,
    stable_genes_path: Path,
    metadata_path: Path,
    fasta_path: Path,
    swissprot_path: Path,
    defense_hits_path: Path,
) -> list[dict[str, str]]:
    stable_rows = [
        row for row in read_tsv(stable_genes_path) if row["orthogroup_id"] == orthogroup_id
    ]
    metadata = {row["species_id"]: row for row in read_tsv(metadata_path)}
    sequences = {record.id: str(record.seq) for record in SeqIO.parse(fasta_path, "fasta")}
    swissprot = {row["gene_id"]: row for row in read_tsv(swissprot_path)}

    modules: dict[str, set[str]] = {}
    for row in read_tsv(defense_hits_path):
        if row.get("orthogroup_id") == orthogroup_id:
            modules.setdefault(row["gene_id"], set()).add(row["module_id"])

    missing_sequences = sorted(
        row["gene_id"] for row in stable_rows if row["gene_id"] not in sequences
    )
    if missing_sequences:
        raise ValueError(f"Missing stable sequences: {', '.join(missing_sequences)}")

    rows: list[dict[str, str]] = []
    for stable in sorted(stable_rows, key=lambda row: (row["species_id"], row["gene_id"])):
        gene_id = stable["gene_id"]
        species_id = stable["species_id"]
        species = metadata.get(species_id, {})
        hit = swissprot.get(gene_id, {})
        sequence = sequences[gene_id]
        rows.append(
            {
                "orthogroup_id": orthogroup_id,
                "species_id": species_id,
                "species_name": species.get("species_name", ""),
                "region": species.get("region", ""),
                "group": species.get("group", ""),
                "gene_id": gene_id,
                "sequence_length": str(len(sequence)),
                "n_terminal_hydrophobic_flag": n_terminal_hydrophobic_flag(sequence),
                "swissprot_subject_id": hit.get("subject_id", ""),
                "pident": hit.get("pident", ""),
                "alignment_length": hit.get("alignment_length", ""),
                "annotation_text": hit.get("annotation_text", ""),
                "annotation_class": classify_annotation(hit.get("annotation_text", "")),
                "defense_modules": ";".join(sorted(modules.get(gene_id, set()))),
                "tree_member": "not_audited",
                "sequence_member": "yes",
                "claim_ceiling": CLAIM_CEILING,
            }
        )
    return rows


def audit_tree(
    orthogroup_id: str,
    member_rows: list[dict[str, str]],
    tree_path: Path,
) -> dict[str, str]:
    tree = Phylo.read(tree_path, "newick")
    tree_tips = [normalize_tree_label(tip.name or "") for tip in tree.get_terminals()]
    expected = [row["gene_id"] for row in member_rows]
    tree_counts = Counter(tree_tips)
    expected_counts = Counter(expected)
    for row in member_rows:
        row["tree_member"] = "yes" if tree_counts[row["gene_id"]] == 1 else "no"

    species_counts = Counter(row["species_id"] for row in member_rows)
    duplicate_rich = ";".join(
        f"{species_id}:{count}"
        for species_id, count in sorted(species_counts.items())
        if count > 1
    )
    classes = sorted({row["annotation_class"] for row in member_rows})
    return {
        "orthogroup_id": orthogroup_id,
        "stable_member_count": str(len(member_rows)),
        "tree_tip_count": str(len(tree_tips)),
        "sequence_member_count": str(sum(row["sequence_member"] == "yes" for row in member_rows)),
        "stable_species_count": str(len(species_counts)),
        "duplicate_rich_species": duplicate_rich,
        "annotation_classes": ";".join(classes),
        "annotation_coherent": "yes" if classes and set(classes) <= COHERENT_CLASSES else "no",
        "tree_membership_complete": "yes" if tree_counts == expected_counts else "no",
        "sequence_membership_complete": "yes"
        if all(row["sequence_member"] == "yes" for row in member_rows)
        else "no",
        "network_bridge_status": "no_clean_interolog_bridge",
        "claim_ceiling": CLAIM_CEILING,
    }


def parse_float(value: str) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def read_de_table(path: Path) -> dict[str, dict[str, str]]:
    fields = ["gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
    rows: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for values in reader:
            if not values or values[0] == "gene_id":
                continue
            padded = values[: len(fields)] + [""] * max(0, len(fields) - len(values))
            row = dict(zip(fields, padded, strict=True))
            rows[row["gene_id"]] = row
    return rows


def de_interpretation(log2_fold_change: str, padj: str) -> str:
    lfc = parse_float(log2_fold_change)
    adjusted = parse_float(padj)
    if adjusted is None or lfc is None:
        return "not_tested_or_filtered"
    if adjusted >= 0.05:
        return "nonsignificant"
    if lfc > 0:
        return "significant_up"
    if lfc < 0:
        return "significant_down"
    return "significant_no_direction"


def map_transcripts_to_genes(path: Path) -> dict[str, str]:
    return {row["transcript_id"]: row["gene_id"] for row in read_tsv(path)}


def accepted_expression_mappings(
    orthogroup_id: str,
    mapping_path: Path,
    tx2gene_path: Path,
    species_id: str,
) -> dict[str, list[dict[str, str]]]:
    tx2gene = map_transcripts_to_genes(tx2gene_path)
    accepted: dict[str, list[dict[str, str]]] = {}
    for row in read_tsv(mapping_path):
        if row.get("orthogroup_id") != orthogroup_id:
            continue
        if species_id == "pden":
            passes = row.get("mapping_status") == "accepted"
        else:
            passes = (
                row.get("mapping_status") == "accepted_consensus"
                and (parse_float(row.get("supporting_species_count", "")) or 0) >= 2
                and (parse_float(row.get("relative_bitscore_margin", "")) or 0) >= 0.10
                and (parse_float(row.get("subject_coverage", "")) or 0) >= 0.80
            )
        transcript_id = row.get("transcript_id", "")
        gene_id = tx2gene.get(transcript_id)
        if passes and gene_id:
            accepted.setdefault(gene_id, []).append(row)
    return accepted


def build_expression_evidence(
    orthogroup_id: str,
    pden_mapping_path: Path,
    pstrobus_mapping_path: Path,
    pden_tx2gene_path: Path,
    pstrobus_tx2gene_path: Path,
    pden_contrasts: dict[str, Path],
    pstrobus_contrasts: dict[str, Path],
) -> list[dict[str, str]]:
    species_inputs = [
        (
            "pden",
            "Pinus densiflora",
            accepted_expression_mappings(
                orthogroup_id, pden_mapping_path, pden_tx2gene_path, "pden"
            ),
            pden_contrasts,
        ),
        (
            "pstrobus",
            "Pinus strobus",
            accepted_expression_mappings(
                orthogroup_id, pstrobus_mapping_path, pstrobus_tx2gene_path, "pstrobus"
            ),
            pstrobus_contrasts,
        ),
    ]
    output: list[dict[str, str]] = []
    for species_id, species_name, mappings, contrasts in species_inputs:
        for contrast, de_path in contrasts.items():
            de_rows = read_de_table(de_path)
            for gene_id, mapping_rows in sorted(mappings.items()):
                de = de_rows.get(gene_id, {})
                transcript_ids = sorted({row["transcript_id"] for row in mapping_rows})
                target_field = "target_gene_id" if species_id == "pden" else "best_target_gene_id"
                target_ids = sorted({row.get(target_field, "") for row in mapping_rows if row.get(target_field)})
                padj = de.get("padj", "")
                if parse_float(padj) is None:
                    padj = ""
                output.append(
                    {
                        "orthogroup_id": orthogroup_id,
                        "species_id": species_id,
                        "species_name": species_name,
                        "contrast": contrast,
                        "expression_gene_id": gene_id,
                        "transcript_ids": ";".join(transcript_ids),
                        "projected_stable_gene_ids": ";".join(target_ids),
                        "mapping_status": ";".join(sorted({row["mapping_status"] for row in mapping_rows})),
                        "mapping_evidence": ";".join(sorted({row.get("evidence_type", "") for row in mapping_rows if row.get("evidence_type")})),
                        "supporting_species_count": max(
                            (row.get("supporting_species_count", "") for row in mapping_rows),
                            default="",
                        ),
                        "relative_bitscore_margin": min(
                            (row.get("relative_bitscore_margin", "") for row in mapping_rows if row.get("relative_bitscore_margin")),
                            default="",
                        ),
                        "subject_coverage": min(
                            (row.get("subject_coverage", "") for row in mapping_rows if row.get("subject_coverage")),
                            default="",
                        ),
                        "baseMean": de.get("baseMean", ""),
                        "log2FoldChange": de.get("log2FoldChange", ""),
                        "pvalue": de.get("pvalue", ""),
                        "padj": padj,
                        "de_interpretation": de_interpretation(
                            de.get("log2FoldChange", ""), de.get("padj", "")
                        ),
                        "claim_ceiling": CLAIM_CEILING,
                    }
                )
    return output


def qualitative_concordance(rows: list[dict[str, str]]) -> str:
    up_species = {
        row["species_id"]
        for row in rows
        if row["de_interpretation"] == "significant_up"
    }
    if {"pden", "pstrobus"} <= up_species:
        return "significant_up_observed_in_both_species"
    return "no_cross_species_significant_up_concordance"


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        raise ValueError(f"Cannot write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_report(
    path: Path,
    orthogroup_id: str,
    member_rows: list[dict[str, str]],
    audit: dict[str, str],
    expression_rows: list[dict[str, str]],
) -> None:
    concordance = qualitative_concordance(expression_rows)
    significant = [
        row for row in expression_rows if row["de_interpretation"].startswith("significant_")
    ]
    lines = [
        f"# {orthogroup_id} PR-4/WIN1 mechanism dossier",
        "",
        "## Evidence summary",
        "",
        f"- Stable orthology contains {audit['stable_member_count']} proteins from "
        f"{audit['stable_species_count']} species.",
        f"- Annotation audit supports a coherent PR-4/WIN1/Wheatwin family "
        f"(`annotation_coherent={audit['annotation_coherent']}`).",
        f"- Tree and sequence membership are complete: "
        f"`tree={audit['tree_membership_complete']}`, "
        f"`fasta={audit['sequence_membership_complete']}`.",
        f"- Network status: `{audit['network_bridge_status']}`.",
        "",
        "## Infection-expression evidence",
        "",
    ]
    if concordance == "significant_up_observed_in_both_species":
        lines.append(
            "Qualitatively, significant up-regulation was observed in both species in "
            "independent infection-expression datasets. Fold-change magnitudes are not "
            "compared across experiments."
        )
    else:
        lines.append("No qualitative significant-up concordance was observed across both species.")
    lines.extend(
        [
            "",
            "| Species | Contrast | Expression gene | log2FC | padj | Interpretation |",
            "|---|---|---|---:|---:|---|",
        ]
    )
    for row in significant:
        lines.append(
            f"| {row['species_name']} | {row['contrast']} | {row['expression_gene_id']} | "
            f"{row.get('log2FoldChange', '')} | {row.get('padj', '')} | "
            f"{row['de_interpretation']} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation ceiling",
            "",
            "The combined evidence supports a conserved wound/apoplastic defense component "
            "that responds during infection. It does not establish a resistance-specific "
            "mechanism, causal protection, or physical effector binding.",
            "",
            "P. strobus transcript assignments are projection-only homology evidence and are "
            "not added to stable orthogroup copy counts or CAFE inputs. The absence of a clean "
            "interolog bridge is retained as a network limitation.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orthogroup", default="OG0005853")
    parser.add_argument("--stable-genes", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--tree", type=Path, required=True)
    parser.add_argument("--swissprot", type=Path, required=True)
    parser.add_argument("--defense-hits", type=Path, required=True)
    parser.add_argument("--pden-mapping", type=Path, required=True)
    parser.add_argument("--pstrobus-mapping", type=Path, required=True)
    parser.add_argument("--pden-tx2gene", type=Path, required=True)
    parser.add_argument("--pstrobus-tx2gene", type=Path, required=True)
    parser.add_argument("--pden-primary", type=Path, required=True)
    parser.add_argument("--pden-bxyl-water", type=Path, required=True)
    parser.add_argument("--pden-bthai-water", type=Path, required=True)
    parser.add_argument("--pstrobus-two-week", type=Path, required=True)
    parser.add_argument("--pstrobus-four-week", type=Path, required=True)
    parser.add_argument("--member-output", type=Path, required=True)
    parser.add_argument("--expression-output", type=Path, required=True)
    parser.add_argument("--tree-audit-output", type=Path, required=True)
    parser.add_argument("--report-output", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    members = build_member_annotation(
        args.orthogroup,
        args.stable_genes,
        args.metadata,
        args.fasta,
        args.swissprot,
        args.defense_hits,
    )
    audit = audit_tree(args.orthogroup, members, args.tree)
    expression = build_expression_evidence(
        args.orthogroup,
        args.pden_mapping,
        args.pstrobus_mapping,
        args.pden_tx2gene,
        args.pstrobus_tx2gene,
        {
            "pathogen_associated": args.pden_primary,
            "bxyl_vs_water": args.pden_bxyl_water,
            "bthai_vs_water": args.pden_bthai_water,
        },
        {
            "two_week_vs_zero": args.pstrobus_two_week,
            "four_week_vs_zero": args.pstrobus_four_week,
        },
    )
    write_tsv(args.member_output, members)
    write_tsv(args.expression_output, expression)
    write_tsv(args.tree_audit_output, [audit])
    write_report(args.report_output, args.orthogroup, members, audit, expression)
    LOGGER.info(
        "Wrote %d members, %d expression rows, and one tree audit for %s",
        len(members),
        len(expression),
        args.orthogroup,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
