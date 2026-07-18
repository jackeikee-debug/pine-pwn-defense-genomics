#!/usr/bin/env python3
"""Audit Route 2 Tier 1 proteins and draft a controlled cleavage-assay design."""

from __future__ import annotations

import argparse
import csv
import logging
import re
from collections import defaultdict
from pathlib import Path


AUDIT_FIELDS = [
    "protein_id", "protein_role", "orthogroup_id", "host_family", "sequence_length",
    "signal_peptide_end", "recommended_recombinant_region", "swissprot_subject_id",
    "swissprot_pident", "swissprot_alignment_length", "alignment_length_to_query_length_ratio",
    "expected_domain", "interpro_support", "interpro_domain_span", "catalytic_motif",
    "proposed_catalytic_dead_change", "alphafold_ranking_score", "alphafold_ptm",
    "alphafold_mean_plddt", "alphafold_catalytic_feature_mean_plddt",
    "alphafold_global_quality_tier", "interaction_log2_fold_change", "interaction_padj",
    "interaction_coefficient", "domain_completeness_status", "direct_interaction_evidence",
    "claim_ceiling",
]
ASSAY_FIELDS = [
    "effector_id", "substrate_id", "host_orthogroup_id", "host_family",
    "effector_construct", "substrate_construct", "condition", "condition_detail",
    "primary_readout", "functional_readout", "interpretation_guardrail",
    "direct_interaction_evidence",
]
CLAIM = (
    "Sequence/domain/expression-supported biochemical hypothesis only; these data do not demonstrate "
    "M8 binding, cleavage, in planta targeting, or causal contribution to pine wilt disease"
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    current = ""
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                records[current] = ""
            else:
                if not current:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                records[current] += line.upper()
    return records


def read_interpro(path: Path) -> dict[str, list[dict[str, str]]]:
    hits: dict[str, list[dict[str, str]]] = defaultdict(list)
    with path.open(encoding="utf-8", newline="") as handle:
        for fields in csv.reader(handle, delimiter="\t"):
            if not fields or fields[0].startswith("#") or len(fields) < 13:
                continue
            hits[fields[0]].append({
                "analysis": fields[3], "signature": fields[4], "description": fields[5],
                "start": fields[6], "end": fields[7], "interpro_id": fields[11],
                "interpro_description": fields[12],
            })
    return hits


def read_signal_ends(path: Path) -> dict[str, str]:
    ends: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 5 and fields[2].lower() == "signal peptide":
                ends[fields[0]] = fields[4]
    return ends


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def motif_for(role: str, family: str, sequence: str, domain_start: int, domain_end: int) -> tuple[str, str]:
    if role == "candidate_effector":
        matches = list(re.finditer(r"HE..H", sequence))
        match = next((m for m in matches if m.start() + 1 >= domain_start and m.end() <= domain_end), None)
        if match:
            start = match.start() + 1
            return f"{match.group()}@{start}-{match.end()}", f"E{start + 1}A_proposed"
        return "canonical_HEXXH_not_found", "not_available"
    if family == "XTH":
        match = re.search(r"DEID.EFL", sequence)
        if match:
            return f"{match.group()}@{match.start() + 1}-{match.end()}", "not_applicable"
        return "canonical_XTH_motif_not_found", "not_applicable"
    return "not_applicable_non_enzyme_candidate", "not_applicable"


def build_outputs(
    fasta_path: Path,
    interpro_path: Path,
    shortlist_path: Path,
    validation_path: Path,
    effector_deepsig_path: Path,
    nematode_hits_path: Path,
    pine_hits_path: Path,
    alphafold_path: Path,
    interaction_path: Path,
    audit_path: Path,
    assay_path: Path,
    report_path: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    sequences = read_fasta(fasta_path)
    interpro = read_interpro(interpro_path)
    shortlist = [row for row in read_tsv(shortlist_path) if row["sequence_resolved_tier"] == "Tier 1_sequence_supported"]
    validation = {row["gene_id"]: row for row in read_tsv(validation_path)}
    signal_ends = read_signal_ends(effector_deepsig_path)
    swiss = {row["gene_id"]: row for row in read_tsv(nematode_hits_path) + read_tsv(pine_hits_path)}
    alphafold = {row["protein_id"]: row for row in read_tsv(alphafold_path)}
    interaction = {row["feature_id"]: row for row in read_tsv(interaction_path)}

    effectors = sorted({row["effector_id"] for row in shortlist})
    metadata: dict[str, dict[str, str]] = {}
    for effector in effectors:
        metadata[effector] = {"role": "candidate_effector", "orthogroup": "", "family": "M8"}
    for row in shortlist:
        metadata[row["expression_selected_gene_id"]] = {
            "role": "candidate_substrate", "orthogroup": row["host_orthogroup_id"],
            "family": row["host_family"],
        }

    audit: list[dict[str, str]] = []
    expected_terms = {
        "M8": ("Peptidase M8/leishmanolysin", ("leishmanolysin", "peptidase m8")),
        "PMEI": ("Pectinesterase inhibitor", ("pectinesterase inhibitor", "pmei")),
        "XTH": ("GH16 XTH plus XET C-terminal", ("xyloglucan endotrans", "gh16_xet", "glycosyl hydrolases family 16")),
    }
    for protein_id, meta in metadata.items():
        sequence = sequences[protein_id]
        family = meta["family"]
        expected, keywords = expected_terms[family]
        relevant = [
            hit for hit in interpro.get(protein_id, [])
            if any(keyword in (hit["description"] + " " + hit["interpro_description"]).lower() for keyword in keywords)
        ]
        starts = [int(hit["start"]) for hit in relevant]
        ends = [int(hit["end"]) for hit in relevant]
        domain_start = min(starts) if starts else 1
        domain_end = max(ends) if ends else len(sequence)
        support = ";".join(sorted({
            f'{hit["analysis"]}:{hit["signature"]}:{hit["description"]}' for hit in relevant
        }))
        motif, mutant = motif_for(meta["role"], family, sequence, domain_start, domain_end)
        if meta["role"] == "candidate_effector":
            signal_end = signal_ends.get(protein_id, "")
        else:
            signal_end = validation.get(protein_id, {}).get("deepsig_signal_end", "")
        mature_start = int(signal_end) + 1 if signal_end.isdigit() else 1
        if meta["role"] == "candidate_effector":
            construct = f"residues_{mature_start}-{len(sequence)}_after_predicted_signal_peptide; activation_boundary_unresolved"
        else:
            construct = f"residues_{mature_start}-{len(sequence)}_after_predicted_signal_peptide"
        hit = swiss.get(protein_id, {})
        ratio = float(hit["alignment_length"]) / len(sequence) if hit.get("alignment_length") else 0.0
        af = alphafold.get(protein_id, {})
        feature_id = protein_id.split("|", 1)[-1]
        expr = interaction.get(feature_id, {})
        audit.append({
            "protein_id": protein_id, "protein_role": meta["role"], "orthogroup_id": meta["orthogroup"],
            "host_family": family, "sequence_length": str(len(sequence)), "signal_peptide_end": signal_end,
            "recommended_recombinant_region": construct, "swissprot_subject_id": hit.get("subject_id", ""),
            "swissprot_pident": hit.get("pident", ""), "swissprot_alignment_length": hit.get("alignment_length", ""),
            "alignment_length_to_query_length_ratio": f"{ratio:.3f}", "expected_domain": expected,
            "interpro_support": support, "interpro_domain_span": f"{domain_start}-{domain_end}" if relevant else "",
            "catalytic_motif": motif, "proposed_catalytic_dead_change": mutant,
            "alphafold_ranking_score": af.get("ranking_score", ""), "alphafold_ptm": af.get("ptm", ""),
            "alphafold_mean_plddt": af.get("mean_ca_plddt", ""),
            "alphafold_catalytic_feature_mean_plddt": af.get("catalytic_feature_mean_plddt", ""),
            "alphafold_global_quality_tier": af.get("global_quality_tier", ""),
            "interaction_log2_fold_change": expr.get("interaction_log2_fold_change", ""),
            "interaction_padj": expr.get("interaction_padj", ""),
            "interaction_coefficient": expr.get("interaction_coefficient", ""),
            "domain_completeness_status": "expected_domain_supported" if relevant else "expected_domain_not_confirmed",
            "direct_interaction_evidence": "none", "claim_ceiling": CLAIM,
        })
    audit.sort(key=lambda row: (row["protein_role"] != "candidate_effector", row["host_family"], row["protein_id"]))

    effector_row = next(row for row in audit if row["protein_role"] == "candidate_effector")
    assay: list[dict[str, str]] = []
    conditions = [
        ("wild_type_M8", "Purified wild-type M8 ectodomain plus substrate; optimize time, pH, and metal conditions"),
        ("catalytic_dead_M8", f'Proposed {effector_row["proposed_catalytic_dead_change"]} M8 plus substrate; verify the residue assignment before synthesis'),
        ("heat_inactivated_M8", "Heat-inactivated wild-type M8 plus substrate"),
        ("metalloprotease_inhibitor", "Wild-type M8 plus substrate with a 1,10-phenanthroline inhibition control"),
        ("substrate_only", "Substrate incubated in the matched buffer without M8"),
    ]
    for substrate in (row for row in audit if row["protein_role"] == "candidate_substrate"):
        functional = (
            "Residual inhibition of a validated plant pectin methylesterase"
            if substrate["host_family"] == "PMEI"
            else "Residual xyloglucan endotransglucosylase/hydrolase activity"
        )
        for condition, detail in conditions:
            assay.append({
                "effector_id": effector_row["protein_id"], "substrate_id": substrate["protein_id"],
                "host_orthogroup_id": substrate["orthogroup_id"], "host_family": substrate["host_family"],
                "effector_construct": effector_row["recommended_recombinant_region"],
                "substrate_construct": substrate["recommended_recombinant_region"],
                "condition": condition, "condition_detail": detail,
                "primary_readout": "Time-course SDS-PAGE/immunoblot followed by intact-mass LC-MS and cleavage-site mapping for positive samples",
                "functional_readout": functional,
                "interpretation_guardrail": "Specific cleavage requires reproducibility with wild-type M8 and loss or reduction in catalytic-dead/inhibitor controls; in vitro cleavage alone does not establish in planta targeting or disease causality",
                "direct_interaction_evidence": "none",
            })

    write_tsv(audit_path, audit, AUDIT_FIELDS)
    write_tsv(assay_path, assay, ASSAY_FIELDS)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        "# Route 2 Tier 1 domain audit and experimental design\n\n"
        f"- Proteins audited: {len(audit)} (one candidate M8 effector and {len(audit) - 1} candidate substrates).\n"
        f"- Expected domains supported by InterProScan: {sum(r['domain_completeness_status'] == 'expected_domain_supported' for r in audit)}/{len(audit)}.\n"
        f"- Cleavage-assay conditions drafted: {len(assay)}.\n\n"
        "The M8 HEXXH and XTH catalytic motifs are sequence features within supported domains. The proposed M8 "
        "catalytic-dead substitution must be checked against a curated M8 active-site alignment before construct "
        "synthesis. The predicted signal-peptide boundary is suitable for an initial ectodomain construct, but the "
        "M8 activation/propeptide boundary remains unresolved.\n\n"
        "Expression interaction coefficients prioritize these substrates biologically; they do not show that M8 "
        "caused the expression pattern or directly contacted either pine protein.\n",
        encoding="utf-8",
    )
    logging.info("Wrote %d audit rows and %d assay rows", len(audit), len(assay))
    return audit, assay


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--interpro", type=Path, required=True)
    parser.add_argument("--shortlist", type=Path, required=True)
    parser.add_argument("--validation", type=Path, required=True)
    parser.add_argument("--effector-deepsig", type=Path, required=True)
    parser.add_argument("--nematode-hits", type=Path, required=True)
    parser.add_argument("--pine-hits", type=Path, required=True)
    parser.add_argument("--alphafold", type=Path, required=True)
    parser.add_argument("--interaction", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--assay", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    build_outputs(
        args.fasta, args.interpro, args.shortlist, args.validation, args.effector_deepsig,
        args.nematode_hits, args.pine_hits, args.alphafold, args.interaction,
        args.audit, args.assay, args.report,
    )


if __name__ == "__main__":
    main()
