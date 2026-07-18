#!/usr/bin/env python3
"""Integrate InterPro, MEROPS, and dbCAN evidence for the standard secretome."""

from __future__ import annotations

import argparse
import csv
import logging
import re
from collections import defaultdict
from pathlib import Path


INTERPRO_COLUMNS = [
    "protein_id",
    "sequence_md5",
    "sequence_length",
    "analysis",
    "signature_accession",
    "signature_description",
    "start",
    "stop",
    "score",
    "status",
    "date",
    "interpro_accession",
    "interpro_description",
    "go_annotations",
    "pathways",
]
NON_FUNCTIONAL_ANALYSES = {"AntiFam", "MobiDBLite"}
MEROPS_FAMILY_RE = re.compile(r"\]#([A-Z][A-Z0-9]+)#")


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or "protein_id" not in rows[0]:
        raise ValueError("Sequence manifest must contain protein_id")
    identifiers = [row["protein_id"] for row in rows]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("Sequence manifest contains duplicate protein identifiers")
    return rows


def parse_interpro(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if not 13 <= len(fields) <= len(INTERPRO_COLUMNS):
                raise ValueError(
                    f"InterPro TSV line {line_number} has {len(fields)} fields; expected 13-15"
                )
            fields.extend(["-"] * (len(INTERPRO_COLUMNS) - len(fields)))
            rows.append(dict(zip(INTERPRO_COLUMNS, fields)))
    return rows


def parse_merops(
    path: Path,
    max_evalue: float = 1e-10,
    min_identity: float = 25.0,
    min_subject_coverage: float = 50.0,
) -> list[dict[str, str]]:
    hits: list[dict[str, str]] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t", maxsplit=10)
            if len(fields) != 11:
                raise ValueError(f"MEROPS TSV line {line_number} must have 11 fields")
            qseqid, sseqid, pident, length, qlen, slen, evalue, bitscore, qcov, scov, stitle = fields
            family_match = MEROPS_FAMILY_RE.search(stitle)
            family = family_match.group(1) if family_match else "unresolved"
            accepted = (
                float(evalue) <= max_evalue
                and float(pident) >= min_identity
                and float(scov) >= min_subject_coverage
            )
            if accepted:
                hits.append(
                    {
                        "protein_id": qseqid,
                        "merops_id": sseqid,
                        "merops_family": family,
                        "pident": pident,
                        "alignment_length": length,
                        "query_length": qlen,
                        "subject_length": slen,
                        "evalue": evalue,
                        "bitscore": bitscore,
                        "query_coverage": qcov,
                        "subject_coverage": scov,
                        "title": stitle,
                    }
                )
    return hits


def parse_dbcan(
    path: Path,
    max_i_evalue: float = 1e-15,
    min_hmm_coverage: float = 0.35,
) -> list[dict[str, str]]:
    hits: list[dict[str, str]] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split(maxsplit=22)
            if len(fields) < 22:
                raise ValueError(f"dbCAN domtblout line {line_number} has too few fields")
            model_length = int(fields[2])
            hmm_from = int(fields[15])
            hmm_to = int(fields[16])
            i_evalue = float(fields[12])
            coverage = (hmm_to - hmm_from + 1) / model_length
            if i_evalue <= max_i_evalue and coverage >= min_hmm_coverage:
                family = re.sub(r"\.hmm$", "", fields[0])
                hits.append(
                    {
                        "protein_id": fields[3],
                        "dbcan_family": family,
                        "i_evalue": fields[12],
                        "domain_score": fields[13],
                        "hmm_coverage": f"{coverage:.6f}",
                        "hmm_from": fields[15],
                        "hmm_to": fields[16],
                    }
                )
    return hits


def join_values(values: set[str]) -> str:
    return ";".join(sorted(value for value in values if value and value != "-")) or "NA"


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def integrate_functional_evidence(
    manifest_path: Path,
    interpro_path: Path,
    merops_path: Path,
    dbcan_path: Path,
    interpro_hits_path: Path,
    output_path: Path,
    report_path: Path,
) -> list[dict[str, str]]:
    manifest = read_manifest(manifest_path)
    expected_ids = {row["protein_id"] for row in manifest}
    interpro_hits = parse_interpro(interpro_path)
    merops_hits = parse_merops(merops_path)
    dbcan_hits = parse_dbcan(dbcan_path)

    observed_ids = {
        row["protein_id"] for row in interpro_hits + merops_hits + dbcan_hits
    }
    unexpected = sorted(observed_ids - expected_ids)
    if unexpected:
        raise ValueError(f"Unexpected query identifiers in functional evidence: {', '.join(unexpected[:10])}")

    interpro_by_id: dict[str, list[dict[str, str]]] = defaultdict(list)
    merops_by_id: dict[str, list[dict[str, str]]] = defaultdict(list)
    dbcan_by_id: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in interpro_hits:
        interpro_by_id[row["protein_id"]].append(row)
    for row in merops_hits:
        merops_by_id[row["protein_id"]].append(row)
    for row in dbcan_hits:
        dbcan_by_id[row["protein_id"]].append(row)

    output_rows: list[dict[str, str]] = []
    for manifest_row in manifest:
        protein_id = manifest_row["protein_id"]
        ipr = interpro_by_id[protein_id]
        functional_ipr = [
            row
            for row in ipr
            if row["analysis"] not in NON_FUNCTIONAL_ANALYSES
            and row["signature_accession"] not in {"", "-"}
        ]
        merops = merops_by_id[protein_id]
        dbcan = dbcan_by_id[protein_id]
        domain_support = bool(functional_ipr)
        support_count = sum((domain_support, bool(merops), bool(dbcan)))
        evidence_tier = (
            "multi_source_functional_support"
            if support_count >= 2
            else "single_source_functional_support"
            if support_count == 1
            else "sequence_defined_only"
        )
        claim_ceiling = (
            "functionally annotated candidate secreted-soluble protein"
            if support_count
            else "candidate secreted-soluble protein"
        )
        output_rows.append(
            {
                "protein_id": protein_id,
                "sequence_length": manifest_row.get("sequence_length", "NA"),
                "interpro_hit_count": str(len(ipr)),
                "interpro_accessions": join_values({row["interpro_accession"] for row in ipr}),
                "pfam_accessions": join_values(
                    {row["signature_accession"] for row in ipr if row["analysis"] == "Pfam"}
                ),
                "signature_applications": join_values({row["analysis"] for row in ipr}),
                "domain_support": "yes" if domain_support else "no",
                "merops_hit_count": str(len(merops)),
                "merops_families": join_values({row["merops_family"] for row in merops}),
                "merops_support": "yes" if merops else "no",
                "dbcan_hit_count": str(len(dbcan)),
                "dbcan_families": join_values({row["dbcan_family"] for row in dbcan}),
                "dbcan_support": "yes" if dbcan else "no",
                "functional_support_count": str(support_count),
                "evidence_tier": evidence_tier,
                "claim_ceiling": claim_ceiling,
            }
        )

    write_tsv(interpro_hits_path, interpro_hits, INTERPRO_COLUMNS)
    write_tsv(output_path, output_rows, list(output_rows[0]))
    tier_counts: dict[str, int] = defaultdict(int)
    for row in output_rows:
        tier_counts[row["evidence_tier"]] += 1
    interpro_proteins = sum(int(row["interpro_hit_count"]) > 0 for row in output_rows)
    pfam_proteins = sum(row["pfam_accessions"] != "NA" for row in output_rows)
    domain_proteins = sum(row["domain_support"] == "yes" for row in output_rows)
    merops_proteins = sum(row["merops_support"] == "yes" for row in output_rows)
    dbcan_proteins = sum(row["dbcan_support"] == "yes" for row in output_rows)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        "# Nematode standard-secretome functional evidence\n\n"
        f"- Standard secreted-soluble proteins: {len(output_rows)}\n"
        f"- InterPro-signature proteins: {interpro_proteins}\n"
        f"- Pfam-supported proteins: {pfam_proteins}\n"
        f"- Functional-domain-supported proteins: {domain_proteins}\n"
        f"- MEROPS-supported proteins: {merops_proteins}\n"
        f"- dbCAN-supported proteins: {dbcan_proteins}\n"
        f"- Multi-source functional support: {tier_counts['multi_source_functional_support']}\n"
        f"- Single-source functional support: {tier_counts['single_source_functional_support']}\n"
        f"- Sequence-defined only: {tier_counts['sequence_defined_only']}\n\n"
        "MEROPS support requires E-value <= 1e-10, identity >= 25%, and "
        "peptidase-unit subject coverage >= 50%. dbCAN support requires "
        "domain i-Evalue <= 1e-15 and HMM coverage >= 0.35. AntiFam and "
        "MobiDBLite matches are retained but do not count as functional-domain support.\n\n"
        "These annotations support candidate prioritization only; they do not establish "
        "secretion, enzymatic activity, effector function, direct host targeting, or causality.\n",
        encoding="utf-8",
    )
    logging.info("Integrated functional evidence for %d proteins", len(output_rows))
    return output_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--interpro", type=Path, required=True)
    parser.add_argument("--merops", type=Path, required=True)
    parser.add_argument("--dbcan", type=Path, required=True)
    parser.add_argument("--interpro-hits", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--log", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    handlers: list[logging.Handler] = [logging.StreamHandler()]
    if args.log:
        args.log.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(args.log, encoding="utf-8"))
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=handlers,
    )
    integrate_functional_evidence(
        args.manifest,
        args.interpro,
        args.merops,
        args.dbcan,
        args.interpro_hits,
        args.output,
        args.report,
    )


if __name__ == "__main__":
    main()
