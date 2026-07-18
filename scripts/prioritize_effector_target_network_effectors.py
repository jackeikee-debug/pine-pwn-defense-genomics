#!/usr/bin/env python3
"""Prioritize functionally annotated route 1 nematode effector candidates."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


FIELDS = [
    "protein_id", "effector_class", "candidate_confidence", "functional_support_level",
    "swissprot_subject_id", "pident", "alignment_length", "query_coverage", "evalue",
    "bitscore", "swissprot_annotation", "secretome_evidence", "length", "cysteine_fraction",
    "structure_followup_status", "claim_ceiling",
]
effector_target_network_CLASSES = {"protease", "detoxification", "mimicry"}
CLAIM = "Computational secretome and homology evidence only; not experimentally validated secretion, effector activity, or host targeting"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def support_level(row: dict[str, str], length: float) -> tuple[str, float]:
    coverage = float(row["alignment_length"]) / length if length else 0.0
    pident, evalue, bitscore = float(row["pident"]), float(row["evalue"]), float(row["bitscore"])
    if evalue <= 1e-20 and bitscore >= 100 and coverage >= 0.5 and pident >= 30:
        return "strong_swissprot_alignment", coverage
    if evalue <= 1e-5 and bitscore >= 50 and coverage >= 0.25:
        return "moderate_swissprot_alignment", coverage
    return "weak_swissprot_alignment", coverage


def prioritize(effectors_path: Path, hits_path: Path, output_path: Path) -> list[dict[str, str]]:
    hits = {row["gene_id"]: row for row in read_tsv(hits_path)}
    rows = []
    for effector in read_tsv(effectors_path):
        if effector.get("confidence") not in {"medium", "high"} or effector.get("effector_class") not in effector_target_network_CLASSES:
            continue
        hit = hits.get(effector["protein_id"])
        if not hit:
            continue
        length = float(effector.get("length", "0") or 0)
        level, coverage = support_level(hit, length)
        rows.append({
            "protein_id": effector["protein_id"],
            "effector_class": effector["effector_class"],
            "candidate_confidence": effector["confidence"],
            "functional_support_level": level,
            "swissprot_subject_id": hit["subject_id"],
            "pident": hit["pident"],
            "alignment_length": hit["alignment_length"],
            "query_coverage": f"{coverage:.4f}",
            "evalue": hit["evalue"],
            "bitscore": hit["bitscore"],
            "swissprot_annotation": hit["annotation_text"],
            "secretome_evidence": effector["evidence"],
            "length": effector["length"],
            "cysteine_fraction": effector["cysteine_fraction"],
            "structure_followup_status": "include_in_effector_target_network_sequence_clustering",
            "claim_ceiling": CLAIM,
        })
    rows.sort(key=lambda row: (row["effector_class"], row["protein_id"]))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--effectors", type=Path, required=True)
    parser.add_argument("--hits", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = prioritize(args.effectors, args.hits, args.output)
    print(f"Wrote {len(rows)} route 1 effector candidates")


if __name__ == "__main__":
    main()
