#!/usr/bin/env python3
"""Build interolog/PPI input rows for current Tier A and Tier B orthogroups."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


ARATH_RE = re.compile(r"\|([A-Z0-9]+_ARATH)\b")
GENE_RE = re.compile(r"\bGN=([A-Za-z0-9_.-]+)")
FIELDS = [
    "orthogroup_id", "evidence_tier", "module_ids", "mechanism_axes", "feature_ids",
    "pine_annotation_gene_ids", "arabidopsis_swissprot_ids", "arabidopsis_homolog_symbols",
    "expression_support", "ppi_eligibility", "claim_ceiling",
]
CLAIM = "Predicted Arabidopsis interolog/PPI evidence only; not direct pine interaction, effector targeting, or causal resistance evidence"


def _read(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_network_input(shortlist_path: Path, hits_path: Path, validation_path: Path, output_path: Path) -> list[dict[str, str]]:
    hits: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in _read(hits_path):
        hits[row["orthogroup_id"]].append(row)
    support: dict[str, set[str]] = defaultdict(set)
    for row in _read(validation_path):
        classification = row.get("evidence_classification", "")
        if classification not in {"unavailable", "not_significant", ""}:
            support[row["orthogroup_id"]].add(f"{row['species']}:{row['contrast']}:{classification}")

    rows = []
    for candidate in sorted(_read(shortlist_path), key=lambda row: row["orthogroup_id"]):
        if candidate["evidence_tier"] not in {"Tier A", "Tier B"}:
            continue
        annotations = hits.get(candidate["orthogroup_id"], [])
        arath_ids = sorted({match for row in annotations for match in ARATH_RE.findall(row.get("annotation_text", ""))})
        symbols = sorted({symbol for row in annotations if ARATH_RE.search(row.get("annotation_text", "")) for symbol in GENE_RE.findall(row.get("annotation_text", ""))})
        rows.append({
            "orthogroup_id": candidate["orthogroup_id"],
            "evidence_tier": candidate["evidence_tier"],
            "module_ids": candidate.get("module_ids", ""),
            "mechanism_axes": candidate.get("mechanism_axes", ""),
            "feature_ids": candidate.get("feature_ids", ""),
            "pine_annotation_gene_ids": ";".join(sorted({row.get("gene_id", "") for row in annotations if row.get("gene_id")})),
            "arabidopsis_swissprot_ids": ";".join(arath_ids),
            "arabidopsis_homolog_symbols": ";".join(symbols),
            "expression_support": ";".join(sorted(support.get(candidate["orthogroup_id"], set()))),
            "ppi_eligibility": "arabidopsis_interolog_ready" if arath_ids else "annotation_gap",
            "claim_ceiling": CLAIM,
        })
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with Path(output_path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shortlist", type=Path, required=True)
    parser.add_argument("--hits", type=Path, required=True)
    parser.add_argument("--validation", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = build_network_input(args.shortlist, args.hits, args.validation, args.output)
    print(f"Wrote {len(rows)} mechanism network candidates; {sum(row['ppi_eligibility'] == 'arabidopsis_interolog_ready' for row in rows)} interolog-ready")


if __name__ == "__main__":
    main()
