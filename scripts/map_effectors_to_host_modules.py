#!/usr/bin/env python3
"""Map nematode effector classes to prioritized host defense modules."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


FIELDS = [
    "effector_class",
    "host_module",
    "link_type",
    "evidence_basis",
    "link_confidence",
    "effector_count",
    "host_candidate_count",
    "effector_ids",
    "host_orthogroups",
    "host_bias_directions",
    "host_matched_keywords",
    "effector_evidence",
    "notes",
]
CONFIDENCE_ORDER = {"low": 0, "medium": 1, "high": 2}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def map_effectors_to_host_modules(
    effectors_path: Path,
    mapping_path: Path,
    host_shortlist_path: Path,
    output_path: Path,
    min_confidence: str = "medium",
) -> list[dict[str, str]]:
    effectors_by_class = collect_effectors(effectors_path, min_confidence)
    host_by_module = collect_host_candidates(host_shortlist_path)
    rows: list[dict[str, str]] = []

    for mapping in read_tsv(mapping_path):
        effector_class = mapping["effector_class"]
        host_module = mapping["host_module"]
        effectors = effectors_by_class.get(effector_class)
        host_candidates = host_by_module.get(host_module)
        if not effectors or not host_candidates:
            continue
        rows.append(
            {
                "effector_class": effector_class,
                "host_module": host_module,
                "link_type": mapping["link_type"],
                "evidence_basis": mapping["evidence_basis"],
                "link_confidence": "candidate",
                "effector_count": str(len(effectors["ids"])),
                "host_candidate_count": str(len(host_candidates["orthogroups"])),
                "effector_ids": join_set(effectors["ids"]),
                "host_orthogroups": join_set(host_candidates["orthogroups"]),
                "host_bias_directions": join_set(host_candidates["bias_directions"]),
                "host_matched_keywords": join_set(host_candidates["matched_keywords"]),
                "effector_evidence": join_set(effectors["evidence"]),
                "notes": mapping.get("notes", ""),
            }
        )

    write_tsv(output_path, FIELDS, rows)
    return rows


def collect_effectors(path: Path, min_confidence: str) -> dict[str, dict[str, set[str]]]:
    threshold = CONFIDENCE_ORDER[min_confidence]
    grouped: dict[str, dict[str, set[str]]] = defaultdict(lambda: {"ids": set(), "evidence": set()})
    for row in read_tsv(path):
        if CONFIDENCE_ORDER.get(row.get("confidence", "low"), 0) < threshold:
            continue
        group = grouped[row["effector_class"]]
        group["ids"].add(row["protein_id"])
        group["evidence"].add(row["evidence"])
    return grouped


def collect_host_candidates(path: Path) -> dict[str, dict[str, set[str]]]:
    grouped: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: {"orthogroups": set(), "bias_directions": set(), "matched_keywords": set()}
    )
    for row in read_tsv(path):
        group = grouped[row["module_id"]]
        group["orthogroups"].add(row["orthogroup_id"])
        group["bias_directions"].add(row["copy_bias_direction"])
        for keyword in split_values(row.get("matched_keywords", "")):
            group["matched_keywords"].add(keyword)
    return grouped


def split_values(value: str) -> list[str]:
    if not value.strip():
        return []
    for separator in [";", ","]:
        if separator in value:
            return [part.strip() for part in value.split(separator) if part.strip()]
    return [value.strip()]


def join_set(values: set[str]) -> str:
    return ";".join(sorted(values))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--effectors", type=Path, required=True)
    parser.add_argument("--mapping", type=Path, default=Path("config/effector_host_module_map.tsv"))
    parser.add_argument("--host-shortlist", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--min-confidence", default="medium", choices=sorted(CONFIDENCE_ORDER))
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = map_effectors_to_host_modules(
        effectors_path=args.effectors,
        mapping_path=args.mapping,
        host_shortlist_path=args.host_shortlist,
        output_path=args.output,
        min_confidence=args.min_confidence,
    )
    print(f"Effector-host module links written: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
