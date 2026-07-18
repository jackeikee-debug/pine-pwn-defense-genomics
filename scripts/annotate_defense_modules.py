#!/usr/bin/env python3
"""Annotate stable orthogroups with defense-module keyword evidence."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


MATRIX_FIELDS = [
    "orthogroup_id",
    "module_id",
    "species_id",
    "gene_count",
    "gene_ids",
    "matched_keywords",
    "annotation_sources",
    "evidence_type",
]
GENE_HIT_FIELDS = [
    "orthogroup_id",
    "species_id",
    "gene_id",
    "module_id",
    "matched_keyword",
    "annotation_source",
    "annotation_text",
    "evidence_type",
]
EVIDENCE_TYPE = "keyword_annotation_match"


def read_defense_modules(path: Path) -> dict[str, list[str]]:
    modules: dict[str, list[str]] = {}
    current_module: str | None = None
    in_keywords = False

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#") or stripped == "modules:":
            continue
        if raw_line.startswith("  ") and not raw_line.startswith("    ") and stripped.endswith(":"):
            current_module = stripped[:-1]
            modules[current_module] = []
            in_keywords = False
            continue
        if current_module and stripped == "keywords:":
            in_keywords = True
            continue
        if current_module and in_keywords and stripped.startswith("- "):
            modules[current_module].append(stripped[2:].strip().strip("'\""))
            continue
        if current_module and raw_line.startswith("    ") and not raw_line.startswith("      "):
            in_keywords = False

    return modules


def read_stable_genes(path: Path) -> dict[str, dict[str, str]]:
    genes: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            genes[row["gene_id"]] = row
    return genes


def read_annotations(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def keyword_matches(keyword: str, annotation_text: str) -> bool:
    has_letter = any(character.isalpha() for character in keyword)
    if has_letter and keyword == keyword.upper():
        pattern = rf"(?<![A-Za-z0-9]){re.escape(keyword)}[0-9]*(?![A-Za-z0-9])"
        return re.search(pattern, annotation_text) is not None
    return keyword.lower() in annotation_text.lower()


def annotate_defense_modules(
    stable_genes_path: Path,
    annotations_path: Path,
    modules_path: Path,
    matrix_output: Path,
    gene_hits_output: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    modules = read_defense_modules(modules_path)
    stable_genes = read_stable_genes(stable_genes_path)
    annotations = read_annotations(annotations_path)

    hit_rows: list[dict[str, str]] = []
    aggregate: dict[tuple[str, str, str], dict[str, set[str]]] = defaultdict(
        lambda: {"gene_ids": set(), "matched_keywords": set(), "annotation_sources": set()}
    )

    for annotation in annotations:
        gene_id = annotation.get("gene_id", "").strip()
        stable_gene = stable_genes.get(gene_id)
        if not stable_gene:
            continue
        annotation_text = annotation.get("annotation_text", "").strip()
        annotation_source = annotation.get("annotation_source", "").strip() or "unspecified"

        for module_id, keywords in modules.items():
            for keyword in keywords:
                if not keyword_matches(keyword, annotation_text):
                    continue
                hit_row = {
                    "orthogroup_id": stable_gene["orthogroup_id"],
                    "species_id": stable_gene["species_id"],
                    "gene_id": gene_id,
                    "module_id": module_id,
                    "matched_keyword": keyword,
                    "annotation_source": annotation_source,
                    "annotation_text": annotation_text,
                    "evidence_type": EVIDENCE_TYPE,
                }
                hit_rows.append(hit_row)
                key = (stable_gene["orthogroup_id"], module_id, stable_gene["species_id"])
                aggregate[key]["gene_ids"].add(gene_id)
                aggregate[key]["matched_keywords"].add(keyword)
                aggregate[key]["annotation_sources"].add(annotation_source)

    matrix_rows = []
    for orthogroup_id, module_id, species_id in sorted(aggregate):
        values = aggregate[(orthogroup_id, module_id, species_id)]
        gene_ids = sorted(values["gene_ids"])
        matrix_rows.append(
            {
                "orthogroup_id": orthogroup_id,
                "module_id": module_id,
                "species_id": species_id,
                "gene_count": str(len(gene_ids)),
                "gene_ids": ",".join(gene_ids),
                "matched_keywords": ";".join(sorted(values["matched_keywords"])),
                "annotation_sources": ";".join(sorted(values["annotation_sources"])),
                "evidence_type": EVIDENCE_TYPE,
            }
        )

    write_tsv(matrix_output, MATRIX_FIELDS, matrix_rows)
    write_tsv(gene_hits_output, GENE_HIT_FIELDS, hit_rows)
    return matrix_rows, hit_rows


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stable-genes", type=Path, required=True)
    parser.add_argument("--annotations", type=Path, required=True)
    parser.add_argument("--modules", type=Path, required=True)
    parser.add_argument("--matrix-output", type=Path, required=True)
    parser.add_argument("--gene-hits-output", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    matrix_rows, hit_rows = annotate_defense_modules(
        stable_genes_path=args.stable_genes,
        annotations_path=args.annotations,
        modules_path=args.modules,
        matrix_output=args.matrix_output,
        gene_hits_output=args.gene_hits_output,
    )
    print(f"Defense module matrix rows written: {len(matrix_rows)}")
    print(f"Defense module gene hits written: {len(hit_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
