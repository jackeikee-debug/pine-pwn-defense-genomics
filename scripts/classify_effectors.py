#!/usr/bin/env python3
"""Classify candidate nematode effectors from predicted secretome features."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


FIELDS = [
    "protein_id",
    "effector_class",
    "evidence",
    "confidence",
    "annotation_text",
    "length",
    "cysteine_fraction",
]


def read_effector_keywords(path: Path) -> dict[str, list[str]]:
    classes: dict[str, list[str]] = {}
    current_class: str | None = None
    in_keywords = False
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#") or stripped == "effector_classes:":
            continue
        if raw_line.startswith("  ") and not raw_line.startswith("    ") and stripped.endswith(":"):
            current_class = stripped[:-1]
            classes[current_class] = []
            in_keywords = False
            continue
        if current_class and stripped == "keywords:":
            in_keywords = True
            continue
        if current_class and in_keywords and stripped.startswith("- "):
            classes[current_class].append(stripped[2:].strip().strip("'\""))
            continue
        if current_class and raw_line.startswith("    ") and not raw_line.startswith("      "):
            in_keywords = False
    return classes


def classify_effectors(
    secretome_path: Path,
    keywords_path: Path,
    output_path: Path,
    annotations_path: Path | None = None,
) -> list[dict[str, str]]:
    effector_keywords = read_effector_keywords(keywords_path)
    annotations = read_annotations(annotations_path)
    rows: list[dict[str, str]] = []
    with secretome_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("secretome_candidate") != "yes":
                continue
            annotation_text = combined_annotation_text(row, annotations)
            matched = keyword_classes(annotation_text, effector_keywords)
            for effector_class, keyword in matched:
                rows.append(build_row(row, effector_class, f"keyword:{keyword};secretome:{row['evidence']}", "medium", annotation_text))
            if not matched and row.get("small_secreted_candidate") == "yes" and row.get("cysteine_rich_candidate") == "yes":
                rows.append(build_row(row, "small_cysteine_rich_secreted", f"secretome:{row['evidence']}", "low", annotation_text))
            elif not matched and row.get("small_secreted_candidate") == "yes":
                rows.append(build_row(row, "small_secreted", f"secretome:{row['evidence']}", "low", annotation_text))

    write_tsv(output_path, FIELDS, rows)
    return rows


def read_annotations(path: Path | None) -> dict[str, str]:
    if path is None or not path.exists():
        return {}
    annotations: dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            annotations[row["gene_id"]] = f"{row['annotation_source']}: {row['annotation_text']}"
    return annotations


def combined_annotation_text(row: dict[str, str], annotations: dict[str, str]) -> str:
    base_text = row.get("annotation_text", "")
    external_text = annotations.get(row["protein_id"])
    if external_text:
        return f"{base_text} | {external_text}" if base_text else external_text
    return base_text


def keyword_classes(annotation_text: str, effector_keywords: dict[str, list[str]]) -> list[tuple[str, str]]:
    lower_text = annotation_text.lower()
    matches: list[tuple[str, str]] = []
    for effector_class, keywords in effector_keywords.items():
        for keyword in keywords:
            if keyword.lower() in lower_text:
                matches.append((effector_class, keyword))
                break
    return matches


def build_row(
    row: dict[str, str],
    effector_class: str,
    evidence: str,
    confidence: str,
    annotation_text: str,
) -> dict[str, str]:
    return {
        "protein_id": row["protein_id"],
        "effector_class": effector_class,
        "evidence": evidence,
        "confidence": confidence,
        "annotation_text": annotation_text,
        "length": row.get("length", ""),
        "cysteine_fraction": row.get("cysteine_fraction", ""),
    }


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--secretome", type=Path, required=True)
    parser.add_argument("--keywords", type=Path, default=Path("config/effector_keywords.yaml"))
    parser.add_argument("--annotations", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = classify_effectors(args.secretome, args.keywords, args.output, args.annotations)
    print(f"Candidate effector rows written: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
