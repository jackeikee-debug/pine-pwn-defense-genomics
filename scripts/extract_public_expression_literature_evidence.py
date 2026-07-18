#!/usr/bin/env python3
"""Extract keyword-level expression evidence from cached public literature text."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


FIELDS = [
    "source_id",
    "mechanism_axis",
    "orthogroup_id",
    "matched_keyword",
    "module_ids",
    "manuscript_section",
    "species_or_panel",
    "evidence_type",
    "support_status",
    "hit_count",
    "example_context",
    "source_url",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def extract_public_expression_literature_evidence(
    candidates_path: Path,
    evidence_plan_path: Path,
    text_dir: Path,
    output_path: Path,
    markdown_path: Path,
) -> list[dict[str, str]]:
    candidates = {row["orthogroup_id"]: row for row in read_tsv(candidates_path)}
    plan_rows = [
        row for row in read_tsv(evidence_plan_path)
        if row.get("source_role") == "candidate_expression_source"
    ]
    rows: list[dict[str, str]] = []
    for plan in plan_rows:
        rows.extend(extract_source_rows(plan, candidates, text_dir))
    write_tsv(output_path, rows)
    write_markdown(markdown_path, rows)
    return rows


def extract_source_rows(
    plan: dict[str, str],
    candidates: dict[str, dict[str, str]],
    text_dir: Path,
) -> list[dict[str, str]]:
    text_path = text_dir / f"{plan['source_id']}.txt"
    if not text_path.exists():
        text = ""
        missing = True
    else:
        text = trim_after_reference_sections(text_path.read_text(encoding="utf-8", errors="replace"))
        missing = False
    sentences = split_sentences(text)
    rows: list[dict[str, str]] = []
    for orthogroup_id in split_multi(plan.get("candidate_orthogroups", "")):
        candidate = candidates.get(orthogroup_id, {})
        keyword = candidate.get("matched_keywords", "")
        hits = [] if missing else find_keyword_sentences(sentences, keyword)
        rows.append(
            {
                "source_id": plan.get("source_id", ""),
                "mechanism_axis": plan.get("mechanism_axis", ""),
                "orthogroup_id": orthogroup_id,
                "matched_keyword": keyword,
                "module_ids": candidate.get("module_ids", ""),
                "manuscript_section": candidate.get("manuscript_section", ""),
                "species_or_panel": plan.get("species_or_panel", ""),
                "evidence_type": plan.get("evidence_type", ""),
                "support_status": classify_status(missing, hits),
                "hit_count": str(len(hits)),
                "example_context": hits[0] if hits else "",
                "source_url": plan.get("url", ""),
            }
        )
    return rows


def find_keyword_sentences(sentences: list[str], keyword: str) -> list[str]:
    if not keyword:
        return []
    pattern = re.compile(rf"\b{re.escape(keyword)}\b", flags=re.IGNORECASE)
    return [sentence for sentence in sentences if pattern.search(sentence)]


def split_sentences(text: str) -> list[str]:
    cleaned = re.sub(r"\s+", " ", text).strip()
    if not cleaned:
        return []
    return [
        sentence.strip()
        for sentence in re.split(r"(?<=[.!?])\s+", cleaned)
        if sentence.strip()
    ]


def trim_after_reference_sections(text: str) -> str:
    marker_patterns = [
        r"\bReferences\b",
        r"\bAssociated Data\b",
        r"\bSupplementary Information\b",
    ]
    matches = [re.search(pattern, text) for pattern in marker_patterns]
    starts = [match.start() for match in matches if match]
    if not starts:
        return text
    return text[: min(starts)].strip()


def classify_status(missing: bool, hits: list[str]) -> str:
    if missing:
        return "source_text_missing"
    if hits:
        return "literature_keyword_hit"
    return "no_keyword_hit"


def split_multi(value: str) -> list[str]:
    return [item for item in value.replace(",", ";").split(";") if item]


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    hit_rows = [row for row in rows if row["support_status"] == "literature_keyword_hit"]
    lines = [
        "# Public Expression Literature Evidence",
        "",
        "This table reports keyword-level evidence from cached public literature text.",
        "It is an evidence index, not a quantitative expression analysis or direct target validation.",
        "",
        f"- Candidate-source rows: {len(rows)}",
        f"- Literature keyword hits: {len(hit_rows)}",
        "",
    ]
    for row in hit_rows[:30]:
        lines.append(
            f"- {row['source_id']} / {row['mechanism_axis']} / {row['orthogroup_id']} "
            f"({row['matched_keyword']}): {row['example_context']}"
        )
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--evidence-plan", type=Path, required=True)
    parser.add_argument("--text-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = extract_public_expression_literature_evidence(
        candidates_path=args.candidates,
        evidence_plan_path=args.evidence_plan,
        text_dir=args.text_dir,
        output_path=args.output,
        markdown_path=args.markdown,
    )
    hit_count = sum(row["support_status"] == "literature_keyword_hit" for row in rows)
    print(f"Expression literature candidate-source rows written: {len(rows)}; keyword hits: {hit_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
