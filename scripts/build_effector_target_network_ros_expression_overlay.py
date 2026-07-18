#!/usr/bin/env python3
"""Overlay route 1 ROS neighbor modules with cached public expression text."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


FIELDS = [
    "source_id",
    "functional_module",
    "neighbor_symbols",
    "module_keywords",
    "support_status",
    "hit_count",
    "matched_terms",
    "example_context",
    "source_url",
    "source_error",
    "claim_ceiling",
]

MODULE_KEYWORDS = {
    "superoxide_dismutase_family": ["superoxide dismutase", "SOD"],
    "ascorbate_peroxidase_ros_scavenging": ["ascorbate peroxidase", "APX"],
    "catalase_ros_scavenging": ["catalase", "CAT"],
    "chloroplast_redox_gene_expression": ["chloroplast", "plastid", "redox", "oxidative stress"],
    "copper_homeostasis_sod_cofactor_context": ["copper", "SPL7"],
    "thioredoxin_peroxiredoxin_redox": ["thioredoxin", "peroxiredoxin", "redox"],
    "unclassified_string_neighbor": ["oxidative stress", "reactive oxygen", "ROS"],
}

CLAIM_CEILING = (
    "module-level infection-expression text context only; not exact pine gene DEG "
    "or validated effector-target evidence"
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_effector_target_network_ros_expression_overlay(
    modules_path: Path,
    manifest_path: Path,
    text_dir: Path,
    output_path: Path,
    markdown_path: Path,
) -> list[dict[str, str]]:
    modules = module_terms(read_tsv(modules_path))
    rows = []
    for source in read_tsv(manifest_path):
        text, missing = load_source_text(source, text_dir)
        sentences = split_sentences(trim_after_reference_sections(text))
        for module_name, terms in modules.items():
            rows.append(build_overlay_row(source, module_name, terms, sentences, missing))
    write_tsv(output_path, rows)
    write_markdown(markdown_path, rows)
    return rows


def module_terms(rows: list[dict[str, str]]) -> dict[str, dict[str, list[str]]]:
    grouped: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        module_name = row.get("functional_module", "")
        partner = row.get("partner_symbol", "")
        if module_name and partner:
            grouped[module_name].add(partner)
    modules = {}
    for module_name, partners in grouped.items():
        keyword_terms = MODULE_KEYWORDS.get(module_name, [])
        modules[module_name] = {
            "neighbor_symbols": sorted(partners),
            "module_keywords": keyword_terms,
            "search_terms": sorted(set(keyword_terms) | set(partners)),
        }
    return dict(sorted(modules.items()))


def load_source_text(source: dict[str, str], text_dir: Path) -> tuple[str, bool]:
    if source.get("status") != "ok":
        return "", True
    text_path = Path(source.get("text_path", ""))
    if not text_path.exists():
        text_path = text_dir / f"{source.get('source_id', '')}.txt"
    if not text_path.exists():
        return "", True
    return text_path.read_text(encoding="utf-8", errors="replace"), False


def build_overlay_row(
    source: dict[str, str],
    module_name: str,
    terms: dict[str, list[str]],
    sentences: list[str],
    missing: bool,
) -> dict[str, str]:
    hits, matched_terms = find_hits(sentences, terms["search_terms"]) if not missing else ([], [])
    return {
        "source_id": source.get("source_id", ""),
        "functional_module": module_name,
        "neighbor_symbols": ";".join(terms["neighbor_symbols"]),
        "module_keywords": ";".join(terms["module_keywords"]),
        "support_status": classify_status(missing, hits),
        "hit_count": str(len(hits)),
        "matched_terms": ";".join(matched_terms),
        "example_context": hits[0] if hits else "",
        "source_url": source.get("url", ""),
        "source_error": source.get("error", ""),
        "claim_ceiling": CLAIM_CEILING,
    }


def find_hits(sentences: list[str], terms: list[str]) -> tuple[list[str], list[str]]:
    hits = []
    matched_terms: set[str] = set()
    patterns = [(term, compile_term(term)) for term in terms if term]
    for sentence in sentences:
        sentence_terms = [term for term, pattern in patterns if pattern.search(sentence)]
        if sentence_terms:
            hits.append(sentence)
            matched_terms.update(sentence_terms)
    return hits, sorted(matched_terms, key=lambda value: value.lower())


def compile_term(term: str) -> re.Pattern:
    if re.search(r"\s", term):
        pattern = re.escape(term)
        flags = re.IGNORECASE
    elif term.isupper() and len(term) <= 4:
        pattern = rf"\b{re.escape(term)}[A-Z0-9_.-]*\b"
        flags = 0
    else:
        pattern = rf"\b{re.escape(term)}\b"
        flags = re.IGNORECASE
    return re.compile(pattern, flags=flags)


def classify_status(missing: bool, hits: list[str]) -> str:
    if missing:
        return "source_text_missing"
    if hits:
        return "module_expression_context_hit"
    return "no_module_keyword_hit"


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


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    hits = [row for row in rows if row["support_status"] == "module_expression_context_hit"]
    lines = [
        "# Route 1 ROS Expression Overlay",
        "",
        "This report overlays FSD1/FSD3 ROS neighbor modules with cached public infection-expression literature text.",
        "",
        f"- Module-source rows: {len(rows)}",
        f"- Module-source text hits: {len(hits)}",
        "",
        "Evidence boundary: hits are module-level text context, not exact pine gene DEG evidence.",
        "",
    ]
    for row in hits[:40]:
        lines.append(
            f"- {row['source_id']} / {row['functional_module']} "
            f"({row['matched_terms']}): {row['example_context']}"
        )
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--modules", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--text-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = build_effector_target_network_ros_expression_overlay(
        modules_path=args.modules,
        manifest_path=args.manifest,
        text_dir=args.text_dir,
        output_path=args.output,
        markdown_path=args.markdown,
    )
    hits = sum(row["support_status"] == "module_expression_context_hit" for row in rows)
    print(f"Route 1 ROS expression overlay rows written: {len(rows)}; text hits: {hits}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
