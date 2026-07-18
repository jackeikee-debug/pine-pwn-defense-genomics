#!/usr/bin/env python3
"""Build a mechanism-axis summary for higher-impact hypothesis development."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


AXIS_FIELDS = [
    "mechanism_axis",
    "route_label",
    "working_hypothesis",
    "primary_candidate_orthogroups",
    "supporting_candidate_orthogroups",
    "deprioritized_candidate_orthogroups",
    "host_modules",
    "effector_classes",
    "effector_link_count",
    "candidate_status",
    "claim_ceiling",
]

GAP_FIELDS = [
    "mechanism_axis",
    "next_evidence_type",
    "why_needed",
    "candidate_orthogroups",
    "priority",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_mechanism_axis_summary(
    candidates_path: Path,
    links_path: Path,
    axes_path: Path,
    axis_output_path: Path,
    gaps_output_path: Path,
    markdown_path: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    candidates = read_tsv(candidates_path)
    links = read_tsv(links_path)
    axes = read_axes(axes_path)
    axis_rows = [build_axis_row(axis, spec, candidates, links) for axis, spec in axes.items()]
    gap_rows = build_gap_rows(axis_rows)
    write_tsv(axis_output_path, AXIS_FIELDS, axis_rows)
    write_tsv(gaps_output_path, GAP_FIELDS, gap_rows)
    write_markdown(markdown_path, axis_rows, gap_rows)
    return axis_rows, gap_rows


def read_axes(path: Path) -> dict[str, dict[str, object]]:
    rows = read_tsv(path)
    axes: dict[str, dict[str, object]] = {}
    for row in rows:
        axis = row["mechanism_axis"]
        axes[axis] = {
            "route_label": row["route_label"],
            "modules": split_multi(row.get("host_modules", "")),
            "keywords": split_multi(row.get("candidate_keywords", "")),
            "hypothesis": row["working_hypothesis"],
            "claim": row["claim_ceiling"],
        }
    return axes


def build_axis_row(
    axis: str,
    spec: dict[str, object],
    candidates: list[dict[str, str]],
    links: list[dict[str, str]],
) -> dict[str, str]:
    matched_candidates = [row for row in candidates if candidate_matches_axis(row, axis, spec)]
    primary = [
        row["orthogroup_id"]
        for row in matched_candidates
        if row.get("manuscript_section") == "main_text_priority_review"
    ]
    supporting = [
        row["orthogroup_id"]
        for row in matched_candidates
        if row.get("manuscript_section") != "main_text_priority_review"
    ]
    deprioritized = [
        row["orthogroup_id"]
        for row in matched_candidates
        if row.get("manuscript_section") == "supplement_background"
    ]
    matched_links = [row for row in links if link_matches_axis(row, axis, spec)]
    return {
        "mechanism_axis": axis,
        "route_label": str(spec["route_label"]),
        "working_hypothesis": str(spec["hypothesis"]),
        "primary_candidate_orthogroups": join_unique(primary),
        "supporting_candidate_orthogroups": join_unique(supporting),
        "deprioritized_candidate_orthogroups": join_unique(deprioritized),
        "host_modules": join_unique(row.get("host_module", "") for row in matched_links),
        "effector_classes": join_unique(row.get("effector_class", "") for row in matched_links),
        "effector_link_count": str(len(matched_links)),
        "candidate_status": classify_axis_status(primary, supporting, deprioritized),
        "claim_ceiling": str(spec["claim"]),
    }


def candidate_matches_axis(row: dict[str, str], axis: str, spec: dict[str, object]) -> bool:
    if axis == "susceptibility_vulnerability":
        return True
    modules = split_multi(row.get("module_ids", ""))
    keywords = {row.get("matched_keywords", "")}
    return bool(modules & spec["modules"] or keywords & spec["keywords"])


def link_matches_axis(row: dict[str, str], axis: str, spec: dict[str, object]) -> bool:
    if axis == "susceptibility_vulnerability":
        return True
    modules = {row.get("host_module", "")}
    keywords = split_multi(row.get("host_matched_keywords", ""))
    return bool(modules & spec["modules"] or keywords & spec["keywords"])


def split_multi(value: str) -> set[str]:
    return {item for item in value.replace(",", ";").split(";") if item}


def join_unique(values) -> str:
    return ";".join(sorted({value for value in values if value}))


def classify_axis_status(primary: list[str], supporting: list[str], deprioritized: list[str]) -> str:
    if primary:
        return "main_text_hypothesis_with_priority_candidates"
    if supporting:
        return "supporting_hypothesis_only"
    if deprioritized:
        return "background_only"
    return "no_current_candidate"


def build_gap_rows(axis_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for row in axis_rows:
        candidates = join_unique(
            split_multi(row["primary_candidate_orthogroups"])
            | split_multi(row["supporting_candidate_orthogroups"])
            | split_multi(row["deprioritized_candidate_orthogroups"])
        )
        rows.extend(
            [
                {
                    "mechanism_axis": row["mechanism_axis"],
                    "next_evidence_type": "infection_expression_support",
                    "why_needed": "Test whether candidate modules respond during PWN infection or resistant/susceptible contrasts.",
                    "candidate_orthogroups": candidates,
                    "priority": "high",
                },
                {
                    "mechanism_axis": row["mechanism_axis"],
                    "next_evidence_type": "target_or_network_support",
                    "why_needed": "Distinguish effector-guided network vulnerability from co-occurrence by annotation.",
                    "candidate_orthogroups": candidates,
                    "priority": "high" if row["mechanism_axis"] == "effector_target_network" else "medium",
                },
            ]
        )
    return rows


def write_markdown(path: Path, axis_rows: list[dict[str, str]], gap_rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Mechanism Strategy",
        "",
        "This route 1 + route 2 + route 3 strategy reframes the project as a mechanism-generating framework.",
        "The rows below are hypotheses and evidence gaps, not validated mechanisms.",
        "",
    ]
    for row in axis_rows:
        lines.extend(
            [
                f"## {row['mechanism_axis']} ({row['route_label']})",
                "",
                f"- Working hypothesis: {row['working_hypothesis']}",
                f"- Primary candidates: {row['primary_candidate_orthogroups'] or 'none'}",
                f"- Supporting candidates: {row['supporting_candidate_orthogroups'] or 'none'}",
                f"- Deprioritized candidates: {row['deprioritized_candidate_orthogroups'] or 'none'}",
                f"- Effector classes: {row['effector_classes'] or 'none'}",
                f"- Claim ceiling: {row['claim_ceiling']}",
                "",
            ]
        )
    lines.extend(["## Evidence gaps", ""])
    for row in gap_rows:
        lines.append(
            f"- {row['mechanism_axis']} / {row['next_evidence_type']} ({row['priority']}): {row['why_needed']}"
        )
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--links", type=Path, required=True)
    parser.add_argument("--axes", type=Path, required=True)
    parser.add_argument("--axis-output", type=Path, required=True)
    parser.add_argument("--gaps-output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    axis_rows, gap_rows = build_mechanism_axis_summary(
        candidates_path=args.candidates,
        links_path=args.links,
        axes_path=args.axes,
        axis_output_path=args.axis_output,
        gaps_output_path=args.gaps_output,
        markdown_path=args.markdown,
    )
    print(f"Mechanism axes written: {len(axis_rows)}; evidence gaps written: {len(gap_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
