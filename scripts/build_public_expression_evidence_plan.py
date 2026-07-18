#!/usr/bin/env python3
"""Map mechanism evidence gaps to public expression evidence sources."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


PLAN_FIELDS = [
    "mechanism_axis",
    "next_evidence_type",
    "gap_priority",
    "candidate_orthogroups",
    "source_id",
    "source_role",
    "source_priority",
    "species_or_panel",
    "evidence_type",
    "url",
    "extraction_goal",
    "notes",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_public_expression_evidence_plan(
    gaps_path: Path,
    sources_path: Path,
    output_path: Path,
    markdown_path: Path,
) -> list[dict[str, str]]:
    gaps = read_tsv(gaps_path)
    sources = read_tsv(sources_path)
    plan_rows = build_plan_rows(gaps, sources)
    write_tsv(output_path, PLAN_FIELDS, plan_rows)
    write_markdown(markdown_path, plan_rows)
    return plan_rows


def build_plan_rows(gaps: list[dict[str, str]], sources: list[dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for gap in gaps:
        if gap.get("next_evidence_type") == "infection_expression_support":
            rows.extend(build_expression_source_rows(gap, sources))
        else:
            rows.append(build_gap_marker_row(gap))
    return rows


def build_expression_source_rows(
    gap: dict[str, str],
    sources: list[dict[str, str]],
) -> list[dict[str, str]]:
    axis = gap.get("mechanism_axis", "")
    matched_sources = [
        source for source in sources if axis in split_multi(source.get("route_relevance", ""))
    ]
    matched_sources = sorted(matched_sources, key=source_sort_key)
    if not matched_sources:
        return [build_no_source_row(gap)]
    return [build_expression_row(gap, source) for source in matched_sources]


def build_expression_row(gap: dict[str, str], source: dict[str, str]) -> dict[str, str]:
    axis = gap.get("mechanism_axis", "")
    candidates = gap.get("candidate_orthogroups", "")
    return {
        "mechanism_axis": axis,
        "next_evidence_type": gap.get("next_evidence_type", ""),
        "gap_priority": gap.get("priority", ""),
        "candidate_orthogroups": candidates,
        "source_id": source.get("source_id", ""),
        "source_role": "candidate_expression_source",
        "source_priority": source.get("priority", ""),
        "species_or_panel": source.get("species_or_panel", ""),
        "evidence_type": source.get("evidence_type", ""),
        "url": source.get("url", ""),
        "extraction_goal": build_expression_goal(axis, candidates),
        "notes": source.get("notes", ""),
    }


def build_no_source_row(gap: dict[str, str]) -> dict[str, str]:
    return {
        "mechanism_axis": gap.get("mechanism_axis", ""),
        "next_evidence_type": gap.get("next_evidence_type", ""),
        "gap_priority": gap.get("priority", ""),
        "candidate_orthogroups": gap.get("candidate_orthogroups", ""),
        "source_id": "no_matching_expression_source",
        "source_role": "gap_marker_no_expression_source",
        "source_priority": "",
        "species_or_panel": "",
        "evidence_type": "",
        "url": "",
        "extraction_goal": "Find additional public infection-expression or resistant/susceptible expression data.",
        "notes": gap.get("why_needed", ""),
    }


def build_gap_marker_row(gap: dict[str, str]) -> dict[str, str]:
    return {
        "mechanism_axis": gap.get("mechanism_axis", ""),
        "next_evidence_type": gap.get("next_evidence_type", ""),
        "gap_priority": gap.get("priority", ""),
        "candidate_orthogroups": gap.get("candidate_orthogroups", ""),
        "source_id": "target_network_curation_needed",
        "source_role": "gap_marker_not_expression_source",
        "source_priority": "",
        "species_or_panel": "cross-system literature and future interaction data",
        "evidence_type": "target_or_network_support",
        "url": "",
        "extraction_goal": (
            "Identify direct effector-target, protein-interaction, or network-centrality evidence; "
            "do not treat expression support alone as target validation."
        ),
        "notes": gap.get("why_needed", ""),
    }


def build_expression_goal(axis: str, candidates: str) -> str:
    if axis == "effector_target_network":
        return (
            "Test whether immune-signaling and ROS candidate orthogroups respond during PWN infection "
            f"or resistant/susceptible contrasts: {candidates}."
        )
    if axis == "hydraulic_xylem_collapse":
        return (
            "Test whether cell-wall, xylem, wound, lignin, and ROS candidates show infection-timecourse "
            f"responses consistent with wilt progression: {candidates}."
        )
    return (
        "Test whether prioritized vulnerability candidates show reproducible infection or phenotype-contrast "
        f"expression support: {candidates}."
    )


def source_sort_key(source: dict[str, str]) -> tuple[int, str]:
    priority_rank = {"high": 0, "medium": 1, "low": 2}
    return (priority_rank.get(source.get("priority", ""), 99), source.get("source_id", ""))


def split_multi(value: str) -> set[str]:
    return {item for item in value.replace(",", ";").split(";") if item}


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Public Expression Evidence Plan",
        "",
        "This plan maps the route 1 + route 2 + route 3 evidence gaps to public expression sources.",
        "Rows marked as target/network curation are gap markers, not expression-derived validation.",
        "",
    ]
    for axis in sorted({row["mechanism_axis"] for row in rows}):
        axis_rows = [row for row in rows if row["mechanism_axis"] == axis]
        lines.extend([f"## {axis}", ""])
        for row in axis_rows:
            lines.append(
                f"- {row['next_evidence_type']} / {row['source_id']} "
                f"({row['source_role']}, {row['source_priority'] or row['gap_priority']}): "
                f"{row['candidate_orthogroups']}"
            )
            lines.append(f"  Goal: {row['extraction_goal']}")
        lines.append("")
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gaps", type=Path, required=True)
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = build_public_expression_evidence_plan(
        gaps_path=args.gaps,
        sources_path=args.sources,
        output_path=args.output,
        markdown_path=args.markdown,
    )
    print(f"Public expression evidence plan rows written: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
