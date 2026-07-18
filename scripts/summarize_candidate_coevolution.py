#!/usr/bin/env python3
"""Summarize candidate effector-host module links for plotting and review."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


SUMMARY_FIELDS = [
    "module_id",
    "effector_class",
    "host_bias_direction",
    "orthogroup_count",
    "unique_effector_count",
    "mean_priority_score",
    "max_priority_score",
    "mean_bitscore",
    "top_orthogroups",
    "matched_keywords",
]

TOP_FIELDS = [
    "top_rank",
    "top_group_rank",
    "module_id",
    "effector_class",
    "host_bias_direction",
    "orthogroup_id",
    "host_priority_score",
    "host_priority_rank",
    "host_shortlist_rank",
    "mean_bitscore",
    "best_bitscore",
    "host_gene_count",
    "susceptible_gene_count",
    "tolerant_gene_count",
    "outgroup_gene_count",
    "matched_keywords",
    "effector_count",
    "effector_ids",
    "effector_confidences",
    "effector_evidence",
    "interpretation",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def summarize_candidate_coevolution(
    candidates_path: Path,
    summary_output: Path,
    top_output: Path,
    max_per_group: int = 5,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    rows = read_tsv(candidates_path)
    summary_rows = build_summary_rows(rows)
    top_rows = build_top_rows(rows, max_per_group=max_per_group)
    write_tsv(summary_output, SUMMARY_FIELDS, summary_rows)
    write_tsv(top_output, TOP_FIELDS, top_rows)
    return summary_rows, top_rows


def build_summary_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row["module_id"], row["effector_class"], row["host_bias_direction"])].append(row)

    summary_rows: list[dict[str, str]] = []
    for (module_id, effector_class, bias), group_rows in sorted(grouped.items()):
        priority_scores = [to_float(row.get("host_priority_score", "")) for row in group_rows]
        bitscores = [to_float(row.get("mean_bitscore", "")) for row in group_rows]
        effector_ids = set()
        keywords = set()
        for row in group_rows:
            effector_ids.update(split_values(row.get("effector_ids", "")))
            keywords.update(split_values(row.get("matched_keywords", "")))
        top_orthogroups = [
            row["orthogroup_id"]
            for row in sorted(group_rows, key=candidate_sort_key)[:5]
        ]
        summary_rows.append(
            {
                "module_id": module_id,
                "effector_class": effector_class,
                "host_bias_direction": bias,
                "orthogroup_count": str(len(group_rows)),
                "unique_effector_count": str(len(effector_ids)),
                "mean_priority_score": format_float(mean(priority_scores)),
                "max_priority_score": format_float(max(priority_scores) if priority_scores else 0.0),
                "mean_bitscore": format_float(mean(bitscores)),
                "top_orthogroups": ";".join(top_orthogroups),
                "matched_keywords": ";".join(sorted(keywords)),
            }
        )
    return summary_rows


def build_top_rows(rows: list[dict[str, str]], max_per_group: int) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row["module_id"], row["effector_class"], row["host_bias_direction"])].append(row)

    selected: list[dict[str, str]] = []
    for group_key in sorted(grouped):
        group_rows = sorted(grouped[group_key], key=candidate_sort_key)
        cap = max_per_group if max_per_group > 0 else len(group_rows)
        for group_rank, row in enumerate(group_rows[:cap], start=1):
            selected.append(
                {
                    "top_rank": "0",
                    "top_group_rank": str(group_rank),
                    "module_id": row.get("module_id", ""),
                    "effector_class": row.get("effector_class", ""),
                    "host_bias_direction": row.get("host_bias_direction", ""),
                    "orthogroup_id": row.get("orthogroup_id", ""),
                    "host_priority_score": row.get("host_priority_score", ""),
                    "host_priority_rank": row.get("host_priority_rank", ""),
                    "host_shortlist_rank": row.get("host_shortlist_rank", ""),
                    "mean_bitscore": row.get("mean_bitscore", ""),
                    "best_bitscore": row.get("best_bitscore", ""),
                    "host_gene_count": row.get("host_gene_count", ""),
                    "susceptible_gene_count": row.get("susceptible_gene_count", ""),
                    "tolerant_gene_count": row.get("tolerant_gene_count", ""),
                    "outgroup_gene_count": row.get("outgroup_gene_count", ""),
                    "matched_keywords": row.get("matched_keywords", ""),
                    "effector_count": row.get("effector_count", ""),
                    "effector_ids": row.get("effector_ids", ""),
                    "effector_confidences": row.get("effector_confidences", ""),
                    "effector_evidence": row.get("effector_evidence", ""),
                    "interpretation": row.get("interpretation", ""),
                }
            )

    selected.sort(key=candidate_sort_key)
    for rank, row in enumerate(selected, start=1):
        row["top_rank"] = str(rank)
    return selected


def candidate_sort_key(row: dict[str, str]) -> tuple[float, float, int, str]:
    return (
        -to_float(row.get("host_priority_score", "")),
        -to_float(row.get("mean_bitscore", "")),
        to_int(row.get("host_shortlist_rank", row.get("top_rank", ""))),
        row.get("orthogroup_id", ""),
    )


def split_values(value: str) -> list[str]:
    if not value.strip():
        return []
    for separator in [";", ","]:
        if separator in value:
            return [part.strip() for part in value.split(separator) if part.strip()]
    return [value.strip()]


def to_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def to_int(value: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return 10**9


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def format_float(value: float) -> str:
    return f"{value:.2f}"


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--summary-output", type=Path, required=True)
    parser.add_argument("--top-output", type=Path, required=True)
    parser.add_argument("--max-per-group", type=int, default=5)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    summary_rows, top_rows = summarize_candidate_coevolution(
        candidates_path=args.candidates,
        summary_output=args.summary_output,
        top_output=args.top_output,
        max_per_group=args.max_per_group,
    )
    print(f"Candidate coevolution summary rows written: {len(summary_rows)}")
    print(f"Candidate coevolution top rows written: {len(top_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
