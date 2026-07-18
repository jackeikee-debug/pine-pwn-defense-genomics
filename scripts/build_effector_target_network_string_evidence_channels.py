#!/usr/bin/env python3
"""Summarize STRING evidence channels for route 1 network seeds."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


FIELDS = [
    "orthogroup_id",
    "query_symbol",
    "network_priority",
    "matched_keyword",
    "effector_classes",
    "string_id",
    "mapping_status",
    "partner_count",
    "combined_weighted_degree",
    "experimental_partner_count",
    "database_partner_count",
    "coexpression_partner_count",
    "textmining_partner_count",
    "top_experimental_partners",
    "top_database_partners",
    "top_coexpression_partners",
    "top_textmining_partners",
    "channel_support_tier",
    "expression_support_level",
    "evidence_scope",
]

EVIDENCE_SCOPE = (
    "STRING evidence-channel decomposition for Arabidopsis local neighborhood; "
    "not pine-specific physical effector-target validation"
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_effector_target_network_string_evidence_channels(
    centrality_path: Path,
    interactions_path: Path,
    output_path: Path,
    markdown_path: Path,
    channel_threshold: float,
) -> list[dict[str, str]]:
    interaction_index = index_interactions(read_tsv(interactions_path))
    rows = [
        build_channel_row(row, interaction_index, channel_threshold)
        for row in read_tsv(centrality_path)
    ]
    write_tsv(output_path, rows)
    write_markdown(markdown_path, rows, channel_threshold)
    return rows


def index_interactions(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    index: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for row in rows:
        node_a = row.get("stringId_A", "")
        node_b = row.get("stringId_B", "")
        if node_a:
            add_best_partner(index[node_a], with_partner(row, partner_side="B"))
        if node_b:
            add_best_partner(index[node_b], with_partner(row, partner_side="A"))
    return {node: list(partners.values()) for node, partners in index.items()}


def add_best_partner(partners: dict[str, dict[str, str]], item: dict[str, str]) -> None:
    partner_id = item.get("partner_id", "")
    previous = partners.get(partner_id)
    if previous is None or parse_float(item.get("score", "0")) > parse_float(previous.get("score", "0")):
        partners[partner_id] = item


def with_partner(row: dict[str, str], partner_side: str) -> dict[str, str]:
    partner_id = row.get(f"stringId_{partner_side}", "")
    partner_name = row.get(f"preferredName_{partner_side}", partner_id)
    item = dict(row)
    item["partner_id"] = partner_id
    item["partner_name"] = partner_name
    return item


def build_channel_row(
    centrality_row: dict[str, str],
    interaction_index: dict[str, list[dict[str, str]]],
    channel_threshold: float,
) -> dict[str, str]:
    string_id = centrality_row.get("string_id", "")
    interactions = interaction_index.get(string_id, [])
    experimental = channel_hits(interactions, "escore", channel_threshold)
    database = channel_hits(interactions, "dscore", channel_threshold)
    coexpression = channel_hits(interactions, "ascore", channel_threshold)
    textmining = channel_hits(interactions, "tscore", channel_threshold)
    return {
        "orthogroup_id": centrality_row.get("orthogroup_id", ""),
        "query_symbol": centrality_row.get("query_symbol", ""),
        "network_priority": centrality_row.get("network_priority", ""),
        "matched_keyword": centrality_row.get("matched_keyword", ""),
        "effector_classes": centrality_row.get("effector_classes", ""),
        "string_id": string_id,
        "mapping_status": centrality_row.get("mapping_status", ""),
        "partner_count": centrality_row.get("partner_count", "0"),
        "combined_weighted_degree": centrality_row.get("weighted_degree_score_sum", "0.000"),
        "experimental_partner_count": str(len(experimental)),
        "database_partner_count": str(len(database)),
        "coexpression_partner_count": str(len(coexpression)),
        "textmining_partner_count": str(len(textmining)),
        "top_experimental_partners": format_partner_scores(experimental, "escore"),
        "top_database_partners": format_partner_scores(database, "dscore"),
        "top_coexpression_partners": format_partner_scores(coexpression, "ascore"),
        "top_textmining_partners": format_partner_scores(textmining, "tscore"),
        "channel_support_tier": support_tier(centrality_row, experimental, database, coexpression, textmining),
        "expression_support_level": centrality_row.get("expression_support_level", ""),
        "evidence_scope": EVIDENCE_SCOPE,
    }


def channel_hits(rows: list[dict[str, str]], field: str, threshold: float) -> list[dict[str, str]]:
    return sorted(
        [row for row in rows if parse_float(row.get(field, "0")) >= threshold],
        key=lambda row: (-parse_float(row.get(field, "0")), row.get("partner_name", "")),
    )


def support_tier(
    centrality_row: dict[str, str],
    experimental: list[dict[str, str]],
    database: list[dict[str, str]],
    coexpression: list[dict[str, str]],
    textmining: list[dict[str, str]],
) -> str:
    if centrality_row.get("mapping_status") != "mapped":
        return "unmapped"
    if parse_int(centrality_row.get("partner_count", "0")) == 0:
        return "mapped_no_neighbors"
    if experimental or database:
        return "experimental_or_database_supported"
    if coexpression or textmining:
        return "coexpression_or_textmining_only"
    return "combined_score_only"


def format_partner_scores(rows: list[dict[str, str]], field: str, limit: int = 10) -> str:
    return ";".join(
        f"{row.get('partner_name', '')}:{parse_float(row.get(field, '0')):.3f}"
        for row in rows[:limit]
    )


def parse_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def parse_int(value: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(path: Path, rows: list[dict[str, str]], channel_threshold: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tier_counts = defaultdict(int)
    for row in rows:
        tier_counts[row["channel_support_tier"]] += 1
    lines = [
        "# Route 1 STRING Evidence Channels",
        "",
        "This report decomposes route 1 STRING neighborhood support by evidence channel.",
        f"- Channel threshold: {channel_threshold:.3f}",
        f"- Seed rows: {len(rows)}",
        f"- Experimental/database supported rows: {tier_counts['experimental_or_database_supported']}",
        f"- Coexpression/textmining-only rows: {tier_counts['coexpression_or_textmining_only']}",
        f"- Mapped rows without neighbors: {tier_counts['mapped_no_neighbors']}",
        "",
        "Evidence boundary: these are STRING Arabidopsis functional association channels, not pine-specific physical effector-target validation.",
        "",
        "| orthogroup | symbol | tier | exp | db | coexp | text | top experimental | top database |",
        "|---|---|---|---:|---:|---:|---:|---|---|",
    ]
    for row in rows:
        lines.append(
            "| {orthogroup_id} | {query_symbol} | {channel_support_tier} | {experimental_partner_count} | {database_partner_count} | {coexpression_partner_count} | {textmining_partner_count} | {top_experimental_partners} | {top_database_partners} |".format(
                **row
            )
        )
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--centrality", type=Path, required=True)
    parser.add_argument("--interactions", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    parser.add_argument("--channel-threshold", type=float, default=0.4)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = build_effector_target_network_string_evidence_channels(
        centrality_path=args.centrality,
        interactions_path=args.interactions,
        output_path=args.output,
        markdown_path=args.markdown,
        channel_threshold=args.channel_threshold,
    )
    supported = sum(row["channel_support_tier"] == "experimental_or_database_supported" for row in rows)
    print(f"Route 1 STRING evidence-channel rows written: {len(rows)}; exp/db supported: {supported}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
