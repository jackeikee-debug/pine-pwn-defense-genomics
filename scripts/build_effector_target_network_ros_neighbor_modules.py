#!/usr/bin/env python3
"""Classify FSD1/FSD3 STRING neighbors into route 1 ROS network modules."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


FIELDS = [
    "orthogroup_id",
    "seed_symbol",
    "partner_symbol",
    "functional_module",
    "effector_target_network_relevance",
    "combined_score",
    "experimental_score",
    "database_score",
    "textmining_score",
    "edge_channel_tier",
    "seed_channel_support_tier",
    "effector_classes",
    "expression_support_level",
    "evidence_scope",
]

SUMMARY_FIELDS = [
    "functional_module",
    "neighbor_count",
    "experimental_or_database_supported_count",
    "mean_combined_score",
    "top_neighbors",
]

EVIDENCE_SCOPE = (
    "Heuristic functional classification of Arabidopsis STRING neighbors for the "
    "route 1 SOD/FSD subnetwork; not pine-specific physical validation"
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_effector_target_network_ros_neighbor_modules(
    channels_path: Path,
    interactions_path: Path,
    output_path: Path,
    summary_path: Path,
    markdown_path: Path,
    channel_threshold: float,
) -> list[dict[str, str]]:
    fsd_seeds = [
        row for row in read_tsv(channels_path)
        if row.get("matched_keyword", "").lower() == "superoxide dismutase"
        and row.get("mapping_status") == "mapped"
    ]
    interaction_index = index_interactions(read_tsv(interactions_path))
    rows = []
    for seed in fsd_seeds:
        for interaction in interaction_index.get(seed["string_id"], []):
            rows.append(build_neighbor_row(seed, interaction, channel_threshold))
    rows = sorted(rows, key=lambda row: (row["seed_symbol"], row["functional_module"], row["partner_symbol"]))
    summary_rows = summarize_modules(rows)
    write_tsv(output_path, rows, FIELDS)
    write_tsv(summary_path, summary_rows, SUMMARY_FIELDS)
    write_markdown(markdown_path, rows, summary_rows)
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


def with_partner(row: dict[str, str], partner_side: str) -> dict[str, str]:
    partner_id = row.get(f"stringId_{partner_side}", "")
    partner_name = row.get(f"preferredName_{partner_side}", partner_id)
    item = dict(row)
    item["partner_id"] = partner_id
    item["partner_name"] = partner_name
    return item


def add_best_partner(partners: dict[str, dict[str, str]], item: dict[str, str]) -> None:
    partner_id = item.get("partner_id", "")
    previous = partners.get(partner_id)
    if previous is None or parse_float(item.get("score", "0")) > parse_float(previous.get("score", "0")):
        partners[partner_id] = item


def build_neighbor_row(seed: dict[str, str], interaction: dict[str, str], channel_threshold: float) -> dict[str, str]:
    partner = interaction.get("partner_name", "")
    module, relevance = classify_partner(partner)
    return {
        "orthogroup_id": seed.get("orthogroup_id", ""),
        "seed_symbol": seed.get("query_symbol", ""),
        "partner_symbol": partner,
        "functional_module": module,
        "effector_target_network_relevance": relevance,
        "combined_score": format_score(interaction.get("score", "0")),
        "experimental_score": format_score(interaction.get("escore", "0")),
        "database_score": format_score(interaction.get("dscore", "0")),
        "textmining_score": format_score(interaction.get("tscore", "0")),
        "edge_channel_tier": edge_channel_tier(interaction, channel_threshold),
        "seed_channel_support_tier": seed.get("channel_support_tier", ""),
        "effector_classes": seed.get("effector_classes", ""),
        "expression_support_level": seed.get("expression_support_level", ""),
        "evidence_scope": EVIDENCE_SCOPE,
    }


def classify_partner(symbol: str) -> tuple[str, str]:
    upper = symbol.upper()
    if upper.startswith(("CSD", "FSD", "MSD")):
        return (
            "superoxide_dismutase_family",
            "Direct ROS detoxification neighborhood around Fe/Cu/Zn/Mn SOD enzymes.",
        )
    if upper.startswith("APX"):
        return (
            "ascorbate_peroxidase_ros_scavenging",
            "Hydrogen peroxide detoxification branch linked to oxidative stress buffering.",
        )
    if upper.startswith("CAT"):
        return (
            "catalase_ros_scavenging",
            "Catalase branch for hydrogen peroxide detoxification.",
        )
    if upper in {"CITRX"} or upper.startswith(("TRX", "PRX")):
        return (
            "thioredoxin_peroxiredoxin_redox",
            "Redox relay context adjacent to ROS homeostasis.",
        )
    if upper.startswith(("PTAC", "MRL", "FLN")):
        return (
            "chloroplast_redox_gene_expression",
            "Chloroplast transcription or plastid maintenance context linked to redox stress.",
        )
    if upper == "SPL7":
        return (
            "copper_homeostasis_sod_cofactor_context",
            "Copper homeostasis context relevant to Cu/Zn SOD activity.",
        )
    return (
        "unclassified_string_neighbor",
        "Functional relevance requires manual review before mechanistic interpretation.",
    )


def edge_channel_tier(row: dict[str, str], threshold: float) -> str:
    if parse_float(row.get("escore", "0")) >= threshold or parse_float(row.get("dscore", "0")) >= threshold:
        return "experimental_or_database_supported"
    if parse_float(row.get("ascore", "0")) >= threshold or parse_float(row.get("tscore", "0")) >= threshold:
        return "coexpression_or_textmining_only"
    return "combined_score_only"


def summarize_modules(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["functional_module"]].append(row)
    summary = []
    for module, module_rows in sorted(grouped.items()):
        scores = [parse_float(row["combined_score"]) for row in module_rows]
        exp_db = [
            row for row in module_rows
            if row["edge_channel_tier"] == "experimental_or_database_supported"
        ]
        top = sorted(
            module_rows,
            key=lambda row: (-parse_float(row["combined_score"]), row["partner_symbol"]),
        )[:10]
        summary.append(
            {
                "functional_module": module,
                "neighbor_count": str(len(module_rows)),
                "experimental_or_database_supported_count": str(len(exp_db)),
                "mean_combined_score": f"{(sum(scores) / len(scores)) if scores else 0.0:.3f}",
                "top_neighbors": ";".join(
                    f"{row['seed_symbol']}-{row['partner_symbol']}:{row['combined_score']}"
                    for row in top
                ),
            }
        )
    return summary


def format_score(value: str) -> str:
    return f"{parse_float(value):.3f}"


def parse_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(path: Path, rows: list[dict[str, str]], summary_rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    exp_db = [row for row in rows if row["edge_channel_tier"] == "experimental_or_database_supported"]
    lines = [
        "# Route 1 ROS Neighbor Modules",
        "",
        "This report classifies the FSD1/FSD3 STRING local-neighborhood partners into ROS-relevant functional modules.",
        "",
        f"- FSD neighbor edges: {len(rows)}",
        f"- Experimental/database-supported edges: {len(exp_db)}",
        f"- Functional modules: {len(summary_rows)}",
        "",
        "Evidence boundary: classifications are heuristic Arabidopsis STRING-neighborhood annotations, not pine-specific physical validation.",
        "",
        "| module | neighbors | exp/db supported | mean combined score | top neighbors |",
        "|---|---:|---:|---:|---|",
    ]
    for row in summary_rows:
        lines.append(
            "| {functional_module} | {neighbor_count} | {experimental_or_database_supported_count} | {mean_combined_score} | {top_neighbors} |".format(
                **row
            )
        )
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--channels", type=Path, required=True)
    parser.add_argument("--interactions", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    parser.add_argument("--channel-threshold", type=float, default=0.4)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = build_effector_target_network_ros_neighbor_modules(
        channels_path=args.channels,
        interactions_path=args.interactions,
        output_path=args.output,
        summary_path=args.summary,
        markdown_path=args.markdown,
        channel_threshold=args.channel_threshold,
    )
    exp_db = sum(row["edge_channel_tier"] == "experimental_or_database_supported" for row in rows)
    print(f"Route 1 ROS neighbor rows written: {len(rows)}; exp/db-supported edges: {exp_db}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
