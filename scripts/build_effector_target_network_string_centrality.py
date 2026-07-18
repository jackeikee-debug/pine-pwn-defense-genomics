#!/usr/bin/env python3
"""Compute STRING network centrality for route 1 Arabidopsis homolog seeds."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict, deque
from pathlib import Path


FIELDS = [
    "orthogroup_id",
    "query_symbol",
    "network_priority",
    "matched_keyword",
    "host_modules",
    "effector_classes",
    "string_id",
    "preferred_name",
    "mapping_status",
    "partner_count",
    "weighted_degree_score_sum",
    "betweenness_centrality",
    "top_partners",
    "ppi_support_status",
    "expression_support_level",
    "evidence_scope",
]

EVIDENCE_SCOPE = (
    "STRING functional association; predicted interolog/PPI overlay, "
    "not validated pine effector-target interaction"
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_effector_target_network_string_centrality(
    seed_path: Path,
    mapping_path: Path,
    interactions_path: Path,
    output_path: Path,
    markdown_path: Path,
) -> list[dict[str, str]]:
    seeds = aggregate_seed_symbols(read_tsv(seed_path))
    mappings = {row["query_symbol"]: row for row in read_tsv(mapping_path)}
    interactions = read_tsv(interactions_path)
    adjacency = build_adjacency(interactions)
    betweenness = betweenness_centrality(adjacency)
    rows = [
        build_seed_row(seed, mappings.get(seed["query_symbol"]), adjacency, betweenness)
        for seed in seeds
    ]
    write_tsv(output_path, rows)
    write_markdown(markdown_path, rows)
    return rows


def aggregate_seed_symbols(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str], dict[str, set[str] | str]] = {}
    for row in rows:
        for symbol in split_multi(row.get("arabidopsis_homolog_symbols", "")):
            key = (row["orthogroup_id"], symbol)
            if key not in grouped:
                grouped[key] = {
                    "orthogroup_id": row["orthogroup_id"],
                    "query_symbol": symbol,
                    "network_priority": row.get("network_priority", row.get("evidence_tier", "")),
                    "matched_keyword": row.get("matched_keyword", ""),
                    "host_modules": set(),
                    "effector_classes": set(),
                    "expression_levels": set(),
                }
            grouped[key]["host_modules"].add(row.get("host_module", row.get("module_ids", "")))
            grouped[key]["effector_classes"].add(row.get("effector_class", row.get("mechanism_axes", "")))
            grouped[key]["expression_levels"].add(row.get("expression_support_level", row.get("expression_support", "")))
    seeds = []
    for item in grouped.values():
        seeds.append(
            {
                "orthogroup_id": str(item["orthogroup_id"]),
                "query_symbol": str(item["query_symbol"]),
                "network_priority": str(item["network_priority"]),
                "matched_keyword": str(item["matched_keyword"]),
                "host_modules": join_sorted(item["host_modules"]),
                "effector_classes": join_sorted(item["effector_classes"]),
                "expression_support_level": join_sorted(item["expression_levels"]),
            }
        )
    return sorted(seeds, key=lambda row: (row["network_priority"], row["orthogroup_id"], row["query_symbol"]))


def build_adjacency(rows: list[dict[str, str]]) -> dict[str, dict[str, dict[str, str]]]:
    adjacency: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for row in rows:
        node_a = row.get("stringId_A", "")
        node_b = row.get("stringId_B", "")
        if not node_a or not node_b:
            continue
        score = parse_float(row.get("score", "0"))
        adjacency[node_a][node_b] = {
            "score": f"{score:.3f}",
            "partner_name": row.get("preferredName_B", node_b),
        }
        adjacency[node_b][node_a] = {
            "score": f"{score:.3f}",
            "partner_name": row.get("preferredName_A", node_a),
        }
    return adjacency


def build_seed_row(
    seed: dict[str, str],
    mapping: dict[str, str] | None,
    adjacency: dict[str, dict[str, dict[str, str]]],
    betweenness: dict[str, float],
) -> dict[str, str]:
    string_id = mapping.get("string_id", "") if mapping else ""
    mapping_status = mapping.get("mapping_status", "unmapped") if mapping else "unmapped"
    preferred_name = mapping.get("preferred_name", "") if mapping else ""
    partners = adjacency.get(string_id, {})
    top_partners = sorted(
        ((data["partner_name"], parse_float(data["score"])) for data in partners.values()),
        key=lambda item: (-item[1], item[0]),
    )
    weighted_degree = sum(score for _, score in top_partners)
    return {
        "orthogroup_id": seed["orthogroup_id"],
        "query_symbol": seed["query_symbol"],
        "network_priority": seed["network_priority"],
        "matched_keyword": seed["matched_keyword"],
        "host_modules": seed["host_modules"],
        "effector_classes": seed["effector_classes"],
        "string_id": string_id,
        "preferred_name": preferred_name,
        "mapping_status": mapping_status,
        "partner_count": str(len(partners)),
        "weighted_degree_score_sum": f"{weighted_degree:.3f}",
        "betweenness_centrality": f"{betweenness.get(string_id, 0.0):.6f}",
        "top_partners": ";".join(f"{name}:{score:.3f}" for name, score in top_partners[:10]),
        "ppi_support_status": ppi_support_status(mapping_status, len(partners)),
        "expression_support_level": seed["expression_support_level"],
        "evidence_scope": EVIDENCE_SCOPE,
    }


def ppi_support_status(mapping_status: str, partner_count: int) -> str:
    if mapping_status != "mapped":
        return "unmapped_no_string_network"
    if partner_count == 0:
        return "mapped_no_string_neighbors_in_subnetwork"
    return "string_neighbors_found"


def betweenness_centrality(adjacency: dict[str, dict[str, dict[str, str]]]) -> dict[str, float]:
    nodes = sorted(adjacency)
    centrality = {node: 0.0 for node in nodes}
    for source in nodes:
        stack: list[str] = []
        predecessors = {node: [] for node in nodes}
        sigma = dict.fromkeys(nodes, 0.0)
        distance = dict.fromkeys(nodes, -1)
        sigma[source] = 1.0
        distance[source] = 0
        queue: deque[str] = deque([source])
        while queue:
            vertex = queue.popleft()
            stack.append(vertex)
            for neighbor in adjacency[vertex]:
                if distance[neighbor] < 0:
                    queue.append(neighbor)
                    distance[neighbor] = distance[vertex] + 1
                if distance[neighbor] == distance[vertex] + 1:
                    sigma[neighbor] += sigma[vertex]
                    predecessors[neighbor].append(vertex)
        dependency = dict.fromkeys(nodes, 0.0)
        while stack:
            node = stack.pop()
            for predecessor in predecessors[node]:
                if sigma[node]:
                    dependency[predecessor] += (sigma[predecessor] / sigma[node]) * (1 + dependency[node])
            if node != source:
                centrality[node] += dependency[node]
    if len(nodes) > 2:
        scale = 1 / ((len(nodes) - 1) * (len(nodes) - 2))
        for node in centrality:
            centrality[node] *= scale
    return centrality


def split_multi(value: str) -> list[str]:
    return [item.strip() for item in value.replace(",", ";").split(";") if item.strip()]


def join_sorted(values) -> str:
    return ";".join(sorted({value for value in values if value}))


def parse_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    mapped = [row for row in rows if row["mapping_status"] == "mapped"]
    supported = [row for row in rows if row["ppi_support_status"] == "string_neighbors_found"]
    mechanism_scope = any(row.get("network_priority", "").startswith("Tier ") for row in rows)
    title = "Mechanism Candidate" if mechanism_scope else "Route 1"
    seed_description = "current Tier A/B mechanism candidates" if mechanism_scope else "route 1 Arabidopsis homolog seeds"
    lines = [
        f"# {title} STRING Centrality",
        "",
        f"This report summarizes STRING functional-association network metrics for {seed_description}.",
        "",
        f"- Seed symbol rows: {len(rows)}",
        f"- Mapped STRING rows: {len(mapped)}",
        f"- Rows with STRING neighbors: {len(supported)}",
        "",
        "Evidence boundary: STRING scores are confidence scores for functional associations, not pine-specific effector-target validation.",
        "Centrality is computed on the fetched STRING neighborhood subnetwork, not on the complete Arabidopsis interactome.",
        "",
        "| orthogroup | symbol | effector classes | partners | weighted degree | betweenness | top partners |",
        "|---|---|---|---:|---:|---:|---|",
    ]
    for row in rows:
        lines.append(
            "| {orthogroup_id} | {query_symbol} | {effector_classes} | {partner_count} | {weighted_degree_score_sum} | {betweenness_centrality} | {top_partners} |".format(
                **row
            )
        )
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=Path, required=True)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--interactions", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = build_effector_target_network_string_centrality(
        seed_path=args.seed,
        mapping_path=args.mapping,
        interactions_path=args.interactions,
        output_path=args.output,
        markdown_path=args.markdown,
    )
    supported = sum(row["ppi_support_status"] == "string_neighbors_found" for row in rows)
    print(f"Route 1 STRING centrality rows written: {len(rows)}; supported rows: {supported}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
