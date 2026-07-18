#!/usr/bin/env python3
"""Build a CAFE-compatible ultrametric Newick tree from node-age metadata."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


SUMMARY_FIELDS = [
    "root_node",
    "root_age_mya",
    "leaf_count",
    "internal_node_count",
    "expected_species_count",
    "is_ultrametric",
    "max_root_to_tip_delta",
    "tree_note",
]


class TreeNode:
    def __init__(
        self,
        node_id: str,
        parent_id: str,
        node_type: str,
        label: str,
        age_mya: float,
        child_order: int,
    ) -> None:
        self.node_id = node_id
        self.parent_id = parent_id
        self.node_type = node_type
        self.label = label
        self.age_mya = age_mya
        self.child_order = child_order


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def prepare_ultrametric_tree(
    calibrations_path: Path,
    manifest_path: Path,
    output_path: Path,
    summary_path: Path,
) -> list[dict[str, str]]:
    nodes = parse_nodes(read_tsv(calibrations_path))
    expected_species = [row["species_id"] for row in read_tsv(manifest_path)]
    validate_species(nodes, expected_species)

    root = find_root(nodes)
    children = build_children(nodes)
    validate_ages(nodes, children)

    newick = build_newick(root.node_id, nodes, children)
    write_text(output_path, f"{newick};\n")

    distances = root_to_tip_distances(root.node_id, nodes, children, 0.0)
    max_delta = max(distances.values()) - min(distances.values()) if distances else 0.0
    summary_rows = [
        {
            "root_node": root.node_id,
            "root_age_mya": format_number(root.age_mya),
            "leaf_count": str(sum(1 for node in nodes.values() if node.node_type == "leaf")),
            "internal_node_count": str(sum(1 for node in nodes.values() if node.node_type != "leaf")),
            "expected_species_count": str(len(expected_species)),
            "is_ultrametric": str(max_delta < 1e-6).lower(),
            "max_root_to_tip_delta": format_number(max_delta),
            "tree_note": "Pilot ultrametric tree assembled from config/cafe_time_calibrations.tsv; replace or refine calibrations before formal inference.",
        }
    ]
    write_tsv(summary_path, SUMMARY_FIELDS, summary_rows)
    return summary_rows


def parse_nodes(rows: list[dict[str, str]]) -> dict[str, TreeNode]:
    nodes: dict[str, TreeNode] = {}
    for index, row in enumerate(rows, start=2):
        node_id = row.get("node_id", "").strip()
        if not node_id:
            raise ValueError(f"Missing node_id in calibration table line {index}")
        if node_id in nodes:
            raise ValueError(f"Duplicate node_id in calibration table: {node_id}")
        nodes[node_id] = TreeNode(
            node_id=node_id,
            parent_id=row.get("parent_id", "").strip(),
            node_type=row.get("node_type", "").strip() or "internal",
            label=row.get("label", "").strip(),
            age_mya=float(row.get("age_mya", "").strip()),
            child_order=int(row.get("child_order", "0").strip() or "0"),
        )
    return nodes


def validate_species(nodes: dict[str, TreeNode], expected_species: list[str]) -> None:
    leaves = {node.label for node in nodes.values() if node.node_type == "leaf"}
    missing = sorted(species for species in expected_species if species not in leaves)
    extra = sorted(leaf for leaf in leaves if leaf not in set(expected_species))
    if missing:
        raise ValueError(f"Species missing from calibration leaves: {', '.join(missing)}")
    if extra:
        raise ValueError(f"Calibration leaves not present in manifest: {', '.join(extra)}")


def find_root(nodes: dict[str, TreeNode]) -> TreeNode:
    roots = [node for node in nodes.values() if not node.parent_id]
    if len(roots) != 1:
        raise ValueError(f"Expected one root node, found {len(roots)}")
    return roots[0]


def build_children(nodes: dict[str, TreeNode]) -> dict[str, list[TreeNode]]:
    children: dict[str, list[TreeNode]] = {node_id: [] for node_id in nodes}
    for node in nodes.values():
        if not node.parent_id:
            continue
        if node.parent_id not in nodes:
            raise ValueError(f"Parent node not found for {node.node_id}: {node.parent_id}")
        children[node.parent_id].append(node)
    for node_children in children.values():
        node_children.sort(key=lambda node: (node.child_order, node.node_id))
    return children


def validate_ages(nodes: dict[str, TreeNode], children: dict[str, list[TreeNode]]) -> None:
    for parent_id, node_children in children.items():
        parent = nodes[parent_id]
        for child in node_children:
            if child.age_mya > parent.age_mya:
                raise ValueError(f"Node {child.node_id} is older than parent {parent.node_id}")


def build_newick(node_id: str, nodes: dict[str, TreeNode], children: dict[str, list[TreeNode]]) -> str:
    node = nodes[node_id]
    node_children = children[node_id]
    if node.node_type == "leaf":
        return sanitize_label(node.label)
    if not node_children:
        raise ValueError(f"Internal node has no children: {node_id}")

    child_parts = []
    for child in node_children:
        branch_length = node.age_mya - child.age_mya
        child_parts.append(f"{build_newick(child.node_id, nodes, children)}:{format_number(branch_length)}")
    return f"({','.join(child_parts)})"


def root_to_tip_distances(
    node_id: str,
    nodes: dict[str, TreeNode],
    children: dict[str, list[TreeNode]],
    distance: float,
) -> dict[str, float]:
    node = nodes[node_id]
    if node.node_type == "leaf":
        return {node.label: distance}
    distances: dict[str, float] = {}
    for child in children[node_id]:
        distances.update(
            root_to_tip_distances(
                child.node_id,
                nodes,
                children,
                distance + nodes[node_id].age_mya - child.age_mya,
            )
        )
    return distances


def sanitize_label(label: str) -> str:
    if not label:
        raise ValueError("Leaf node is missing label")
    return label.replace(" ", "_")


def format_number(value: float) -> str:
    formatted = f"{value:.6f}".rstrip("0").rstrip(".")
    return formatted or "0"


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calibrations", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    summary_rows = prepare_ultrametric_tree(
        calibrations_path=args.calibrations,
        manifest_path=args.manifest,
        output_path=args.output,
        summary_path=args.summary,
    )
    print(f"Ultrametric tree written: {args.output}")
    print(f"Summary rows written: {len(summary_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
