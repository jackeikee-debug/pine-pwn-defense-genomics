#!/usr/bin/env python3
"""Extract conservative tree-clustering metrics for model-supported candidates."""

from __future__ import annotations

import argparse
import csv
from io import StringIO
from pathlib import Path

try:
    from Bio import Phylo
except ImportError:  # pragma: no cover - exercised only when Biopython is absent.
    Phylo = None


REGION_SPECIES = {
    "East_Asia": {"pden", "pmas", "ptab"},
    "North_America": {"plam", "plon", "ptae"},
}

FIELDS = [
    "check_rank",
    "orthogroup_id",
    "evidence_tier",
    "regional_direction",
    "focal_region",
    "tree_tip_count",
    "newick_parser_backend",
    "focal_region_tip_count",
    "focal_region_mrca_tip_count",
    "focal_region_mrca_purity",
    "focal_region_cluster_status",
    "focal_max_species",
    "dominant_species_tip_count",
    "dominant_species_mrca_tip_count",
    "dominant_species_mrca_purity",
    "dominant_species_cluster_status",
    "outgroup_tip_count",
    "non_focal_region_tip_count",
    "manual_tree_review_priority",
    "tree_path",
]


class Node:
    def __init__(self, name: str = "", children: list["Node"] | None = None):
        self.name = name
        self.children = children or []

    @property
    def is_leaf(self) -> bool:
        return not self.children


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_candidate_tree_metrics(checklist_path: Path, output_path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for source in read_tsv(checklist_path):
        tree_path = Path(source["tree_path"])
        tree = parse_newick(tree_path.read_text(encoding="utf-8").strip())
        leaves = leaf_names(tree)
        species = [species_id(name) for name in leaves]
        focal_species = REGION_SPECIES.get(source.get("focal_region", ""), set())
        focal_leaves = [name for name in leaves if species_id(name) in focal_species]
        dominant_species = source.get("focal_max_species", "")
        dominant_leaves = [name for name in leaves if species_id(name) == dominant_species]
        focal_mrca_leaves = mrca_leaf_names(tree, set(focal_leaves)) if focal_leaves else []
        dominant_mrca_leaves = mrca_leaf_names(tree, set(dominant_leaves)) if dominant_leaves else []
        focal_purity = purity(focal_mrca_leaves, focal_species)
        dominant_purity = purity(dominant_mrca_leaves, {dominant_species})
        outgroup_count = sum(1 for item in species if item not in REGION_SPECIES["East_Asia"] | REGION_SPECIES["North_America"])
        non_focal_count = sum(1 for item in species if item not in focal_species)
        focal_status = classify_region_cluster(focal_leaves, focal_mrca_leaves, focal_purity, len(leaves))
        dominant_status = classify_dominant_cluster(dominant_leaves, dominant_mrca_leaves, dominant_purity)
        rows.append(
            {
                "check_rank": source.get("check_rank", ""),
                "orthogroup_id": source.get("orthogroup_id", ""),
                "evidence_tier": source.get("evidence_tier", ""),
                "regional_direction": source.get("regional_direction", ""),
                "focal_region": source.get("focal_region", ""),
                "tree_tip_count": str(len(leaves)),
                "newick_parser_backend": newick_parser_backend(),
                "focal_region_tip_count": str(len(focal_leaves)),
                "focal_region_mrca_tip_count": str(len(focal_mrca_leaves)),
                "focal_region_mrca_purity": f"{focal_purity:.3f}",
                "focal_region_cluster_status": focal_status,
                "focal_max_species": dominant_species,
                "dominant_species_tip_count": str(len(dominant_leaves)),
                "dominant_species_mrca_tip_count": str(len(dominant_mrca_leaves)),
                "dominant_species_mrca_purity": f"{dominant_purity:.3f}",
                "dominant_species_cluster_status": dominant_status,
                "outgroup_tip_count": str(outgroup_count),
                "non_focal_region_tip_count": str(non_focal_count),
                "manual_tree_review_priority": classify_review_priority(
                    focal_status,
                    dominant_status,
                    source.get("driver_risk", ""),
                ),
                "tree_path": source.get("tree_path", ""),
            }
        )
    write_tsv(output_path, FIELDS, rows)
    return rows


def parse_newick(text: str) -> Node:
    if Phylo is not None:
        tree = Phylo.read(StringIO(text), "newick")
        return from_biopython_clade(tree.root)
    parser = NewickParser(text.rstrip(";"))
    return parser.parse_node()


def newick_parser_backend() -> str:
    return "biopython" if Phylo is not None else "fallback"


def from_biopython_clade(clade) -> Node:
    children = [from_biopython_clade(child) for child in clade.clades]
    return Node(name=clade.name or "", children=children)


class NewickParser:
    def __init__(self, text: str):
        self.text = text
        self.index = 0

    def parse_node(self) -> Node:
        if self.peek() == "(":
            self.index += 1
            children = [self.parse_node()]
            while self.peek() == ",":
                self.index += 1
                children.append(self.parse_node())
            self.expect(")")
            name = self.read_label()
            self.skip_branch_length()
            return Node(name=name, children=children)
        name = self.read_label()
        self.skip_branch_length()
        return Node(name=name)

    def peek(self) -> str:
        self.skip_space()
        if self.index >= len(self.text):
            return ""
        return self.text[self.index]

    def expect(self, value: str) -> None:
        self.skip_space()
        if self.index >= len(self.text) or self.text[self.index] != value:
            raise ValueError(f"Expected {value!r} at Newick offset {self.index}")
        self.index += 1

    def read_label(self) -> str:
        self.skip_space()
        start = self.index
        while self.index < len(self.text) and self.text[self.index] not in ",():":
            self.index += 1
        return self.text[start:self.index].strip()

    def skip_branch_length(self) -> None:
        self.skip_space()
        if self.index >= len(self.text) or self.text[self.index] != ":":
            return
        self.index += 1
        while self.index < len(self.text) and self.text[self.index] not in ",()":
            self.index += 1

    def skip_space(self) -> None:
        while self.index < len(self.text) and self.text[self.index].isspace():
            self.index += 1


def leaf_names(node: Node) -> list[str]:
    if node.is_leaf:
        return [node.name]
    names: list[str] = []
    for child in node.children:
        names.extend(leaf_names(child))
    return names


def species_id(leaf_name: str) -> str:
    return leaf_name.split("|", 1)[0]


def mrca_leaf_names(node: Node, targets: set[str]) -> list[str]:
    best = find_mrca(node, targets)
    return leaf_names(best) if best else []


def find_mrca(node: Node, targets: set[str]) -> Node | None:
    leaves = set(leaf_names(node))
    if not targets.issubset(leaves):
        return None
    for child in node.children:
        child_leaves = set(leaf_names(child))
        if targets.issubset(child_leaves):
            return find_mrca(child, targets)
    return node


def purity(leaves: list[str], target_species: set[str]) -> float:
    if not leaves:
        return 0.0
    target_count = sum(1 for name in leaves if species_id(name) in target_species)
    return target_count / len(leaves)


def classify_region_cluster(
    focal_leaves: list[str],
    mrca_leaves: list[str],
    mrca_purity: float,
    total_tip_count: int,
) -> str:
    if not focal_leaves:
        return "no_focal_region_tips"
    if len(focal_leaves) == len(mrca_leaves) and mrca_purity == 1.0:
        return "region_coherent_clade"
    if len(mrca_leaves) < total_tip_count and mrca_purity >= 0.75:
        return "region_enriched_subclade"
    if len(mrca_leaves) == total_tip_count and mrca_purity >= 0.75:
        return "region_biased_treewide_distribution"
    return "region_intermingled"


def classify_dominant_cluster(dominant_leaves: list[str], mrca_leaves: list[str], mrca_purity: float) -> str:
    if len(dominant_leaves) <= 1:
        return "single_tip_or_absent"
    if len(dominant_leaves) == len(mrca_leaves) and mrca_purity == 1.0:
        return "dominant_species_coherent_clade"
    return "dominant_species_split_or_mixed"


def classify_review_priority(focal_status: str, dominant_status: str, driver_risk: str) -> str:
    if focal_status in {"region_coherent_clade", "region_enriched_subclade"} and driver_risk in {"low", "moderate"}:
        return "inspect_as_best_tree_supported_candidate"
    if dominant_status == "dominant_species_coherent_clade" or driver_risk == "high":
        return "inspect_for_species_specific_duplication"
    return "inspect_for_mixed_or_ancient_family_structure"


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checklist", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = build_candidate_tree_metrics(args.checklist, args.output)
    print(f"Candidate tree metric rows written: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
