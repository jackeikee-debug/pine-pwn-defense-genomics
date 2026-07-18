#!/usr/bin/env python3
"""Write a readable tree-review report for model-supported candidate families."""

from __future__ import annotations

import argparse
import csv
from io import StringIO
from pathlib import Path

try:
    from Bio import Phylo
except ImportError:  # pragma: no cover - fallback supports minimal environments.
    Phylo = None


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_candidate_tree_review_report(
    checklist_path: Path,
    metrics_path: Path,
    markdown_path: Path,
    ascii_dir: Path,
) -> list[str]:
    checklist_rows = read_tsv(checklist_path)
    metrics_by_og = {row["orthogroup_id"]: row for row in read_tsv(metrics_path)}
    ascii_dir.mkdir(parents=True, exist_ok=True)
    for stale in ascii_dir.glob("*.ascii.txt"):
        stale.unlink()
    markdown_path.parent.mkdir(parents=True, exist_ok=True)

    written_ids: list[str] = []
    lines = [
        "# Model-Supported Candidate Tree Review",
        "",
        "This report is a manual-review aid for model-supported candidate families.",
        "ASCII trees are topology snapshots, not final phylogenetic evidence.",
        "",
    ]
    for row in checklist_rows:
        orthogroup_id = row["orthogroup_id"]
        metric = metrics_by_og.get(orthogroup_id, {})
        ascii_path = ascii_dir / f"{orthogroup_id}.ascii.txt"
        tree_text = ascii_tree_text(Path(row["tree_path"]))
        ascii_path.write_text(tree_text, encoding="utf-8")
        written_ids.append(orthogroup_id)
        lines.extend(candidate_markdown_block(row, metric, ascii_path))

    markdown_path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    return written_ids


def ascii_tree_text(tree_path: Path) -> str:
    if Phylo is not None:
        tree = Phylo.read(str(tree_path), "newick")
        handle = StringIO()
        Phylo.draw_ascii(tree, file=handle)
        return handle.getvalue()
    return fallback_tip_list(tree_path)


def fallback_tip_list(tree_path: Path) -> str:
    text = tree_path.read_text(encoding="utf-8")
    tips = []
    token = []
    for char in text:
        if char in "(),:;":
            if token:
                label = "".join(token).strip()
                if "|" in label:
                    tips.append(label)
                token = []
        else:
            token.append(char)
    if token:
        label = "".join(token).strip()
        if "|" in label:
            tips.append(label)
    return "\n".join(tips) + "\n"


def candidate_markdown_block(row: dict[str, str], metric: dict[str, str], ascii_path: Path) -> list[str]:
    return [
        f"## {row['check_rank']}. {row['orthogroup_id']} ({row['evidence_tier']})",
        "",
        f"- Direction: {row['regional_direction']}",
        f"- Modules: {row['module_ids']}",
        f"- Effector classes: {row['effector_classes']}",
        f"- Matched keywords: {row['matched_keywords']}",
        f"- Driver risk: {row['driver_risk']} ({row['focal_max_species']} share {row['focal_max_share']})",
        f"- Tree status: {metric.get('focal_region_cluster_status', '')}",
        f"- Focal MRCA purity: {metric.get('focal_region_mrca_purity', '')}",
        f"- Tree review priority: {metric.get('manual_tree_review_priority', '')}",
        f"- Review question: {row['review_question']}",
        f"- ASCII tree: `{ascii_path}`",
        "",
    ]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checklist", type=Path, required=True)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    parser.add_argument("--ascii-dir", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    ids = build_candidate_tree_review_report(
        checklist_path=args.checklist,
        metrics_path=args.metrics,
        markdown_path=args.markdown,
        ascii_dir=args.ascii_dir,
    )
    print(f"Candidate tree review entries written: {len(ids)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
