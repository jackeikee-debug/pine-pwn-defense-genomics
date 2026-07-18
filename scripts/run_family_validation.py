#!/usr/bin/env python3
"""Run MAFFT and FastTree for regional family review FASTA files."""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
from pathlib import Path


SUMMARY_FIELDS = [
    "validation_rank",
    "orthogroup_rank",
    "orthogroup_id",
    "review_tier",
    "regional_copy_directions",
    "module_ids",
    "matched_keywords",
    "best_abs_log2_ratio",
    "best_log2_ratio",
    "total_selected_genes",
    "review_flags",
    "fasta_path",
    "alignment_path",
    "tree_path",
    "alignment_sequence_count",
    "alignment_length",
    "status",
]


TIER_ORDER = {
    "priority_low_risk": 0,
    "priority_with_outgroup_context": 1,
    "region_specific_candidate": 2,
    "caution_single_species_driver": 3,
    "caution_large_single_species_family": 4,
    "caution_large_family": 5,
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def classify_review_tier(review_flags: str) -> str:
    flags = set(flag for flag in review_flags.split(";") if flag and flag != "none")
    if "large_family" in flags and "single_species_driver" in flags:
        return "caution_large_single_species_family"
    if "single_species_driver" in flags:
        return "caution_single_species_driver"
    if "large_family" in flags:
        return "caution_large_family"
    if "North_America_only" in flags or "East_Asia_only" in flags:
        return "region_specific_candidate"
    if "outgroup_present" in flags:
        return "priority_with_outgroup_context"
    return "priority_low_risk"


def build_validation_rows(
    review_summary_path: Path,
    alignment_dir: Path,
    tree_dir: Path,
    max_orthogroups: int | None = None,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for row in read_tsv(review_summary_path):
        orthogroup_id = row["orthogroup_id"]
        review_tier = classify_review_tier(row.get("review_flags", ""))
        rows.append(
            {
                "validation_rank": "0",
                "orthogroup_rank": row.get("orthogroup_rank", ""),
                "orthogroup_id": orthogroup_id,
                "review_tier": review_tier,
                "regional_copy_directions": row.get("regional_copy_directions", ""),
                "module_ids": row.get("module_ids", ""),
                "matched_keywords": row.get("matched_keywords", ""),
                "best_abs_log2_ratio": row.get("best_abs_log2_ratio", ""),
                "best_log2_ratio": row.get("best_log2_ratio", ""),
                "total_selected_genes": row.get("total_selected_genes", ""),
                "review_flags": row.get("review_flags", ""),
                "fasta_path": row.get("fasta_path", ""),
                "alignment_path": str(alignment_dir / f"{orthogroup_id}.aln.faa"),
                "tree_path": str(tree_dir / f"{orthogroup_id}.fasttree.nwk"),
                "alignment_sequence_count": "0",
                "alignment_length": "0",
                "status": "pending",
            }
        )

    rows.sort(key=validation_sort_key)
    if max_orthogroups is not None:
        rows = rows[:max_orthogroups]
    for index, row in enumerate(rows, start=1):
        row["validation_rank"] = str(index)
    return rows


def validation_sort_key(row: dict[str, str]) -> tuple[int, float, int, str]:
    return (
        TIER_ORDER.get(row["review_tier"], 99),
        -to_float(row.get("best_abs_log2_ratio", "")),
        to_int(row.get("orthogroup_rank", "")),
        row["orthogroup_id"],
    )


def run_family_validation(
    review_summary_path: Path,
    alignment_dir: Path,
    tree_dir: Path,
    summary_output: Path,
    mafft_bin: str = "mafft",
    fasttree_bin: str = "FastTree",
    max_orthogroups: int | None = None,
) -> list[dict[str, str]]:
    alignment_dir.mkdir(parents=True, exist_ok=True)
    tree_dir.mkdir(parents=True, exist_ok=True)
    clean_directory(alignment_dir, "*.aln.faa")
    clean_directory(tree_dir, "*.fasttree.nwk")
    resolved_mafft = resolve_executable(mafft_bin)
    resolved_fasttree = resolve_executable(fasttree_bin)

    rows = build_validation_rows(review_summary_path, alignment_dir, tree_dir, max_orthogroups)
    for row in rows:
        fasta_path = Path(row["fasta_path"])
        alignment_path = Path(row["alignment_path"])
        tree_path = Path(row["tree_path"])
        run_mafft(resolved_mafft, fasta_path, alignment_path)
        run_fasttree(resolved_fasttree, alignment_path, tree_path)
        sequence_count, alignment_length = alignment_stats(alignment_path)
        row["alignment_sequence_count"] = str(sequence_count)
        row["alignment_length"] = str(alignment_length)
        row["status"] = "aligned_and_tree_built"

    write_tsv(summary_output, SUMMARY_FIELDS, rows)
    return rows


def resolve_executable(name: str, search_path: str | None = None) -> str:
    resolved = shutil.which(name, path=search_path)
    if resolved:
        return resolved
    path_ext = os.environ.get("PATHEXT", ".COM;.EXE;.BAT;.CMD").split(";")
    for directory in (search_path or os.environ.get("PATH", "")).split(os.pathsep):
        if not directory:
            continue
        for extension in path_ext:
            candidate = Path(directory) / f"{name}{extension.lower()}"
            if candidate.exists():
                return str(candidate)
            candidate = Path(directory) / f"{name}{extension.upper()}"
            if candidate.exists():
                return str(candidate)
    raise FileNotFoundError(f"Executable not found on PATH: {name}")


def run_mafft(mafft_bin: str, fasta_path: Path, alignment_path: Path) -> None:
    with alignment_path.open("w", encoding="utf-8") as handle:
        subprocess.run(
            [mafft_bin, "--auto", "--thread", "1", str(fasta_path)],
            stdout=handle,
            text=True,
            check=True,
        )


def run_fasttree(fasttree_bin: str, alignment_path: Path, tree_path: Path) -> None:
    with tree_path.open("w", encoding="utf-8") as handle:
        subprocess.run(
            [fasttree_bin, str(alignment_path)],
            stdout=handle,
            text=True,
            check=True,
        )


def alignment_stats(path: Path) -> tuple[int, int]:
    count = 0
    first_length = 0
    current_length = 0
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if count == 1:
                    first_length = current_length
                count += 1
                current_length = 0
            else:
                current_length += len(line)
    if count == 1:
        first_length = current_length
    elif count > 1 and first_length == 0:
        first_length = current_length
    return count, first_length


def clean_directory(path: Path, pattern: str) -> None:
    for item in path.glob(pattern):
        item.unlink()


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def to_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def to_int(value: str) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 10**9


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--review-summary", type=Path, required=True)
    parser.add_argument("--alignment-dir", type=Path, required=True)
    parser.add_argument("--tree-dir", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--mafft-bin", default="mafft")
    parser.add_argument("--fasttree-bin", default="FastTree")
    parser.add_argument("--max-orthogroups", type=int)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = run_family_validation(
        review_summary_path=args.review_summary,
        alignment_dir=args.alignment_dir,
        tree_dir=args.tree_dir,
        summary_output=args.summary,
        mafft_bin=args.mafft_bin,
        fasttree_bin=args.fasttree_bin,
        max_orthogroups=args.max_orthogroups,
    )
    print(f"Regional family alignments and trees written: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
