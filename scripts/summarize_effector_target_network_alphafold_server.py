#!/usr/bin/env python3
"""Summarize AlphaFold Server outputs for high-priority Route 1 candidates."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import re
import statistics
from collections import defaultdict
from pathlib import Path


CATALYTIC_PATTERN = re.compile(r"active site|catalytic domain|catalytic chain", re.I)
ACTIVE_SITE_PATTERN = re.compile(r"active site", re.I)
PROTEIN_ID_PATTERN = re.compile(r"cag(\d+)_1$", re.I)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-zip-dir", type=Path, required=True)
    parser.add_argument("--extracted-dir", type=Path, required=True)
    parser.add_argument("--interproscan", type=Path, required=True)
    parser.add_argument("--archive-manifest", type=Path, required=True)
    parser.add_argument("--model-table", type=Path, required=True)
    parser.add_argument("--best-model-table", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    return parser.parse_args()


def long_path(path: Path) -> str:
    resolved = str(path.resolve())
    if resolved.startswith("\\\\?\\"):
        return resolved
    return "\\\\?\\" + resolved


def read_json(path: Path):
    with open(long_path(path), encoding="utf-8") as handle:
        return json.load(handle)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(long_path(path), "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def protein_id_from_directory(path: Path) -> str:
    match = PROTEIN_ID_PATTERN.search(path.name)
    if not match:
        raise ValueError(f"Cannot derive protein ID from {path.name}")
    return f"CAG{match.group(1)}.1"


def parse_interpro_ranges(path: Path) -> tuple[dict[str, list[tuple[int, int]]], dict[str, list[tuple[int, int]]]]:
    catalytic: dict[str, list[tuple[int, int]]] = defaultdict(list)
    active_site: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            row.extend(["-"] * (15 - len(row)))
            description = " ".join((row[5], row[12]))
            if CATALYTIC_PATTERN.search(description):
                interval = (int(row[6]), int(row[7]))
                catalytic[row[0]].append(interval)
                if ACTIVE_SITE_PATTERN.search(description):
                    active_site[row[0]].append(interval)
    return catalytic, active_site


def parse_ca_plddt(path: Path) -> dict[int, float]:
    scores: dict[int, float] = {}
    with open(long_path(path), encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("ATOM "):
                continue
            fields = line.split()
            if len(fields) < 15 or fields[3] != "CA":
                continue
            scores[int(fields[8])] = float(fields[14])
    if not scores:
        raise ValueError(f"No CA pLDDT values found in {path}")
    return scores


def interval_scores(scores: dict[int, float], intervals: list[tuple[int, int]]) -> list[float]:
    residues = {position for start, end in intervals for position in range(start, end + 1)}
    return [score for position, score in scores.items() if position in residues]


def format_float(value: float | None, digits: int = 3) -> str:
    return "" if value is None else f"{value:.{digits}f}"


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def quality_tier(mean_plddt: float, ptm: float) -> str:
    if mean_plddt >= 90 and ptm >= 0.80:
        return "high_global_confidence"
    if mean_plddt >= 70 and ptm >= 0.50:
        return "moderate_global_confidence"
    return "low_or_mixed_global_confidence"


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    catalytic_ranges, active_site_ranges = parse_interpro_ranges(args.interproscan)

    archive_rows: list[dict[str, str]] = []
    for archive in sorted(args.raw_zip_dir.glob("*.zip")):
        archive_rows.append(
            {
                "archive_name": archive.name,
                "archive_bytes": str(archive.stat().st_size),
                "sha256": sha256(archive),
                "source": "AlphaFold_Server_manual_download",
                "download_date": "2026-07-13",
            }
        )

    model_rows: list[dict[str, str]] = []
    for directory in sorted(path for path in args.extracted_dir.glob("fold_*") if path.is_dir()):
        protein_id = protein_id_from_directory(directory)
        job_path = next(directory.glob("*job_request.json"))
        job = read_json(job_path)[0]
        sequence = job["sequences"][0]["proteinChain"]["sequence"]
        seed = str(job.get("modelSeeds", [""])[0])
        summaries = sorted(directory.glob("*summary_confidences_*.json"))
        if len(summaries) != 5:
            raise ValueError(f"Expected five models for {protein_id}, found {len(summaries)}")
        for summary_path in summaries:
            model_match = re.search(r"_(\d+)\.json$", summary_path.name)
            if not model_match:
                raise ValueError(f"Cannot parse model index from {summary_path.name}")
            model_index = int(model_match.group(1))
            cif_path = next(directory.glob(f"*model_{model_index}.cif"))
            summary = read_json(summary_path)
            scores = parse_ca_plddt(cif_path)
            all_scores = list(scores.values())
            catalytic_scores = interval_scores(scores, catalytic_ranges.get(protein_id, []))
            active_scores = interval_scores(scores, active_site_ranges.get(protein_id, []))
            mean_plddt = statistics.fmean(all_scores)
            ptm = float(summary["ptm"])
            model_rows.append(
                {
                    "protein_id": protein_id,
                    "job_name": job["name"],
                    "model_seed": seed,
                    "model_index": str(model_index),
                    "sequence_length": str(len(sequence)),
                    "ranking_score": format_float(float(summary["ranking_score"])),
                    "ptm": format_float(ptm),
                    "mean_ca_plddt": format_float(mean_plddt, 2),
                    "median_ca_plddt": format_float(statistics.median(all_scores), 2),
                    "fraction_ca_plddt_ge_90": format_float(sum(score >= 90 for score in all_scores) / len(all_scores)),
                    "fraction_ca_plddt_ge_70": format_float(sum(score >= 70 for score in all_scores) / len(all_scores)),
                    "fraction_ca_plddt_lt_50": format_float(sum(score < 50 for score in all_scores) / len(all_scores)),
                    "fraction_disordered": format_float(float(summary["fraction_disordered"])),
                    "has_clash": str(bool(summary["has_clash"])).lower(),
                    "num_recycles": str(int(summary["num_recycles"])),
                    "catalytic_feature_residue_count": str(len(catalytic_scores)),
                    "catalytic_feature_mean_plddt": format_float(statistics.fmean(catalytic_scores) if catalytic_scores else None, 2),
                    "active_site_residue_count": str(len(active_scores)),
                    "active_site_mean_plddt": format_float(statistics.fmean(active_scores) if active_scores else None, 2),
                    "global_quality_tier": quality_tier(mean_plddt, ptm),
                    "cif_path": cif_path.as_posix(),
                    "summary_json_path": summary_path.as_posix(),
                }
            )

    best_rows: list[dict[str, str]] = []
    by_protein: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in model_rows:
        by_protein[row["protein_id"]].append(row)
    for protein_id, rows in sorted(by_protein.items()):
        best = max(rows, key=lambda row: float(row["ranking_score"]))
        best_rows.append(
            {
                "protein_id": protein_id,
                "job_name": best["job_name"],
                "model_seed": best["model_seed"],
                "best_model_index": best["model_index"],
                "sequence_length": best["sequence_length"],
                "ranking_score": best["ranking_score"],
                "ptm": best["ptm"],
                "mean_ca_plddt": best["mean_ca_plddt"],
                "median_ca_plddt": best["median_ca_plddt"],
                "fraction_ca_plddt_ge_90": best["fraction_ca_plddt_ge_90"],
                "fraction_ca_plddt_ge_70": best["fraction_ca_plddt_ge_70"],
                "fraction_ca_plddt_lt_50": best["fraction_ca_plddt_lt_50"],
                "catalytic_feature_mean_plddt": best["catalytic_feature_mean_plddt"],
                "active_site_mean_plddt": best["active_site_mean_plddt"],
                "global_quality_tier": best["global_quality_tier"],
                "best_cif_path": best["cif_path"],
                "claim_ceiling": "AlphaFold 3 Server structure prediction only; does not validate secretion, enzymatic activity, host binding, or causal effector function",
            }
        )

    if len(archive_rows) != 8 or len(model_rows) != 40 or len(best_rows) != 8:
        raise ValueError(
            f"Unexpected counts: archives={len(archive_rows)}, models={len(model_rows)}, proteins={len(best_rows)}"
        )
    write_tsv(args.archive_manifest, archive_rows)
    write_tsv(args.model_table, model_rows)
    write_tsv(args.best_model_table, best_rows)

    tier_counts: dict[str, int] = defaultdict(int)
    for row in best_rows:
        tier_counts[row["global_quality_tier"]] += 1
    report = [
        "# Route 1 AlphaFold Server Structure Audit",
        "",
        "## Provenance",
        "",
        "- Eight independent monomer jobs were run with AlphaFold Server (AlphaFold 3) and downloaded on 2026-07-13.",
        "- Each archive contains five predicted models, confidence JSON files, mmCIF coordinates, MSA/template evidence, the submitted job request, and the server output terms.",
        "- Archive SHA-256 checksums are recorded in the companion manifest.",
        "",
        "## Best-model summary",
        "",
    ]
    report.extend(f"- {tier}: {count}" for tier, count in sorted(tier_counts.items()))
    report.extend(
        [
            "",
            "The best model for each candidate is selected by AlphaFold Server ranking score. Global CA pLDDT and pTM summarize overall confidence; InterProScan active-site or catalytic-feature intervals are summarized separately when present.",
            "",
            "## Interpretation boundary",
            "",
            "Predicted structural confidence does not demonstrate secretion, enzymatic activity, physical binding to a pine protein, or causal contribution to pine wilt disease. AlphaFold Server outputs remain subject to its output terms of use.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(report) + "\n", encoding="utf-8")
    logging.info("Wrote %d archive, %d model, and %d best-model records", len(archive_rows), len(model_rows), len(best_rows))


if __name__ == "__main__":
    main()
