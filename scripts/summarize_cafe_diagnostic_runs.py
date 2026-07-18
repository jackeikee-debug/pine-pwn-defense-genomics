#!/usr/bin/env python3
"""Summarize CAFE5 diagnostic models, failures, boundaries, and convergence."""

from __future__ import annotations

import argparse
import csv
import itertools
import logging
import math
import re
import statistics
from collections import defaultdict
from pathlib import Path


RUN_FIELDS = [
    "run_id", "replicate_group", "panel_id", "panel_scope", "filter_level", "model",
    "use_error_model", "n_gamma_categories", "execution_status", "final_negative_log_likelihood",
    "lambda", "maximum_possible_lambda", "lambda_to_max_ratio", "boundary_status", "alpha",
    "epsilon", "gamma_multipliers", "attempted_values", "rejected_percent", "iterations",
    "input_families", "root_retained_families", "significant_families",
    "significant_family_fraction", "families_with_high_failure_rate",
    "high_failure_family_fraction", "maximum_family_failure_fraction",
    "numerical_quality_status", "elapsed_seconds", "result_file", "family_result_file",
    "log_path", "claim_ceiling",
]
CONVERGENCE_FIELDS = [
    "replicate_group", "panel_id", "model", "replicate_count", "successful_replicates",
    "negative_log_likelihood_range", "lambda_mean", "lambda_cv",
    "minimum_significant_family_jaccard", "all_interior", "all_numerically_clean",
    "convergence_status",
]
CLAIM = (
    "CAFE5 diagnostic sensitivity evidence only; provisional time tree and heterogeneous proteome annotations "
    "preclude final expansion/contraction, regional resistance, coevolution, or effector-target claims"
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def first_match(pattern: str, text: str) -> str:
    match = re.search(pattern, text, flags=re.IGNORECASE)
    return match.group(1) if match else ""


def finite_float(value: str) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def resolve(path: str) -> Path:
    return Path(path)


def read_text_auto(path: Path) -> str:
    data = path.read_bytes()
    if data.startswith((b"\xff\xfe", b"\xfe\xff")):
        return data.decode("utf-16", errors="replace")
    return data.decode("utf-8", errors="replace")


def significant_family_set(path: Path) -> set[str]:
    if not path.exists():
        return set()
    lines = path.read_text(encoding="utf-8").splitlines()
    if lines and lines[0].startswith("#"):
        lines[0] = lines[0][1:]
    rows = csv.DictReader(lines, delimiter="\t")
    return {
        row.get("FamilyID", "") for row in rows
        if row.get("Significant at 0.05", "").strip().lower() == "y" and row.get("FamilyID", "")
    }


def parse_run(run: dict[str, str]) -> dict[str, str]:
    model = run["model"]
    run_dir = resolve(run["run_dir"])
    result_file = run_dir / f"{model}_results.txt"
    family_file = run_dir / f"{model}_family_results.txt"
    log_path = resolve(run["log_path"])
    base = {field: run.get(field, "") for field in RUN_FIELDS}
    base.update({
        "execution_status": "missing_result_file", "boundary_status": "unavailable",
        "numerical_quality_status": "unavailable",
        "result_file": str(result_file), "family_result_file": str(family_file),
        "log_path": str(log_path), "claim_ceiling": CLAIM,
    })
    if not result_file.exists():
        return base

    text = read_text_auto(result_file)
    log_text = read_text_auto(log_path) if log_path.exists() else ""
    nll = first_match(r"Final Likelihood \(-lnL\):\s*(\S+)", text)
    rate = first_match(r"\bLambda(?:s)?\s*:\s*([0-9.eE+\-]+)", text)
    maximum = first_match(r"Maximum possible lambda for this topology:\s*(\S+)", text)
    alpha = first_match(r"\bAlpha\s*:\s*(\S+)", text)
    epsilon = first_match(r"\b(?:Epsilon|Error)\s*:\s*(\S+)", text)
    multipliers_text = first_match(r"Gamma category multipliers\s*:\s*([^\r\n]+)", text)
    attempted = re.search(r"(\d+) values were attempted \(([^%]+)% rejected\)", text)
    family_failures = [
        int(value) for value in re.findall(r"^\S+ had (\d+) failures\r?$", text, flags=re.MULTILINE)
    ]
    nll_number = finite_float(nll)
    rate_number = finite_float(rate)
    maximum_number = finite_float(maximum)
    ratio = rate_number / maximum_number if rate_number is not None and maximum_number not in (None, 0.0) else None
    status = "ok"
    if nll_number is None:
        status = "nonfinite_or_missing_likelihood"
    elif rate_number is None or maximum_number is None:
        status = "incomplete_model_parameters"
    significant = significant_family_set(family_file)
    input_families = ""
    summary_path = resolve(run.get("input_summary", "")) if run.get("input_summary") else None
    if summary_path and summary_path.exists():
        summaries = read_tsv(summary_path)
        input_families = summaries[0].get("written_families", "") if summaries else ""
    root_match = re.search(r"Filtering families not present at the root from:\s*(\d+)\s+to\s+(\d+)", log_text)
    root_retained = root_match.group(2) if root_match else ""
    denominator = int(root_retained) if root_retained.isdigit() and int(root_retained) else 0
    attempted_count = int(attempted.group(1)) if attempted else 0
    time_match = re.search(r"Time:\s*(\d+)H\s*(\d+)M\s*(\d+)S", log_text)
    elapsed = str(int(time_match.group(1)) * 3600 + int(time_match.group(2)) * 60 + int(time_match.group(3))) if time_match else ""
    iterations = first_match(r"Completed\s+(\d+)\s+iterations", log_text)
    multipliers = ";".join(re.findall(r"[0-9.eE+\-]+", multipliers_text)) if multipliers_text else ""
    base.update({
        "execution_status": status, "final_negative_log_likelihood": nll,
        "lambda": rate, "maximum_possible_lambda": maximum,
        "lambda_to_max_ratio": f"{ratio:.6f}" if ratio is not None else "",
        "boundary_status": "interior" if ratio is not None and ratio < 0.95 else ("boundary_or_near_boundary" if ratio is not None else "unavailable"),
        "alpha": alpha, "epsilon": epsilon, "gamma_multipliers": multipliers,
        "attempted_values": attempted.group(1) if attempted else "",
        "rejected_percent": attempted.group(2) if attempted else "", "iterations": iterations,
        "input_families": input_families, "root_retained_families": root_retained,
        "significant_families": str(len(significant)),
        "significant_family_fraction": f"{len(significant) / denominator:.6f}" if denominator else "",
        "families_with_high_failure_rate": str(len(family_failures)),
        "high_failure_family_fraction": f"{len(family_failures) / denominator:.6f}" if denominator else "",
        "maximum_family_failure_fraction": (
            f"{max(family_failures) / attempted_count:.6f}" if family_failures and attempted_count else ""
        ),
        "numerical_quality_status": "high_family_failure_warning" if family_failures else "clean",
        "elapsed_seconds": elapsed,
    })
    return base


def jaccard(left: set[str], right: set[str]) -> float:
    union = left | right
    return len(left & right) / len(union) if union else 1.0


def summarize_convergence(
    rows: list[dict[str, str]],
    significant_sets: dict[str, set[str]],
) -> list[dict[str, str]]:
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[row["replicate_group"]].append(row)
    output: list[dict[str, str]] = []
    for group_id, members in sorted(groups.items()):
        successful = [row for row in members if row["execution_status"] == "ok"]
        nlls = [float(row["final_negative_log_likelihood"]) for row in successful]
        rates = [float(row["lambda"]) for row in successful]
        comparisons = [
            jaccard(significant_sets.get(left["run_id"], set()), significant_sets.get(right["run_id"], set()))
            for left, right in itertools.combinations(successful, 2)
        ]
        nll_range = max(nlls) - min(nlls) if nlls else math.nan
        rate_mean = statistics.mean(rates) if rates else math.nan
        rate_cv = statistics.pstdev(rates) / rate_mean if len(rates) > 1 and rate_mean else (0.0 if len(rates) == 1 else math.nan)
        minimum_jaccard = min(comparisons) if comparisons else (1.0 if len(successful) == 1 else math.nan)
        all_interior = bool(successful) and all(row["boundary_status"] == "interior" for row in successful)
        all_numerically_clean = bool(successful) and all(
            row.get("numerical_quality_status") == "clean" for row in successful
        )
        if len(members) < 3:
            convergence = "not_replicated"
        elif len(successful) != len(members):
            convergence = "failed_replicate"
        elif not all_numerically_clean:
            convergence = "numerical_warning"
        elif nll_range <= 1.0 and rate_cv <= 0.05 and all_interior and minimum_jaccard >= 0.95:
            convergence = "converged"
        else:
            convergence = "not_converged"
        output.append({
            "replicate_group": group_id, "panel_id": members[0]["panel_id"],
            "model": members[0]["model"], "replicate_count": str(len(members)),
            "successful_replicates": str(len(successful)),
            "negative_log_likelihood_range": f"{nll_range:.6f}" if math.isfinite(nll_range) else "",
            "lambda_mean": f"{rate_mean:.12g}" if math.isfinite(rate_mean) else "",
            "lambda_cv": f"{rate_cv:.6f}" if math.isfinite(rate_cv) else "",
            "minimum_significant_family_jaccard": f"{minimum_jaccard:.6f}" if math.isfinite(minimum_jaccard) else "",
            "all_interior": str(all_interior).lower(),
            "all_numerically_clean": str(all_numerically_clean).lower(),
            "convergence_status": convergence,
        })
    return output


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def summarize_runs(
    manifest_path: Path,
    output_path: Path,
    convergence_path: Path,
    report_path: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    definitions = read_tsv(manifest_path)
    rows = [parse_run(run) for run in definitions]
    definitions_by_id = {row["run_id"]: row for row in definitions}
    sets = {
        row["run_id"]: significant_family_set(Path(row["family_result_file"]))
        for row in rows if row["run_id"] in definitions_by_id
    }
    convergence = summarize_convergence(rows, sets)
    write_tsv(output_path, RUN_FIELDS, rows)
    write_tsv(convergence_path, CONVERGENCE_FIELDS, convergence)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    ok = sum(row["execution_status"] == "ok" for row in rows)
    interior = sum(row["boundary_status"] == "interior" for row in rows)
    report_path.write_text(
        "# CAFE5 diagnostic model comparison\n\n"
        f"- Predeclared fitted runs: {len(rows)}\n"
        f"- Successfully parsed runs: {ok}\n"
        f"- Interior-rate runs: {interior}\n\n"
        "This table is a convergence and sensitivity audit using a provisional time tree. It does not establish "
        "final gene-family expansion/contraction or geographic resistance differences.\n",
        encoding="utf-8",
    )
    logging.info("Parsed %d diagnostic runs", len(rows))
    return rows, convergence


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--convergence", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    summarize_runs(args.manifest, args.output, args.convergence, args.report)


if __name__ == "__main__":
    main()
