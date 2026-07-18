#!/usr/bin/env python3
"""Build a candidate-restricted CAFE5 turnover sensitivity table."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path


INTERPRETATION_STATUS = "sensitivity_only_model_not_accepted"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def safe_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def expected_direction(change: float) -> str:
    if change > 0:
        return "expansion"
    if change < 0:
        return "contraction"
    return "no_change"


def candidate_label(row: dict[str, str]) -> str:
    return row.get("candidate_label") or row.get("mapped_symbols") or ""


def diagnostic_summary(rows: list[dict[str, str]]) -> str:
    if not rows:
        return "diagnostic_models_unavailable"
    boundary = sorted({row.get("boundary_status", "unknown") or "unknown" for row in rows})
    quality = sorted(
        {row.get("numerical_quality_status", "unknown") or "unknown" for row in rows}
    )
    return f"boundary={','.join(boundary)};numerical_quality={','.join(quality)}"


def accepted_model_ids(
    diagnostics: list[dict[str, str]],
    convergence: list[dict[str, str]],
) -> set[str]:
    converged_groups = {
        row.get("replicate_group", "")
        for row in convergence
        if row.get("convergence_status") == "converged"
    }
    return {
        row.get("run_id", "")
        for row in diagnostics
        if row.get("replicate_group", "") in converged_groups
        and row.get("execution_status") == "ok"
        and row.get("boundary_status") == "interior"
        and row.get("numerical_quality_status") == "clean"
    }


def build_turnover_rows(
    candidates: list[dict[str, str]],
    families: list[dict[str, str]],
    branches: list[dict[str, str]],
    model_summary: list[dict[str, str]],
    diagnostics: list[dict[str, str]],
    convergence: list[dict[str, str]],
) -> list[dict[str, object]]:
    candidate_ids = [row.get("orthogroup_id", "") for row in candidates]
    if len(candidate_ids) != 12 or len(set(candidate_ids)) != 12 or any(not value for value in candidate_ids):
        raise ValueError("Expected exactly 12 unique candidate orthogroups")
    if not diagnostics or not convergence:
        raise ValueError("Diagnostic model and convergence tables must be nonempty")
    diagnostic_groups = {row.get("replicate_group", "") for row in diagnostics}
    convergence_groups = {row.get("replicate_group", "") for row in convergence}
    if "" in diagnostic_groups or not diagnostic_groups.issubset(convergence_groups):
        missing = sorted(diagnostic_groups - convergence_groups)
        raise ValueError(
            "Every diagnostic replicate group must have convergence evidence; "
            f"missing={','.join(missing) if missing else 'blank_replicate_group'}"
        )
    accepted = accepted_model_ids(diagnostics, convergence)
    if accepted:
        raise ValueError(
            "At least one accepted diagnostic model is present; rejected-model sensitivity labeling "
            f"is invalid: {','.join(sorted(accepted))}"
        )

    family_by_id = {row.get("family_id", ""): row for row in families}
    branches_by_family: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in branches:
        branches_by_family[row.get("family_id", "")].append(row)

    model = model_summary[0] if model_summary else {}
    warning = diagnostic_summary(diagnostics)
    output: list[dict[str, object]] = []

    for order, candidate in enumerate(candidates, 1):
        orthogroup = candidate["orthogroup_id"]
        family = family_by_id.get(orthogroup)
        family_branches = sorted(
            branches_by_family.get(orthogroup, []),
            key=lambda row: (row.get("taxon_label", ""), row.get("branch_id", "")),
        )
        common = {
            "candidate_order": order,
            "orthogroup_id": orthogroup,
            "candidate_label": candidate_label(candidate),
            "model": model.get("model", "not_available"),
            "lambda": model.get("lambda", ""),
            "maximum_possible_lambda": model.get("maximum_possible_lambda", ""),
            "model_diagnostic_summary": warning,
            "interpretation_status": INTERPRETATION_STATUS,
        }
        if family is None or not family_branches:
            output.append(
                {
                    **common,
                    "availability_state": "not_available",
                    "family_pvalue": "",
                    "branch_id": "not_available",
                    "taxon_label": "not_available",
                    "change": "",
                    "direction": "not_available",
                    "branch_pvalue": "",
                }
            )
            continue

        for branch in family_branches:
            change = safe_float(branch.get("change"))
            if not math.isfinite(change):
                raise ValueError(f"Nonfinite change for {orthogroup} {branch.get('branch_id', '')}")
            direction = branch.get("direction", "")
            if direction != expected_direction(change):
                raise ValueError(
                    f"Branch direction does not match change for {orthogroup} "
                    f"{branch.get('branch_id', '')}: {direction} versus {change}"
                )
            output.append(
                {
                    **common,
                    "availability_state": "available",
                    "family_pvalue": family.get("family_pvalue", ""),
                    "branch_id": branch.get("branch_id", ""),
                    "taxon_label": branch.get("taxon_label", ""),
                    "change": branch.get("change", ""),
                    "direction": direction,
                    "branch_pvalue": branch.get("branch_pvalue", ""),
                }
            )
    return output


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError("No turnover rows were generated")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--family-results", required=True, type=Path)
    parser.add_argument("--branch-results", required=True, type=Path)
    parser.add_argument("--model-summary", required=True, type=Path)
    parser.add_argument("--diagnostics", required=True, type=Path)
    parser.add_argument("--convergence", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rows = build_turnover_rows(
        read_tsv(args.candidates),
        read_tsv(args.family_results),
        read_tsv(args.branch_results),
        read_tsv(args.model_summary),
        read_tsv(args.diagnostics),
        read_tsv(args.convergence),
    )
    write_tsv(args.output, rows)
    available = len({row["orthogroup_id"] for row in rows if row["availability_state"] == "available"})
    print(
        f"Wrote {len(rows)} turnover rows for {len({row['orthogroup_id'] for row in rows})} "
        f"candidates ({available} available) to {args.output}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
