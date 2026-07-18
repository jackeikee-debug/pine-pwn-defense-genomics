#!/usr/bin/env python3
"""Assess predeclared mechanism candidates across accepted CAFE5 diagnostics."""

from __future__ import annotations

import argparse
import csv
import logging
from collections import defaultdict
from pathlib import Path


FIELDS = [
    "orthogroup_id", "mechanism_axes", "regional_direction", "driver_risk",
    "accepted_host_models", "accepted_conifer_models", "accepted_model_count",
    "significant_accepted_models", "significant_model_fraction", "cross_panel_support",
    "single_species_driver_flag", "stability_status", "claim_ceiling",
]
CLAIM = (
    "CAFE5 sensitivity support only; does not establish geographic resistance, "
    "coevolution, or direct effector targeting"
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def significant_families(path: Path) -> set[str]:
    if not path.exists():
        return set()
    lines = path.read_text(encoding="utf-8").splitlines()
    if lines and lines[0].startswith("#"):
        lines[0] = lines[0][1:]
    return {
        row.get("FamilyID", "")
        for row in csv.DictReader(lines, delimiter="\t")
        if row.get("FamilyID") and row.get("Significant at 0.05", "").lower() == "y"
    }


def candidate_axes(rows: list[dict[str, str]]) -> dict[str, set[str]]:
    axes: dict[str, set[str]] = defaultdict(set)
    columns = (
        "primary_candidate_orthogroups",
        "supporting_candidate_orthogroups",
        "deprioritized_candidate_orthogroups",
    )
    for row in rows:
        for column in columns:
            for orthogroup in row.get(column, "").split(";"):
                if orthogroup:
                    axes[orthogroup].add(row.get("mechanism_axis", ""))
    return axes


def build_candidate_stability(
    model_table: Path,
    convergence_table: Path,
    run_manifest: Path,
    candidates: Path,
    review_decisions: Path,
    output_path: Path,
    report_path: Path,
) -> tuple[list[dict[str, str]], str]:
    models = read_tsv(model_table)
    converged_groups = {
        row["replicate_group"] for row in read_tsv(convergence_table)
        if row.get("convergence_status") == "converged"
    }
    accepted = {
        row["run_id"]: row for row in models
        if row.get("replicate_group") in converged_groups
        and row.get("execution_status") == "ok"
        and row.get("boundary_status") == "interior"
        and row.get("numerical_quality_status") == "clean"
    }
    manifests = {row["run_id"]: row for row in read_tsv(run_manifest)}
    significant = {}
    for run_id in accepted:
        run = manifests.get(run_id, {})
        path = Path(run.get("run_dir", "")) / f"{run.get('model', '')}_family_results.txt"
        significant[run_id] = significant_families(path)

    axes = candidate_axes(read_tsv(candidates))
    reviews = {row["orthogroup_id"]: row for row in read_tsv(review_decisions)}
    host_models = {run_id for run_id, row in accepted.items() if row.get("panel_scope") == "pinus_host_focused"}
    conifer_models = {run_id for run_id, row in accepted.items() if row.get("panel_scope") == "conifer_wide"}
    output: list[dict[str, str]] = []
    for orthogroup in sorted(axes):
        review = reviews.get(orthogroup, {})
        supported = {run_id for run_id, families in significant.items() if orthogroup in families}
        host_supported = bool(supported & host_models)
        conifer_supported = bool(supported & conifer_models)
        high_driver = review.get("driver_risk", "").lower() == "high"
        if not accepted:
            status = "no_accepted_models"
        elif high_driver:
            status = "blocked_by_single_species_driver"
        elif host_supported and conifer_supported and review.get("regional_direction"):
            status = "stable_cross_panel_candidate"
        elif supported:
            status = "panel_specific_or_unstable"
        else:
            status = "not_supported_by_accepted_models"
        output.append({
            "orthogroup_id": orthogroup,
            "mechanism_axes": ";".join(sorted(axis for axis in axes[orthogroup] if axis)),
            "regional_direction": review.get("regional_direction", ""),
            "driver_risk": review.get("driver_risk", ""),
            "accepted_host_models": str(len(host_models)),
            "accepted_conifer_models": str(len(conifer_models)),
            "accepted_model_count": str(len(accepted)),
            "significant_accepted_models": str(len(supported)),
            "significant_model_fraction": f"{len(supported) / len(accepted):.6f}" if accepted else "",
            "cross_panel_support": str(host_supported and conifer_supported).lower(),
            "single_species_driver_flag": str(high_driver).lower(),
            "stability_status": status,
            "claim_ceiling": CLAIM,
        })

    stable = sum(row["stability_status"] == "stable_cross_panel_candidate" for row in output)
    boundary_models = sum(row.get("boundary_status") != "interior" for row in models)
    warned_interior = [
        row for row in models
        if row.get("boundary_status") == "interior"
        and row.get("numerical_quality_status") != "clean"
    ]
    warning_lines = "".join(
        f"- `{row['run_id']}`: {row.get('families_with_high_failure_rate', 'NA')} families "
        f"({row.get('high_failure_family_fraction', 'NA')}) exceeded the CAFE5 warning threshold.\n"
        for row in warned_interior
    )
    decision = (
        "proceed_to_formal_cafe"
        if host_models and conifer_models and stable
        else "retain_cafe_as_supplementary"
    )
    write_tsv(output_path, output)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        "# CAFE5 strengthening diagnostic\n\n"
        "## Decision\n\n"
        f"- Decision: `{decision}`\n"
        f"- First-pass models parsed: {len(models)}\n"
        f"- Boundary or near-boundary models: {boundary_models}\n"
        f"- Interior but numerically warned models: {len(warned_interior)}\n"
        f"- Accepted host-focused models: {len(host_models)}\n"
        f"- Accepted conifer-wide models: {len(conifer_models)}\n"
        f"- Stable cross-panel candidates: {stable}\n\n"
        "## Diagnostic outcome\n\n"
        + (warning_lines if warning_lines else "- No interior model carried a high-family-failure warning.\n")
        + "\n"
        "Base and error-aware Base fits that remain at the topology-specific lambda ceiling are not "
        "treated as final expansion/contraction models. Gamma fits may move lambda into the interior, "
        "but models with widespread family-level calculation warnings are not promoted or replicated.\n\n"
        "## Interpretation ceiling\n\n"
        "Acceptance requires successful execution, an interior lambda, clean numerical quality, "
        "and replicated convergence. Lineage-level CAFE results are not interpreted as geographic "
        "resistance or direct effector-target evidence.\n",
        encoding="utf-8",
    )
    logging.info("Assessed %d candidates; decision=%s", len(output), decision)
    return output, decision


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-table", type=Path, required=True)
    parser.add_argument("--convergence-table", type=Path, required=True)
    parser.add_argument("--run-manifest", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--review-decisions", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    build_candidate_stability(
        model_table=args.model_table,
        convergence_table=args.convergence_table,
        run_manifest=args.run_manifest,
        candidates=args.candidates,
        review_decisions=args.review_decisions,
        output_path=args.output,
        report_path=args.report,
    )


if __name__ == "__main__":
    main()
