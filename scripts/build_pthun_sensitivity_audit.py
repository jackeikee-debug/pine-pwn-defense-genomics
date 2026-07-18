#!/usr/bin/env python3
"""Compare full and low-mapping-sample-excluded P. thunbergii evidence."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


FIELDS = [
    "orthogroup_id", "contrast", "mechanism_axes", "full_classification",
    "sensitivity_classification", "full_log2_fold_change", "sensitivity_log2_fold_change",
    "full_padj", "sensitivity_padj", "robustness_status", "claim_ceiling",
]
SPECIES = "Pinus thunbergii"
CONCORDANT = "directionally_concordant"
CLAIM = "Surrogate-reference sensitivity evidence only; not exact P. thunbergii gene evidence or causal resistance"


def read_pthun(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = [row for row in csv.DictReader(handle, delimiter="\t") if row.get("species") == SPECIES]
    return {(row["orthogroup_id"], row["contrast"]): row for row in rows}


def classify(full: str, sensitivity: str) -> str:
    if full == "unavailable" and sensitivity == "unavailable":
        return "unavailable_in_both"
    if full == CONCORDANT and sensitivity == CONCORDANT:
        return "robust_directionally_concordant"
    if full == CONCORDANT and sensitivity != CONCORDANT:
        return "full_model_concordance_not_robust"
    if full != CONCORDANT and sensitivity == CONCORDANT:
        return "sensitivity_only_concordance"
    if full == sensitivity:
        return "stable_nonconcordant_classification"
    return "nonconcordant_classification_changed"


def build_audit(full_path: Path, sensitivity_path: Path, output_path: Path, report_path: Path) -> list[dict[str, str]]:
    full, sensitivity = read_pthun(full_path), read_pthun(sensitivity_path)
    if set(full) != set(sensitivity):
        raise ValueError("Full and sensitivity P. thunbergii candidate keys differ")
    rows = []
    for key in sorted(full):
        original, reduced = full[key], sensitivity[key]
        rows.append({
            "orthogroup_id": key[0],
            "contrast": key[1],
            "mechanism_axes": original.get("mechanism_axes", ""),
            "full_classification": original.get("evidence_classification", ""),
            "sensitivity_classification": reduced.get("evidence_classification", ""),
            "full_log2_fold_change": original.get("cross_species_log2_fold_change", ""),
            "sensitivity_log2_fold_change": reduced.get("cross_species_log2_fold_change", ""),
            "full_padj": original.get("cross_species_padj", ""),
            "sensitivity_padj": reduced.get("cross_species_padj", ""),
            "robustness_status": classify(original.get("evidence_classification", ""), reduced.get("evidence_classification", "")),
            "claim_ceiling": CLAIM,
        })
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    robust = [row for row in rows if row["robustness_status"] == "robust_directionally_concordant"]
    dependent = [row for row in rows if row["robustness_status"] == "full_model_concordance_not_robust"]
    robust_ogs = sorted({row["orthogroup_id"] for row in robust})
    dependent_ogs = sorted({row["orthogroup_id"] for row in dependent})
    lines = [
        "# P. thunbergii sensitivity audit", "",
        "The sensitivity model excludes low-surrogate-mapping baseline run SRR29499878.", "",
        f"- Candidate-contrast rows: {len(rows)}",
        f"- Robust directionally concordant rows: {len(robust)} across {len(robust_ogs)} orthogroups",
        f"- Full-model concordant rows not retained: {len(dependent)} across {len(dependent_ogs)} orthogroups",
        f"- Robust orthogroups: {'; '.join(robust_ogs) or 'none'}",
        f"- Sensitivity-dependent orthogroups: {'; '.join(dependent_ogs) or 'none'}", "",
        "Conservative interpretation: use sensitivity-retained concordance for primary claims; retain lost significance as secondary context because effect directions did not reverse.",
    ]
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", type=Path, required=True)
    parser.add_argument("--sensitivity", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    rows = build_audit(args.full, args.sensitivity, args.output, args.report)
    print(f"Wrote {len(rows)} P. thunbergii sensitivity audit rows")


if __name__ == "__main__":
    main()
