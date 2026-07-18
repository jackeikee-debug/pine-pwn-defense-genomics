#!/usr/bin/env python3
"""Summarize MMseqs2 clusters for route 1 effector structure follow-up."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


FIELDS = [
    "representative_id", "cluster_size", "member_ids", "effector_classes",
    "functional_support_levels", "representative_swissprot_subject_id",
    "representative_query_coverage", "representative_bitscore",
    "structure_followup_status", "claim_ceiling",
]
CLAIM = "Sequence-similarity cluster for representative selection only; does not establish shared effector activity or host target"


def read_tsv(path: Path, fieldnames=None) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t", fieldnames=fieldnames))


def summarize(clusters_path: Path, priority_path: Path, output_path: Path, report_path: Path) -> list[dict[str, str]]:
    metadata = {row["protein_id"]: row for row in read_tsv(priority_path)}
    members: dict[str, list[str]] = defaultdict(list)
    for row in read_tsv(clusters_path, ["representative_id", "member_id"]):
        members[row["representative_id"]].append(row["member_id"])
    rows = []
    for representative in sorted(members):
        cluster_members = sorted(set(members[representative]))
        rep = metadata[representative]
        rows.append({
            "representative_id": representative,
            "cluster_size": str(len(cluster_members)),
            "member_ids": ";".join(cluster_members),
            "effector_classes": ";".join(sorted({metadata[item]["effector_class"] for item in cluster_members})),
            "functional_support_levels": ";".join(sorted({metadata[item]["functional_support_level"] for item in cluster_members})),
            "representative_swissprot_subject_id": rep["swissprot_subject_id"],
            "representative_query_coverage": rep["query_coverage"],
            "representative_bitscore": rep["bitscore"],
            "structure_followup_status": "representative_for_structure_or_domain_followup",
            "claim_ceiling": CLAIM,
        })
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    class_counts: dict[str, int] = defaultdict(int)
    for row in rows:
        for effector_class in row["effector_classes"].split(";"):
            class_counts[effector_class] += 1
    lines = [
        "# Route 1 effector sequence clusters", "",
        f"- Candidate sequences: {sum(int(row['cluster_size']) for row in rows)}",
        f"- MMseqs2 clusters and representatives: {len(rows)}",
        "- Parameters: minimum sequence identity 0.50; bidirectional coverage 0.70; greedy cluster mode 2", "",
    ]
    for effector_class in sorted(class_counts):
        lines.append(f"- {effector_class} representative clusters: {class_counts[effector_class]}")
    lines += ["", "Representatives are selected for manageable structure/domain follow-up; clustering does not imply identical activity or host targeting."]
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clusters", type=Path, required=True)
    parser.add_argument("--priority", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    rows = summarize(args.clusters, args.priority, args.output, args.report)
    print(f"Wrote {len(rows)} route 1 effector cluster summaries")


if __name__ == "__main__":
    main()
