#!/usr/bin/env python3
"""Build route 1 effector-target network seed table."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


FIELDS = [
    "orthogroup_id",
    "network_priority",
    "effector_target_network_role",
    "matched_keyword",
    "host_module",
    "host_network_role",
    "effector_class",
    "effector_count",
    "effector_ids",
    "effector_evidence",
    "link_type",
    "link_confidence",
    "pine_gene_count",
    "pine_gene_ids",
    "arabidopsis_homolog_symbols",
    "arabidopsis_swissprot_ids",
    "other_plant_swissprot_ids",
    "expression_support_level",
    "literature_hit_sources",
    "claim_ceiling",
    "network_next_step",
]

FOCAL_NETWORK_OGS = {"OG0003130", "OG0003639", "OG0002098"}
ARATH_ID_RE = re.compile(r"\|([A-Z0-9]+_ARATH)\b")
GENE_RE = re.compile(r"\bGN=([A-Za-z0-9_.-]+)")
SWISSPROT_ID_RE = re.compile(r"\|([A-Z0-9]+_[A-Z0-9]+)\b")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_effector_target_network_network_seed(
    effector_target_network_path: Path,
    links_path: Path,
    output_path: Path,
    markdown_path: Path,
) -> list[dict[str, str]]:
    candidates = [
        row for row in read_tsv(effector_target_network_path)
        if row.get("deprioritized_in_effector_target_network") != "yes"
    ]
    links = read_tsv(links_path)
    rows = []
    for candidate in candidates:
        rows.extend(build_candidate_rows(candidate, links))
    rows = sorted(rows, key=lambda row: (row["network_priority"], row["orthogroup_id"], row["effector_class"]))
    write_tsv(output_path, rows)
    write_markdown(markdown_path, rows)
    return rows


def build_candidate_rows(candidate: dict[str, str], links: list[dict[str, str]]) -> list[dict[str, str]]:
    rows = []
    orthogroup_id = candidate["orthogroup_id"]
    modules = split_multi(candidate.get("module_ids", ""))
    for link in links:
        if link.get("host_module") not in modules:
            continue
        if orthogroup_id not in split_multi(link.get("host_orthogroups", "")):
            continue
        rows.append(build_row(candidate, link))
    return rows


def build_row(candidate: dict[str, str], link: dict[str, str]) -> dict[str, str]:
    pine_gene_ids = candidate_gene_ids(candidate)
    arath_ids, arath_symbols, other_ids = homolog_summary(candidate.get("candidate_gene_annotations", ""))
    return {
        "orthogroup_id": candidate["orthogroup_id"],
        "network_priority": network_priority(candidate),
        "effector_target_network_role": candidate.get("effector_target_network_role", ""),
        "matched_keyword": candidate.get("matched_keyword", ""),
        "host_module": link.get("host_module", ""),
        "host_network_role": host_network_role(candidate.get("matched_keyword", ""), link.get("host_module", "")),
        "effector_class": link.get("effector_class", ""),
        "effector_count": link.get("effector_count", ""),
        "effector_ids": link.get("effector_ids", ""),
        "effector_evidence": link.get("effector_evidence", ""),
        "link_type": link.get("link_type", ""),
        "link_confidence": link.get("link_confidence", ""),
        "pine_gene_count": str(len(pine_gene_ids)),
        "pine_gene_ids": ";".join(pine_gene_ids),
        "arabidopsis_homolog_symbols": ";".join(arath_symbols),
        "arabidopsis_swissprot_ids": ";".join(arath_ids),
        "other_plant_swissprot_ids": ";".join(other_ids),
        "expression_support_level": candidate.get("expression_support_level", ""),
        "literature_hit_sources": candidate.get("literature_hit_sources", ""),
        "claim_ceiling": claim_ceiling(candidate),
        "network_next_step": "map Arabidopsis homologs into plant PPI/interolog network and compute hub/bottleneck centrality",
    }


def candidate_gene_ids(candidate: dict[str, str]) -> list[str]:
    values = []
    for field in ["pden_gene_ids", "pmas_gene_ids", "ptab_gene_ids"]:
        values.extend(split_multi(candidate.get(field, "")))
    return sorted(set(values))


def homolog_summary(annotation_blob: str) -> tuple[list[str], list[str], list[str]]:
    annotations = split_semicolon(annotation_blob)
    arath_ids: set[str] = set()
    arath_symbols: set[str] = set()
    other_ids: set[str] = set()
    for annotation in annotations:
        swissprot_ids = set(SWISSPROT_ID_RE.findall(annotation))
        arath_matches = set(ARATH_ID_RE.findall(annotation))
        arath_ids.update(arath_matches)
        other_ids.update(swissprot_ids - arath_matches)
        if arath_matches:
            arath_symbols.update(GENE_RE.findall(annotation))
    return sorted(arath_ids), sorted(arath_symbols), sorted(other_ids)


def network_priority(candidate: dict[str, str]) -> str:
    if candidate.get("orthogroup_id") in FOCAL_NETWORK_OGS:
        return "primary_network_seed"
    if candidate.get("effector_target_network_role") == "primary":
        return "primary_network_seed"
    return "secondary_network_context"


def host_network_role(keyword: str, host_module: str) -> str:
    lowered = keyword.lower()
    if any(term in lowered for term in ["wrky", "erf", "myb", "nac"]):
        return "transcriptional_regulatory_hub_candidate"
    if "superoxide dismutase" in lowered or host_module == "ros_detoxification":
        return "oxidative_response_hub_candidate"
    return "defense_module_candidate"


def claim_ceiling(candidate: dict[str, str]) -> str:
    pieces = [
        "predicted_interolog_or_network_overlay only",
        "no direct effector-target interaction evidence",
    ]
    existing = candidate.get("claim_ceiling", "")
    if existing:
        pieces.append(existing)
    return "; ".join(pieces)


def split_multi(value: str) -> list[str]:
    return [item.strip() for item in value.replace(",", ";").split(";") if item.strip()]


def split_semicolon(value: str) -> list[str]:
    return [item.strip() for item in value.split(";") if item.strip()]


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    primary = [row for row in rows if row["network_priority"] == "primary_network_seed"]
    arath_supported = [row for row in rows if row["arabidopsis_swissprot_ids"]]
    lines = [
        "# Route 1 Effector-Target Network Seed",
        "",
        "This table is the first network-analysis layer for the effector-target network hypothesis.",
        "It links predicted PWN effector classes to route 1 pine defense candidates through existing module mappings and SwissProt/Arabidopsis homolog annotations.",
        "",
        f"- Network seed rows: {len(rows)}",
        f"- Primary WRKY/ERF/SOD seed rows: {len(primary)}",
        f"- Rows with Arabidopsis homolog IDs: {len(arath_supported)}",
        "",
        "Interpretation boundary: these are predicted network seeds for interolog/PPI overlay, not validated effector targets.",
        "",
        "| orthogroup | effector class | host module | keyword | Arabidopsis homologs | priority |",
        "|---|---|---|---|---|---|",
    ]
    for row in rows:
        lines.append(
            "| {orthogroup_id} | {effector_class} | {host_module} | {matched_keyword} | {arabidopsis_homolog_symbols} | {network_priority} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "Next step: map the Arabidopsis homolog symbols or SwissProt IDs into a plant PPI/interolog network, then compute degree, betweenness, and bottleneck centrality before overlaying infection-expression evidence.",
        ]
    )
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--effector_target_network", type=Path, required=True)
    parser.add_argument("--links", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = build_effector_target_network_network_seed(
        effector_target_network_path=args.effector_target_network,
        links_path=args.links,
        output_path=args.output,
        markdown_path=args.markdown,
    )
    print(f"Route 1 network seed rows written: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
