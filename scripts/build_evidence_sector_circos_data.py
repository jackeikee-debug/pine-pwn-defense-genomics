#!/usr/bin/env python3
"""Freeze isolated plotting tables for the evidence-sector Circos preview."""

from __future__ import annotations

import argparse
import logging
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd


LOGGER = logging.getLogger("build_evidence_sector_circos_data")

SPECIES = [
    ("pden", "P. densiflora", "East Asia"),
    ("pmas", "P. massoniana", "East Asia"),
    ("ptab", "P. tabuliformis", "East Asia"),
    ("plam", "P. lambertiana", "North America"),
    ("plon", "P. longaeva", "North America"),
    ("ptae", "P. taeda", "North America"),
]

PMAS_CONTRASTS = [
    ("pmas_r", "P. massoniana R", "resistant_pwn_vs_water"),
    ("pmas_s", "P. massoniana S", "susceptible_pwn_vs_water"),
    ("pmas_gxi", "P. massoniana GxI", "genotype_by_inoculum"),
]

TRACK_COLUMNS = [
    "sector_id",
    "item_id",
    "item_label",
    "item_order",
    "track_id",
    "geometry",
    "value",
    "plot_value",
    "fdr",
    "pvalue",
    "state",
    "highlight_tier",
    "module_ids",
    "source_table",
    "evidence_boundary",
    "subitem_id",
    "subitem_order",
    "subitem_count",
    "is_representative",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--effector-context", required=True, type=Path)
    parser.add_argument("--attack-edges", required=True, type=Path)
    parser.add_argument("--cross-species-expression", required=True, type=Path)
    parser.add_argument("--member-audit", required=True, type=Path)
    parser.add_argument("--sectors-output", required=True, type=Path)
    parser.add_argument("--tracks-output", required=True, type=Path)
    parser.add_argument("--links-output", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=True)


def number(value: object) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return math.nan
    return parsed if math.isfinite(parsed) else math.nan


def clean_text(value: object, default: str = "") -> str:
    if pd.isna(value):
        return default
    return str(value).strip()


def candidate_label(row: pd.Series) -> str:
    orthogroup = clean_text(row.get("orthogroup_id"))
    symbols = [part for part in clean_text(row.get("mapped_symbols")).split(";") if part]
    return f"{orthogroup} | {symbols[0]}" if symbols else orthogroup


def split_orthogroups(value: object) -> list[str]:
    return [part.strip() for part in re.split(r"[;,]", clean_text(value)) if part.strip()]


def append_track(
    rows: list[dict[str, object]],
    *,
    sector_id: str,
    item_id: str,
    item_label: str,
    item_order: int,
    track_id: str,
    geometry: str,
    value: object,
    fdr: object = math.nan,
    pvalue: object = math.nan,
    state: str = "observed",
    highlight_tier: str = "context",
    module_ids: str = "",
    source_table: str,
    evidence_boundary: str,
    subitem_id: str = "",
    subitem_order: object = math.nan,
    subitem_count: object = math.nan,
    is_representative: str = "no",
) -> None:
    rows.append(
        {
            "sector_id": sector_id,
            "item_id": item_id,
            "item_label": item_label,
            "item_order": item_order,
            "track_id": track_id,
            "geometry": geometry,
            "value": number(value),
            "plot_value": math.nan,
            "fdr": number(fdr),
            "pvalue": number(pvalue),
            "state": state,
            "highlight_tier": highlight_tier,
            "module_ids": module_ids,
            "source_table": source_table,
            "evidence_boundary": evidence_boundary,
            "subitem_id": subitem_id,
            "subitem_order": number(subitem_order),
            "subitem_count": number(subitem_count),
            "is_representative": is_representative,
        }
    )


def add_expression_tracks(
    rows: list[dict[str, object]],
    *,
    sector_id: str,
    item_id: str,
    item_label: str,
    item_order: int,
    log2fc: object,
    fdr: object,
    pvalue: object,
    highlight_tier: str,
    module_ids: str,
    source_table: str,
    evidence_boundary: str,
) -> None:
    fold_change = number(log2fc)
    adjusted = number(fdr)
    raw_p = number(pvalue)
    state = "significant" if math.isfinite(adjusted) and adjusted < 0.05 else "observed"
    common = dict(
        sector_id=sector_id,
        item_id=item_id,
        item_label=item_label,
        item_order=item_order,
        fdr=adjusted,
        pvalue=raw_p,
        state=state,
        highlight_tier=highlight_tier,
        module_ids=module_ids,
        source_table=source_table,
        evidence_boundary=evidence_boundary,
    )
    append_track(rows, track_id="expression_log2fc", geometry="heatmap", value=fold_change, **common)
    append_track(rows, track_id="expression_abs_log2fc", geometry="bar", value=abs(fold_change), **common)
    append_track(rows, track_id="expression_profile", geometry="line", value=fold_change, **common)
    significance = -math.log10(max(adjusted, 1e-300)) if math.isfinite(adjusted) and adjusted > 0 else math.nan
    append_track(rows, track_id="expression_fdr", geometry="point", value=significance, **common)


def best_rows_by_orthogroup(
    frame: pd.DataFrame,
    *,
    orthogroup_column: str,
    fold_change_column: str,
    fdr_column: str,
    candidate_ids: set[str],
) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for _, row in frame.iterrows():
        for orthogroup in split_orthogroups(row.get(orthogroup_column)):
            if orthogroup not in candidate_ids:
                continue
            record = row.to_dict()
            record["orthogroup_id"] = orthogroup
            record["_fold_change"] = number(row.get(fold_change_column))
            record["_fdr"] = number(row.get(fdr_column))
            records.append(record)
    if not records:
        return pd.DataFrame(columns=["orthogroup_id", "_fold_change", "_fdr"])
    expanded = pd.DataFrame(records)
    expanded["_fdr_sort"] = expanded["_fdr"].fillna(np.inf)
    expanded["_abs_fc_sort"] = expanded["_fold_change"].abs().fillna(-np.inf)
    expanded = expanded.sort_values(
        ["orthogroup_id", "_fdr_sort", "_abs_fc_sort"],
        ascending=[True, True, False],
        kind="mergesort",
    )
    return expanded.groupby("orthogroup_id", sort=False, as_index=False).first()


def best_pstrobus_rows(frame: pd.DataFrame, candidate_ids: set[str]) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for _, row in frame.iterrows():
        orthogroup = clean_text(row.get("orthogroup_id"))
        if orthogroup not in candidate_ids:
            continue
        choices = []
        for label, fc_col, fdr_col in [
            ("2w vs 0w", "two_week_log2_fold_change", "two_week_padj"),
            ("4w vs 0w", "four_week_log2_fold_change", "four_week_padj"),
        ]:
            fc = number(row.get(fc_col))
            fdr = number(row.get(fdr_col))
            choices.append((math.inf if not math.isfinite(fdr) else fdr, -abs(fc) if math.isfinite(fc) else math.inf, label, fc, fdr))
        _, _, label, fc, fdr = sorted(choices)[0]
        records.append(
            {
                "orthogroup_id": orthogroup,
                "_fold_change": fc,
                "_fdr": fdr,
                "_contrast": label,
                "claim_ceiling": clean_text(row.get("claim_ceiling")),
            }
        )
    if not records:
        return pd.DataFrame(columns=["orthogroup_id", "_fold_change", "_fdr", "_contrast", "claim_ceiling"])
    expanded = pd.DataFrame(records)
    expanded["_fdr_sort"] = expanded["_fdr"].fillna(np.inf)
    expanded["_abs_fc_sort"] = expanded["_fold_change"].abs().fillna(-np.inf)
    expanded = expanded.sort_values(
        ["orthogroup_id", "_fdr_sort", "_abs_fc_sort"],
        ascending=[True, True, False],
        kind="mergesort",
    )
    return expanded.groupby("orthogroup_id", sort=False, as_index=False).first()


def scale_plot_values(tracks: pd.DataFrame) -> pd.DataFrame:
    tracks = tracks.copy()
    for track_id, indices in tracks.groupby("track_id").groups.items():
        values = tracks.loc[indices, "value"].astype(float)
        finite = values[np.isfinite(values)]
        if finite.empty:
            continue
        if track_id in {"copy_centered", "expression_log2fc", "expression_profile", "member_log2fc"}:
            limit = max(float(finite.abs().max()), 1e-9)
            tracks.loc[indices, "plot_value"] = (values / (2 * limit) + 0.5).clip(0, 1)
        elif track_id.endswith("support") or track_id == "priority_highlight":
            tracks.loc[indices, "plot_value"] = values.clip(0, 1)
        else:
            maximum = max(float(finite.max()), 1e-9)
            tracks.loc[indices, "plot_value"] = (values / maximum).clip(0, 1)
    return tracks


def build_tables(
    candidates: pd.DataFrame,
    effectors: pd.DataFrame,
    attack_edges: pd.DataFrame,
    cross_species_expression: pd.DataFrame,
    member_audit: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    candidates = candidates.copy()
    candidates["display_order"] = pd.to_numeric(candidates["display_order"], errors="coerce")
    candidates = candidates.sort_values("display_order", kind="mergesort")
    candidates["item_label"] = candidates.apply(candidate_label, axis=1)
    candidate_ids = set(candidates["orthogroup_id"].astype(str))
    candidate_lookup = candidates.set_index("orthogroup_id")

    track_rows: list[dict[str, object]] = []
    sector_rows: list[dict[str, object]] = []

    for sector_order, (species_id, species_label, region) in enumerate(SPECIES, start=1):
        sector_id = f"{species_id}_proteome"
        sector_rows.append(
            {
                "sector_order": sector_order,
                "sector_id": sector_id,
                "sector_label": species_label,
                "sector_group": region,
                "evidence_type": "proteome_copy_context",
                "item_count": len(candidates),
                "evidence_scope": "descriptive orthogroup copy context",
            }
        )
        for _, row in candidates.iterrows():
            common = dict(
                sector_id=sector_id,
                item_id=clean_text(row["orthogroup_id"]),
                item_label=clean_text(row["item_label"]),
                item_order=int(row["display_order"]),
                state=clean_text(row.get("cross_species_state"), "observed"),
                highlight_tier="lead" if clean_text(row.get("priority_class")) == "mechanism_leading" else "context",
                module_ids=clean_text(row.get("module_ids")),
                source_table="manuscript_multilayer_circos_candidates.tsv",
                evidence_boundary="descriptive copy context; not formal expansion or adaptation",
            )
            append_track(track_rows, track_id="copy_centered", geometry="heatmap", value=row.get(f"centered_log2_copy_{species_id}"), **common)
            append_track(track_rows, track_id="copy_count", geometry="bar", value=row.get(f"copy_{species_id}"), **common)
            for support in ["network_support", "structure_support", "sequence_support"]:
                append_track(track_rows, track_id=support, geometry="glyph", value=1 if clean_text(row.get(support)).lower() == "yes" else 0, **common)
            append_track(track_rows, track_id="priority_highlight", geometry="glyph", value=1 if common["highlight_tier"] == "lead" else 0, **common)

    for offset, (sector_id, sector_label, prefix) in enumerate(PMAS_CONTRASTS, start=7):
        sector_rows.append(
            {
                "sector_order": offset,
                "sector_id": sector_id,
                "sector_label": sector_label,
                "sector_group": "P. massoniana expression",
                "evidence_type": "controlled_expression",
                "item_count": len(candidates),
                "evidence_scope": "1 dpi controlled infection contrast",
            }
        )
        for _, row in candidates.iterrows():
            add_expression_tracks(
                track_rows,
                sector_id=sector_id,
                item_id=clean_text(row["orthogroup_id"]),
                item_label=clean_text(row["item_label"]),
                item_order=int(row["display_order"]),
                log2fc=row.get(f"{prefix}_representative_log2_fold_change"),
                fdr=row.get(f"{prefix}_simes_padj"),
                pvalue=row.get(f"{prefix}_simes_pvalue"),
                highlight_tier="lead" if clean_text(row.get("priority_class")) == "mechanism_leading" else "context",
                module_ids=clean_text(row.get("module_ids")),
                source_table="manuscript_multilayer_circos_candidates.tsv",
                evidence_boundary="orthogroup candidate-set expression summary; not direct targeting or causal resistance",
            )

    member_audit = member_audit.copy()
    member_audit = member_audit[
        (member_audit["species"] == "Pinus massoniana")
        & (member_audit["evidence_unit"] == "orthogroup_member")
        & member_audit["orthogroup_id"].isin(candidate_ids)
        & member_audit["feature_id"].notna()
    ].copy()
    member_sector_map = {
        "resistant_pwn_vs_water": "pmas_r",
        "susceptible_pwn_vs_water": "pmas_s",
        "genotype_by_inoculum": "pmas_gxi",
    }
    member_audit = member_audit[member_audit["contrast"].isin(member_sector_map)].copy()
    member_keys = member_audit[
        ["orthogroup_id", "feature_id", "is_displayed_member"]
    ].copy()
    member_keys["member_representative"] = (
        member_keys["is_displayed_member"].fillna("no").str.lower() == "yes"
    )
    member_keys = member_keys.groupby(
        ["orthogroup_id", "feature_id"], as_index=False, sort=False
    )["member_representative"].max()
    member_keys["_representative_sort"] = (~member_keys["member_representative"]).astype(int)
    member_keys = member_keys.sort_values(
        ["orthogroup_id", "_representative_sort", "feature_id"], kind="mergesort"
    )
    member_keys["subitem_order"] = member_keys.groupby("orthogroup_id").cumcount() + 1
    member_keys["subitem_count"] = member_keys.groupby("orthogroup_id")["feature_id"].transform("size")
    member_audit = member_audit.merge(
        member_keys[
            [
                "orthogroup_id",
                "feature_id",
                "subitem_order",
                "subitem_count",
                "member_representative",
            ]
        ],
        on=["orthogroup_id", "feature_id"],
        how="left",
        validate="many_to_one",
    )
    for _, member in member_audit.iterrows():
        orthogroup = clean_text(member["orthogroup_id"])
        candidate = candidate_lookup.loc[orthogroup]
        sector_id = member_sector_map[clean_text(member["contrast"])]
        adjusted = number(member.get("whole_transcriptome_padj"))
        raw_p = number(member.get("raw_pvalue"))
        state = "significant" if math.isfinite(adjusted) and adjusted < 0.05 else "observed"
        common = dict(
            sector_id=sector_id,
            item_id=orthogroup,
            item_label=clean_text(candidate["item_label"]),
            item_order=int(candidate["display_order"]),
            fdr=adjusted,
            pvalue=raw_p,
            state=state,
            highlight_tier="lead" if clean_text(candidate.get("priority_class")) == "mechanism_leading" else "context",
            module_ids=clean_text(candidate.get("module_ids")),
            source_table="manuscript_expression_member_audit.tsv",
            evidence_boundary="P. massoniana orthogroup-member expression; not direct targeting or causal resistance",
            subitem_id=clean_text(member["feature_id"]),
            subitem_order=member["subitem_order"],
            subitem_count=member["subitem_count"],
            is_representative="yes" if bool(member["member_representative"]) else "no",
        )
        append_track(
            track_rows,
            track_id="member_log2fc",
            geometry="heatmap",
            value=member.get("log2_fold_change"),
            **common,
        )
        significance = -math.log10(max(adjusted, 1e-300)) if math.isfinite(adjusted) and adjusted > 0 else math.nan
        append_track(
            track_rows,
            track_id="member_fdr",
            geometry="point",
            value=significance,
            **common,
        )

    cross_species_expression = cross_species_expression.copy()
    cross_species_expression = cross_species_expression[
        cross_species_expression["orthogroup_id"].isin(candidate_ids)
    ].copy()
    cross_species_expression["_fold_change"] = pd.to_numeric(
        cross_species_expression["cross_species_log2_fold_change"], errors="coerce"
    )
    cross_species_expression["_fdr"] = pd.to_numeric(
        cross_species_expression["cross_species_padj"], errors="coerce"
    )
    cross_species_expression["_fdr_sort"] = cross_species_expression["_fdr"].fillna(np.inf)
    cross_species_expression["_abs_fc_sort"] = cross_species_expression["_fold_change"].abs().fillna(-np.inf)
    cross_species_expression = cross_species_expression.sort_values(
        ["species", "orthogroup_id", "_fdr_sort", "_abs_fc_sort"],
        ascending=[True, True, True, False],
        kind="mergesort",
    )
    cross_species_expression = cross_species_expression.groupby(
        ["species", "orthogroup_id"], sort=False, as_index=False
    ).first()
    cross_expression = [
        ("pden_expr", "P. densiflora PWN", "East Asia expression", "Pinus densiflora"),
        ("pthun_expr", "P. thunbergii time", "East Asia expression", "Pinus thunbergii"),
        ("pstrobus_expr", "P. strobus time", "North America expression", "Pinus strobus"),
        (
            "prigtae_expr",
            "P. rigida x taeda",
            "North America expression",
            "Pinus rigida x Pinus taeda",
        ),
    ]
    for offset, (sector_id, label, group, species_name) in enumerate(cross_expression, start=10):
        frame = cross_species_expression[cross_species_expression["species"] == species_name].copy()
        boundary = (
            "cross-species infection-expression summary; homology or time-series context; "
            "not direct targeting or causal resistance"
        )
        sector_rows.append(
            {
                "sector_order": offset,
                "sector_id": sector_id,
                "sector_label": label,
                "sector_group": group,
                "evidence_type": "cross_species_expression",
                "item_count": len(frame),
                "evidence_scope": boundary,
            }
        )
        for _, expression in frame.iterrows():
            orthogroup = clean_text(expression["orthogroup_id"])
            candidate = candidate_lookup.loc[orthogroup]
            add_expression_tracks(
                track_rows,
                sector_id=sector_id,
                item_id=orthogroup,
                item_label=clean_text(candidate["item_label"]),
                item_order=int(candidate["display_order"]),
                log2fc=expression.get("_fold_change"),
                fdr=expression.get("_fdr"),
                pvalue=math.nan,
                highlight_tier="lead" if clean_text(candidate.get("priority_class")) == "mechanism_leading" else "context",
                module_ids=clean_text(candidate.get("module_ids")),
                source_table="pmas_mechanism_cross_species_validation.tsv",
                evidence_boundary=boundary,
            )

    effectors = effectors.copy()
    effectors["display_order"] = pd.to_numeric(effectors["display_order"], errors="coerce")
    effectors = effectors.sort_values("display_order", kind="mergesort")
    sector_rows.append(
        {
            "sector_order": 14,
            "sector_id": "pwn_secretome",
            "sector_label": "PWN secretome",
            "sector_group": "Pine wood nematode",
            "evidence_type": "secretome_functional_context",
            "item_count": len(effectors),
            "evidence_scope": "predicted candidate and class-level annotation context",
        }
    )
    for _, row in effectors.iterrows():
        effector_class = clean_text(row["effector_class"])
        annotation_state = clean_text(row.get("functional_annotation_state"), "unavailable")
        common = dict(
            sector_id="pwn_secretome",
            item_id=effector_class,
            item_label=effector_class.replace("_", " "),
            item_order=int(row["display_order"]),
            state=annotation_state,
            highlight_tier="lead" if annotation_state == "supported" and number(row.get("standard_candidate_count")) > 0 else "context",
            module_ids="",
            source_table="manuscript_circos_effector_class_context.tsv",
            evidence_boundary=clean_text(row.get("claim_boundary"), "class-level context; not direct host targeting"),
        )
        append_track(track_rows, track_id="effector_candidate_count", geometry="bar", value=row.get("standard_candidate_count"), **common)
        append_track(track_rows, track_id="effector_supported_fraction", geometry="heatmap", value=row.get("functional_supported_fraction"), **common)
        append_track(track_rows, track_id="functional_support", geometry="glyph", value=1 if annotation_state == "supported" else 0, **common)
        append_track(track_rows, track_id="priority_highlight", geometry="glyph", value=1 if common["highlight_tier"] == "lead" else 0, **common)

    tracks = pd.DataFrame(track_rows, columns=TRACK_COLUMNS)
    tracks = scale_plot_values(tracks)
    tracks = tracks.sort_values(
        ["sector_id", "item_order", "track_id", "subitem_order"],
        kind="mergesort",
        na_position="last",
    )

    present_items = set(zip(tracks["sector_id"], tracks["item_id"]))
    link_rows: list[dict[str, object]] = []
    for _, edge in attack_edges.iterrows():
        effector_class = clean_text(edge.get("source_id"))
        module = clean_text(edge.get("target_id"))
        if ("pwn_secretome", effector_class) not in present_items:
            continue
        for _, candidate in candidates.iterrows():
            modules = split_orthogroups(candidate.get("module_ids"))
            if module not in modules:
                continue
            link_rows.append(
                {
                    "source_sector": "pwn_secretome",
                    "source_item": effector_class,
                    "target_sector": "pmas_proteome",
                    "target_item": clean_text(candidate["orthogroup_id"]),
                    "link_type": "functional_hypothesis",
                    "weight": number(edge.get("weight")),
                    "highlight_tier": "lead" if clean_text(candidate.get("priority_class")) == "mechanism_leading" else "context",
                    "color_group": module,
                    "evidence_boundary": clean_text(edge.get("evidence_boundary"), "functional hypothesis; not interaction"),
                }
            )

    expression_sectors = [
        "pmas_r",
        "pmas_s",
        "pmas_gxi",
        "pden_expr",
        "pthun_expr",
        "pstrobus_expr",
        "prigtae_expr",
    ]
    for orthogroup in candidates["orthogroup_id"].astype(str):
        for sector_id in expression_sectors:
            subset = tracks[
                (tracks["sector_id"] == sector_id)
                & (tracks["item_id"] == orthogroup)
                & (tracks["track_id"] == "expression_abs_log2fc")
            ]
            if subset.empty:
                continue
            row = subset.iloc[0]
            plotted = number(row["plot_value"])
            if not math.isfinite(plotted):
                continue
            weight = max(plotted, 0.05)
            link_rows.append(
                {
                    "source_sector": "pmas_proteome",
                    "source_item": orthogroup,
                    "target_sector": sector_id,
                    "target_item": orthogroup,
                    "link_type": "orthogroup_evidence",
                    "weight": weight,
                    "highlight_tier": clean_text(row["highlight_tier"], "context"),
                    "color_group": clean_text(row["state"], "observed"),
                    "evidence_boundary": "same orthogroup evidence bridge; not a molecular interaction",
                }
            )

    sectors = pd.DataFrame(sector_rows).sort_values("sector_order", kind="mergesort")
    links = pd.DataFrame(
        link_rows,
        columns=[
            "source_sector",
            "source_item",
            "target_sector",
            "target_item",
            "link_type",
            "weight",
            "highlight_tier",
            "color_group",
            "evidence_boundary",
        ],
    )
    return sectors, tracks, links


def write_table(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, na_rep="NA")


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    sectors, tracks, links = build_tables(
        read_tsv(args.candidates),
        read_tsv(args.effector_context),
        read_tsv(args.attack_edges),
        read_tsv(args.cross_species_expression),
        read_tsv(args.member_audit),
    )
    write_table(sectors, args.sectors_output)
    write_table(tracks, args.tracks_output)
    write_table(links, args.links_output)
    LOGGER.info("Wrote %d sectors, %d track rows, and %d links", len(sectors), len(tracks), len(links))


if __name__ == "__main__":
    main()
