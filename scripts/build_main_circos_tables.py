#!/usr/bin/env python3
"""Filter frozen full-Circos tables into the readable main-figure subset."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd


MAIN_SECTORS = (
    "pden_proteome",
    "pmas_proteome",
    "ptab_proteome",
    "plam_proteome",
    "plon_proteome",
    "ptae_proteome",
    "pmas_r",
    "pmas_s",
    "pmas_gxi",
    "pden_expr",
    "pstrobus_expr",
    "pwn_secretome",
)
MOSAIC_SECTORS = ("pmas_r", "pmas_s", "pmas_gxi", "pden_expr", "pstrobus_expr")


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def require_columns(frame: pd.DataFrame, required: set[str], label: str) -> None:
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{label} is missing columns: {', '.join(sorted(missing))}")


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, sep="\t", index=False, na_rep="NA")
    temporary.replace(path)


def build_main_tables(
    sectors: pd.DataFrame,
    tracks: pd.DataFrame,
    links: pd.DataFrame,
    mosaic: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    require_columns(sectors, {"sector_id", "sector_order", "item_count"}, "sectors")
    require_columns(tracks, {"sector_id", "item_id", "track_id"}, "tracks")
    require_columns(
        links,
        {"source_sector", "target_sector", "link_type", "highlight_tier"},
        "links",
    )
    require_columns(
        mosaic,
        {"sector_id", "orthogroup_id", "record_type", "fdr", "significant"},
        "mosaic",
    )

    main_sector_set = set(MAIN_SECTORS)
    filtered_sectors = sectors[sectors["sector_id"].isin(MAIN_SECTORS)].copy()
    if set(filtered_sectors["sector_id"]) != main_sector_set:
        missing = main_sector_set - set(filtered_sectors["sector_id"])
        raise ValueError(f"Full Circos table lacks main sectors: {', '.join(sorted(missing))}")
    filtered_sectors = filtered_sectors.sort_values("sector_order", key=lambda x: x.astype(int))
    filtered_sectors["sector_order"] = range(1, len(filtered_sectors) + 1)

    filtered_tracks = tracks[tracks["sector_id"].isin(main_sector_set)].copy()
    if set(filtered_tracks["sector_id"]) != main_sector_set:
        raise ValueError("Every main sector must retain plotting tracks")

    filtered_links = links[
        links["source_sector"].isin(main_sector_set)
        & links["target_sector"].isin(main_sector_set)
    ].copy()
    filtered_links = filtered_links[
        (filtered_links["link_type"] != "functional_hypothesis")
        | (filtered_links["highlight_tier"] == "lead")
    ].copy()

    filtered_mosaic = mosaic[mosaic["sector_id"].isin(MOSAIC_SECTORS)].copy()
    if set(filtered_mosaic["sector_id"]) != set(MOSAIC_SECTORS):
        raise ValueError("Every core expression dataset must retain gene-level mosaic records")
    if set(filtered_mosaic["record_type"]) != {"sample_expression", "contrast_de"}:
        raise ValueError("Main mosaic requires both sample_expression and contrast_de records")

    core_ogs = set(
        filtered_tracks.loc[filtered_tracks["sector_id"] == "pmas_proteome", "item_id"]
    )
    if len(core_ogs) != 12:
        raise ValueError(f"Expected 12 prioritized orthogroups, found {len(core_ogs)}")
    if not set(filtered_mosaic["orthogroup_id"]).issubset(core_ogs):
        raise ValueError("Main mosaic contains an orthogroup outside the prioritized set")

    return filtered_sectors, filtered_tracks, filtered_links, filtered_mosaic


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sectors", type=Path, required=True)
    parser.add_argument("--tracks", type=Path, required=True)
    parser.add_argument("--links", type=Path, required=True)
    parser.add_argument("--mosaic", type=Path, required=True)
    parser.add_argument("--sectors-output", type=Path, required=True)
    parser.add_argument("--tracks-output", type=Path, required=True)
    parser.add_argument("--links-output", type=Path, required=True)
    parser.add_argument("--mosaic-output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    outputs = build_main_tables(
        read_tsv(args.sectors),
        read_tsv(args.tracks),
        read_tsv(args.links),
        read_tsv(args.mosaic),
    )
    for frame, path in zip(
        outputs,
        (args.sectors_output, args.tracks_output, args.links_output, args.mosaic_output),
        strict=True,
    ):
        write_tsv(frame, path)
        logging.info("Wrote %s rows to %s", len(frame), path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
