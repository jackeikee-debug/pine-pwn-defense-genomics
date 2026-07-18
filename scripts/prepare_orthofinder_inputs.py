#!/usr/bin/env python3
"""Prepare a frozen FASTA directory for an OrthoFinder run."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def prepare_inputs(
    source_dir: Path,
    output_dir: Path,
    species_ids: list[str],
    manifest_path: Path,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    missing = [species_id for species_id in species_ids if not (source_dir / f"{species_id}.faa").exists()]
    if missing:
        raise FileNotFoundError(f"Missing longest-isoform FASTA for species: {', '.join(missing)}")

    manifest_lines = ["species_id\tinput_fasta\torthofinder_fasta"]
    selected_names = {f"{species_id}.faa" for species_id in species_ids}
    for existing in output_dir.glob("*.faa"):
        if existing.name not in selected_names:
            existing.unlink()

    for species_id in species_ids:
        input_fasta = source_dir / f"{species_id}.faa"
        output_fasta = output_dir / f"{species_id}.faa"
        shutil.copy2(input_fasta, output_fasta)
        manifest_lines.append(f"{species_id}\t{input_fasta}\t{output_fasta}")

    manifest_path.write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True, help="Directory with processed longest-isoform FASTA files.")
    parser.add_argument("--output-dir", type=Path, required=True, help="Output directory for OrthoFinder FASTA inputs.")
    parser.add_argument("--species", nargs="+", required=True, help="Species IDs to include.")
    parser.add_argument("--manifest", type=Path, required=True, help="Manifest TSV to write.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    prepare_inputs(
        source_dir=args.source_dir,
        output_dir=args.output_dir,
        species_ids=args.species,
        manifest_path=args.manifest,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
