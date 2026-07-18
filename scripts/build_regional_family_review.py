#!/usr/bin/env python3
"""Build per-orthogroup FASTA files and review tables for top regional contrasts."""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path


SUMMARY_FIELDS = [
    "orthogroup_rank",
    "orthogroup_id",
    "regional_copy_directions",
    "module_ids",
    "matched_keywords",
    "best_abs_log2_ratio",
    "best_log2_ratio",
    "east_asia_gene_count",
    "north_america_gene_count",
    "outgroup_gene_count",
    "total_selected_genes",
    "species_gene_counts",
    "review_flags",
    "fasta_path",
]

MANIFEST_FIELDS = [
    "orthogroup_id",
    "species_id",
    "region",
    "gene_id",
    "sequence_length",
    "source_fasta",
    "orthogroup_fasta",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    header: str | None = None
    parts: list[str] = []
    with opener(path, "rt", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header = line[1:]
                parts = []
            else:
                parts.append(line)
    if header is not None:
        yield header, "".join(parts)


def build_family_review(
    top_orthogroups_path: Path,
    stable_genes_path: Path,
    species_metadata_path: Path,
    proteome_dir: Path,
    output_dir: Path,
    summary_output: Path,
    manifest_output: Path,
    max_orthogroups: int | None = None,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    top_rows = read_tsv(top_orthogroups_path)
    if max_orthogroups is not None:
        top_rows = top_rows[:max_orthogroups]

    top_by_og = {row["orthogroup_id"]: row for row in top_rows}
    species_regions = read_species_regions(species_metadata_path)
    genes_by_og = collect_top_genes(stable_genes_path, set(top_by_og))
    target_ids_by_species = collect_target_ids_by_species(genes_by_og)
    sequences_by_species = read_target_sequences(proteome_dir, target_ids_by_species)

    prepare_output_dir(output_dir)
    summary_rows: list[dict[str, str]] = []
    manifest_rows: list[dict[str, str]] = []

    for top_row in top_rows:
        orthogroup_id = top_row["orthogroup_id"]
        direction = top_row.get("regional_copy_directions", "")
        output_fasta = output_dir / f"{orthogroup_id}.{safe_name(direction)}.faa"
        ordered_genes = order_genes_for_review(genes_by_og.get(orthogroup_id, []), species_regions, direction)
        region_counts = count_by_region(ordered_genes, species_regions)
        species_counts = count_by_species(ordered_genes)

        with output_fasta.open("w", encoding="utf-8") as handle:
            for gene in ordered_genes:
                species_id = gene["species_id"]
                gene_id = gene["gene_id"]
                sequence = sequences_by_species.get(species_id, {}).get(gene_id)
                if sequence is None:
                    raise ValueError(f"Missing sequence for {gene_id} in {proteome_dir / (species_id + '.faa')}")
                region = species_regions.get(species_id, "Unknown")
                handle.write(
                    f">{gene_id} orthogroup={orthogroup_id} species={species_id} region={region}\n"
                    f"{wrap_sequence(sequence)}\n"
                )
                manifest_rows.append(
                    {
                        "orthogroup_id": orthogroup_id,
                        "species_id": species_id,
                        "region": region,
                        "gene_id": gene_id,
                        "sequence_length": str(len(sequence)),
                        "source_fasta": str(proteome_dir / f"{species_id}.faa"),
                        "orthogroup_fasta": str(output_fasta),
                    }
                )

        summary_rows.append(
            {
                "orthogroup_rank": top_row.get("orthogroup_rank", ""),
                "orthogroup_id": orthogroup_id,
                "regional_copy_directions": direction,
                "module_ids": top_row.get("module_ids", ""),
                "matched_keywords": top_row.get("matched_keywords", ""),
                "best_abs_log2_ratio": top_row.get("best_abs_log2_ratio", ""),
                "best_log2_ratio": top_row.get("best_log2_ratio", ""),
                "east_asia_gene_count": str(region_counts["East_Asia"]),
                "north_america_gene_count": str(region_counts["North_America"]),
                "outgroup_gene_count": str(region_counts["Outgroup"]),
                "total_selected_genes": str(len(ordered_genes)),
                "species_gene_counts": format_counts(species_counts),
                "review_flags": build_review_flags(direction, region_counts, species_counts, species_regions),
                "fasta_path": str(output_fasta),
            }
        )

    write_tsv(summary_output, SUMMARY_FIELDS, summary_rows)
    write_tsv(manifest_output, MANIFEST_FIELDS, manifest_rows)
    return summary_rows, manifest_rows


def read_species_regions(path: Path) -> dict[str, str]:
    return {row["species_id"]: row.get("region", "Unknown") for row in read_tsv(path)}


def prepare_output_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)
    for fasta_path in path.glob("*.faa"):
        fasta_path.unlink()


def collect_top_genes(path: Path, orthogroup_ids: set[str]) -> dict[str, list[dict[str, str]]]:
    genes_by_og: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(path):
        if row["orthogroup_id"] in orthogroup_ids:
            genes_by_og[row["orthogroup_id"]].append(row)
    return genes_by_og


def collect_target_ids_by_species(genes_by_og: dict[str, list[dict[str, str]]]) -> dict[str, set[str]]:
    target_ids: dict[str, set[str]] = defaultdict(set)
    for genes in genes_by_og.values():
        for gene in genes:
            target_ids[gene["species_id"]].add(gene["gene_id"])
    return target_ids


def read_target_sequences(proteome_dir: Path, target_ids_by_species: dict[str, set[str]]) -> dict[str, dict[str, str]]:
    sequences: dict[str, dict[str, str]] = {}
    for species_id, target_ids in target_ids_by_species.items():
        fasta_path = proteome_dir / f"{species_id}.faa"
        species_sequences: dict[str, str] = {}
        for header, sequence in read_fasta(fasta_path):
            record_id = header.split()[0]
            if record_id in target_ids:
                species_sequences[record_id] = sequence
        missing = sorted(target_ids - set(species_sequences))
        if missing:
            raise ValueError(f"{fasta_path} is missing {len(missing)} target records, including {missing[0]}")
        sequences[species_id] = species_sequences
    return sequences


def order_genes_for_review(
    genes: list[dict[str, str]],
    species_regions: dict[str, str],
    direction: str,
) -> list[dict[str, str]]:
    region_order = ordered_regions(direction)
    region_rank = {region: index for index, region in enumerate(region_order)}
    return sorted(
        genes,
        key=lambda gene: (
            region_rank.get(species_regions.get(gene["species_id"], "Unknown"), len(region_order)),
            gene["species_id"],
            gene["gene_id"],
        ),
    )


def ordered_regions(direction: str) -> list[str]:
    if direction == "North_America_enriched":
        return ["North_America", "East_Asia", "Outgroup"]
    return ["East_Asia", "North_America", "Outgroup"]


def count_by_region(genes: list[dict[str, str]], species_regions: dict[str, str]) -> dict[str, int]:
    counts = {"East_Asia": 0, "North_America": 0, "Outgroup": 0}
    for gene in genes:
        region = species_regions.get(gene["species_id"], "Unknown")
        if region in counts:
            counts[region] += 1
    return counts


def count_by_species(genes: list[dict[str, str]]) -> dict[str, int]:
    counts: dict[str, int] = defaultdict(int)
    for gene in genes:
        counts[gene["species_id"]] += 1
    return dict(counts)


def build_review_flags(
    direction: str,
    region_counts: dict[str, int],
    species_counts: dict[str, int],
    species_regions: dict[str, str],
) -> str:
    flags: list[str] = []
    if region_counts["East_Asia"] == 0 and region_counts["North_America"] > 0:
        flags.append("North_America_only")
    if region_counts["North_America"] == 0 and region_counts["East_Asia"] > 0:
        flags.append("East_Asia_only")
    if sum(region_counts.values()) >= 30:
        flags.append("large_family")
    if region_counts["Outgroup"] > 0:
        flags.append("outgroup_present")

    focal_region = "North_America" if direction == "North_America_enriched" else "East_Asia"
    focal_counts = [
        count
        for species_id, count in species_counts.items()
        if species_regions.get(species_id) == focal_region
    ]
    focal_total = sum(focal_counts)
    if focal_total >= 3 and focal_counts and max(focal_counts) / focal_total >= 0.75:
        flags.append("single_species_driver")
    return ";".join(flags) if flags else "none"


def format_counts(counts: dict[str, int]) -> str:
    return ";".join(f"{key}:{value}" for key, value in sorted(counts.items()))


def safe_name(value: str) -> str:
    cleaned = "".join(char if char.isalnum() or char in "._-" else "_" for char in value)
    return cleaned or "regional_contrast"


def wrap_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[index : index + width] for index in range(0, len(sequence), width))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--top-orthogroups", type=Path, required=True)
    parser.add_argument("--stable-genes", type=Path, required=True)
    parser.add_argument("--species-metadata", type=Path, required=True)
    parser.add_argument("--proteome-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--max-orthogroups", type=int)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    summary_rows, manifest_rows = build_family_review(
        top_orthogroups_path=args.top_orthogroups,
        stable_genes_path=args.stable_genes,
        species_metadata_path=args.species_metadata,
        proteome_dir=args.proteome_dir,
        output_dir=args.output_dir,
        summary_output=args.summary,
        manifest_output=args.manifest,
        max_orthogroups=args.max_orthogroups,
    )
    print(f"Regional family review orthogroups written: {len(summary_rows)}")
    print(f"Regional family review sequences written: {len(manifest_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
