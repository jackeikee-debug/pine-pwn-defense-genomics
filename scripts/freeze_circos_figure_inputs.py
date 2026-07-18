#!/usr/bin/env python3
"""Create or verify a deterministic checksum manifest for Circos figure inputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


DATASET_ORDER = ("candidate_tracks", "effector_context", "edges", "sectors")
MANIFEST_FIELDS = (
    "freeze_version",
    "dataset_id",
    "path",
    "sha256",
    "schema_sha256",
    "row_count",
    "column_count",
)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def inspect_tsv(path: Path) -> dict[str, str]:
    data = path.read_bytes()
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as error:
            raise ValueError(f"Frozen dataset is empty: {path}") from error
        if not header or any(not field for field in header) or len(header) != len(set(header)):
            raise ValueError(f"Frozen dataset has an invalid schema: {path}")
        row_count = sum(1 for row in reader if any(cell != "" for cell in row))
    return {
        "path": path.as_posix(),
        "sha256": sha256_bytes(data),
        "schema_sha256": sha256_bytes("\t".join(header).encode("utf-8")),
        "row_count": str(row_count),
        "column_count": str(len(header)),
    }


def build_manifest(paths: dict[str, Path], freeze_version: str) -> list[dict[str, str]]:
    if set(paths) != set(DATASET_ORDER):
        raise ValueError(f"Expected frozen datasets: {', '.join(DATASET_ORDER)}")
    rows = []
    for dataset_id in DATASET_ORDER:
        row = {"freeze_version": freeze_version, "dataset_id": dataset_id}
        row.update(inspect_tsv(paths[dataset_id]))
        rows.append(row)
    return rows


def verify_manifest(
    manifest_rows: list[dict[str, str]],
    paths: dict[str, Path],
    freeze_version: str,
) -> None:
    counts = {dataset_id: 0 for dataset_id in DATASET_ORDER}
    for row in manifest_rows:
        dataset_id = row.get("dataset_id", "")
        if dataset_id in counts:
            counts[dataset_id] += 1
    if any(counts[dataset_id] != 1 for dataset_id in DATASET_ORDER) or len(manifest_rows) != len(DATASET_ORDER):
        raise ValueError("Freeze manifest must contain exactly one row for each dataset")

    expected = {row["dataset_id"]: row for row in manifest_rows}
    observed = {row["dataset_id"]: row for row in build_manifest(paths, freeze_version)}
    for dataset_id in DATASET_ORDER:
        expected_row = expected[dataset_id]
        observed_row = observed[dataset_id]
        mismatched = [field for field in MANIFEST_FIELDS if expected_row.get(field, "") != observed_row[field]]
        if mismatched:
            raise ValueError(
                f"Frozen dataset drift: {dataset_id} ({', '.join(mismatched)}); "
                "update only this dataset and its manifest row after review"
            )


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("create", "verify"), required=True)
    parser.add_argument("--freeze-version", default="circos_v1")
    parser.add_argument("--candidate-tracks", type=Path, required=True)
    parser.add_argument("--effector-context", type=Path, required=True)
    parser.add_argument("--edges", type=Path, required=True)
    parser.add_argument("--sectors", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    paths = {
        "candidate_tracks": args.candidate_tracks,
        "effector_context": args.effector_context,
        "edges": args.edges,
        "sectors": args.sectors,
    }
    if args.mode == "create":
        write_manifest(args.manifest, build_manifest(paths, args.freeze_version))
        print(f"Created Circos figure-data freeze: {args.manifest}")
    else:
        verify_manifest(read_manifest(args.manifest), paths, args.freeze_version)
        print(f"Verified Circos figure-data freeze: {args.manifest}")


if __name__ == "__main__":
    main()
