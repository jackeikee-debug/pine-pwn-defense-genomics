#!/usr/bin/env python3
"""Download versioned secretome functional databases with SHA-256 verification."""

from __future__ import annotations

import argparse
import csv
import hashlib
import logging
import shutil
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


REQUIRED_COLUMNS = {
    "database_id",
    "release",
    "url",
    "filename",
    "expected_sha256",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_config(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not REQUIRED_COLUMNS <= set(reader.fieldnames):
            raise ValueError(f"Database config must contain {sorted(REQUIRED_COLUMNS)}")
        rows = list(reader)
    if not rows:
        raise ValueError("Database config contains no records")
    filenames = [row["filename"] for row in rows]
    if len(filenames) != len(set(filenames)):
        raise ValueError("Database config contains duplicate filenames")
    return rows


def download_databases(
    config_path: Path,
    output_dir: Path,
    manifest_path: Path,
) -> list[dict[str, str]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict[str, str]] = []

    for record in read_config(config_path):
        destination = output_dir / record["filename"]
        partial = destination.with_name(destination.name + ".part")
        expected = record["expected_sha256"].strip().lower()
        partial.unlink(missing_ok=True)
        observed = sha256_file(destination) if destination.exists() else ""
        if expected and observed == expected:
            logging.info("Reusing verified %s release %s", record["database_id"], record["release"])
            destination.touch()
        else:
            logging.info("Retrieving %s release %s", record["database_id"], record["release"])
            try:
                with urllib.request.urlopen(record["url"]) as response, partial.open("wb") as handle:
                    shutil.copyfileobj(response, handle, length=1024 * 1024)
                observed = sha256_file(partial)
                if expected and observed != expected:
                    raise ValueError(
                        f"SHA-256 mismatch for {record['database_id']}: "
                        f"expected {expected}, observed {observed}"
                    )
                partial.replace(destination)
            except Exception:
                partial.unlink(missing_ok=True)
                destination.unlink(missing_ok=True)
                raise

        manifest_rows.append(
            {
                "database_id": record["database_id"],
                "release": record["release"],
                "url": record["url"],
                "local_path": destination.as_posix(),
                "bytes": str(destination.stat().st_size),
                "sha256": observed,
                "checksum_status": "verified" if expected else "recorded",
                "retrieved_utc": datetime.now(timezone.utc).isoformat(),
            }
        )

    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)
    logging.info("Retrieved %d functional databases", len(manifest_rows))
    return manifest_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--log", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    handlers: list[logging.Handler] = [logging.StreamHandler()]
    if args.log:
        args.log.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(args.log, encoding="utf-8"))
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=handlers,
    )
    download_databases(args.config, args.output_dir, args.manifest)


if __name__ == "__main__":
    main()
