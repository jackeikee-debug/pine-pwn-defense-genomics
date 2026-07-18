#!/usr/bin/env python3
"""Download versioned published-secretome tables and legacy protein sequences."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import logging
import time
import urllib.request
import zipfile
from pathlib import Path


FIELDS = ["source_id", "url", "archive_member", "output_path", "sha256", "status"]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def fetch(url: str, attempts: int = 3) -> bytes:
    error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(url, timeout=120) as response:
                return response.read()
        except Exception as exc:  # pragma: no cover - exercised by live downloads
            error = exc
            if attempt < attempts:
                time.sleep(attempt * 2)
    raise RuntimeError(f"Could not download {url}: {error}")


def resolve_output(raw_path: str, config_path: Path) -> Path:
    path = Path(raw_path)
    return path if path.is_absolute() else Path.cwd() / path


def download_sources(config_path: Path, manifest_path: Path) -> list[dict[str, str]]:
    rows = read_tsv(config_path)
    cache: dict[str, bytes] = {}
    manifest: list[dict[str, str]] = []
    for row in rows:
        output = resolve_output(row["output_path"], config_path)
        expected = row["sha256"].lower()
        if output.exists() and digest(output.read_bytes()) == expected:
            status = "reused"
            output.touch()
        else:
            payload = cache.get(row["url"])
            if payload is None:
                payload = fetch(row["url"])
                cache[row["url"]] = payload
            member = row.get("archive_member", "").strip()
            if member and member.upper() != "NA":
                with zipfile.ZipFile(io.BytesIO(payload)) as archive:
                    try:
                        data = archive.read(member)
                    except KeyError as exc:
                        raise ValueError(f"Missing archive member {member} in {row['url']}") from exc
            else:
                data = payload
            observed = digest(data)
            if observed != expected:
                raise ValueError(
                    f"Checksum mismatch for {row['source_id']}: expected {expected}, observed {observed}"
                )
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_bytes(data)
            status = "downloaded"
        manifest.append({**row, "status": status})
        logging.info("%s %s", status, output)

    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(manifest)
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--log", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    handlers: list[logging.Handler] = [logging.StreamHandler()]
    if args.log:
        args.log.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(args.log, encoding="utf-8"))
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", handlers=handlers)
    download_sources(args.config, args.manifest)


if __name__ == "__main__":
    main()
