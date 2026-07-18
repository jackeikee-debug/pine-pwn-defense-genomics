#!/usr/bin/env python3
"""Fetch public expression supplementary files and record a manifest."""

from __future__ import annotations

import argparse
import csv
import io
import tarfile
from pathlib import Path
from http.client import IncompleteRead
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse
from urllib.request import Request, urlopen


FIELDS = [
    "source_id",
    "url",
    "status",
    "local_path",
    "file_size_bytes",
    "error",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fetch_public_expression_supplements(
    sources_path: Path,
    output_dir: Path,
    manifest_path: Path,
    timeout_seconds: int = 60,
    downloader=None,
) -> list[dict[str, str]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    download = downloader or download_url
    rows = []
    for source in read_tsv(sources_path):
        rows.append(fetch_one(source, output_dir, timeout_seconds, download))
    write_tsv(manifest_path, rows)
    return rows


def fetch_one(source: dict[str, str], output_dir: Path, timeout_seconds: int, downloader) -> dict[str, str]:
    source_id = source["source_id"]
    url = source["url"]
    local_path = output_dir / local_filename(source_id, url)
    try:
        payload = downloader(url, timeout_seconds)
        validation_error = validate_payload(url, payload)
        if validation_error:
            package_payload = fetch_from_oa_package(source, timeout_seconds, downloader)
            if package_payload:
                package_error = validate_payload(url, package_payload)
                if not package_error:
                    local_path.write_bytes(package_payload)
                    return manifest_row(source_id, url, "ok_from_oa_package", local_path, len(package_payload), "")
            local_path.unlink(missing_ok=True)
            return manifest_row(source_id, url, "invalid_content", local_path, 0, validation_error)
        local_path.write_bytes(payload)
        return manifest_row(source_id, url, "ok", local_path, len(payload), "")
    except (HTTPError, URLError, TimeoutError, IncompleteRead, ValueError, OSError, tarfile.TarError) as exc:
        local_path.unlink(missing_ok=True)
        return manifest_row(source_id, url, "error", local_path, 0, str(exc))


def download_url(url: str, timeout_seconds: int) -> bytes:
    request = Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urlopen(request, timeout=timeout_seconds) as response:
        return response.read()


def fetch_from_oa_package(source: dict[str, str], timeout_seconds: int, downloader) -> bytes | None:
    package_url = source.get("oa_package_url", "")
    member_name = source.get("package_member", "")
    local_package_path = source.get("oa_package_local_path", "")
    if not member_name:
        return None
    if local_package_path and Path(local_package_path).exists():
        package = Path(local_package_path).read_bytes()
    elif package_url:
        package = downloader(package_url, timeout_seconds)
    else:
        return None
    with tarfile.open(fileobj=io.BytesIO(package), mode="r:gz") as archive:
        member = archive.extractfile(member_name)
        if member is None:
            raise ValueError(f"Package member not found: {member_name}")
        return member.read()


def local_filename(source_id: str, url: str) -> str:
    suffix = Path(urlparse(url).path).suffix or ".dat"
    return f"{source_id}{suffix}"


def validate_payload(url: str, payload: bytes) -> str:
    suffix = Path(urlparse(url).path).suffix.lower()
    if suffix in {".xlsx", ".docx"} and not payload.startswith(b"PK"):
        return f"Expected {suffix.lstrip('.')} ZIP payload, received non-Office content"
    if payload.lstrip().lower().startswith(b"<html"):
        return "Received HTML page instead of supplementary file"
    return ""


def manifest_row(
    source_id: str,
    url: str,
    status: str,
    local_path: Path,
    file_size_bytes: int,
    error: str,
) -> dict[str, str]:
    return {
        "source_id": source_id,
        "url": url,
        "status": status,
        "local_path": str(local_path),
        "file_size_bytes": str(file_size_bytes),
        "error": error,
    }


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--timeout-seconds", type=int, default=60)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = fetch_public_expression_supplements(
        sources_path=args.sources,
        output_dir=args.output_dir,
        manifest_path=args.manifest,
        timeout_seconds=args.timeout_seconds,
    )
    ok_count = sum(row["status"] == "ok" for row in rows)
    print(f"Expression supplementary files fetched: {ok_count}/{len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
