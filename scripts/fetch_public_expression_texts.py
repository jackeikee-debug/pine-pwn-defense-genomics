#!/usr/bin/env python3
"""Fetch public expression evidence pages and cache plain text."""

from __future__ import annotations

import argparse
import csv
import html
import re
from datetime import date
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


MANIFEST_FIELDS = [
    "source_id",
    "url",
    "status",
    "text_path",
    "raw_bytes",
    "text_chars",
    "fetch_date",
    "error",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fetch_public_expression_texts(
    sources_path: Path,
    output_dir: Path,
    manifest_path: Path,
    timeout_seconds: int = 30,
) -> list[dict[str, str]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for source in read_tsv(sources_path):
        rows.append(fetch_one_source(source, output_dir, timeout_seconds))
    write_manifest(manifest_path, rows)
    return rows


def fetch_one_source(source: dict[str, str], output_dir: Path, timeout_seconds: int) -> dict[str, str]:
    source_id = source["source_id"]
    url = source.get("url", "")
    text_path = output_dir / f"{source_id}.txt"
    try:
        raw = download_url(url, timeout_seconds)
        text = html_to_text(raw.decode("utf-8", errors="replace"))
        text_path.write_text(text, encoding="utf-8")
        return manifest_row(source_id, url, "ok", text_path, len(raw), len(text), "")
    except (HTTPError, URLError, TimeoutError, ValueError, OSError) as exc:
        return manifest_row(source_id, url, "error", text_path, 0, 0, str(exc))


def download_url(url: str, timeout_seconds: int) -> bytes:
    if not url:
        raise ValueError("missing URL")
    request = Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urlopen(request, timeout=timeout_seconds) as response:
        return response.read()


def html_to_text(markup: str) -> str:
    text = re.sub(r"(?is)<(script|style).*?</\1>", " ", markup)
    text = re.sub(r"(?s)<[^>]+>", " ", text)
    text = html.unescape(text)
    return re.sub(r"\s+", " ", text).strip()


def manifest_row(
    source_id: str,
    url: str,
    status: str,
    text_path: Path,
    raw_bytes: int,
    text_chars: int,
    error: str,
) -> dict[str, str]:
    return {
        "source_id": source_id,
        "url": url,
        "status": status,
        "text_path": str(text_path),
        "raw_bytes": str(raw_bytes),
        "text_chars": str(text_chars),
        "fetch_date": date.today().isoformat(),
        "error": error,
    }


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--timeout-seconds", type=int, default=30)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = fetch_public_expression_texts(
        sources_path=args.sources,
        output_dir=args.output_dir,
        manifest_path=args.manifest,
        timeout_seconds=args.timeout_seconds,
    )
    ok_count = sum(row["status"] == "ok" for row in rows)
    print(f"Expression source texts fetched: {ok_count}/{len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
