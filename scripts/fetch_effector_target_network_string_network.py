#!/usr/bin/env python3
"""Fetch STRING mappings and interaction partners for route 1 network seeds."""

from __future__ import annotations

import argparse
import csv
import io
import time
import urllib.parse
import urllib.request
from pathlib import Path


MAPPING_FIELDS = [
    "query_symbol",
    "string_id",
    "preferred_name",
    "annotation",
    "mapping_status",
]

INTERACTION_FIELDS = [
    "stringId_A",
    "stringId_B",
    "preferredName_A",
    "preferredName_B",
    "ncbiTaxonId",
    "score",
    "nscore",
    "fscore",
    "pscore",
    "ascore",
    "escore",
    "dscore",
    "tscore",
]

MANIFEST_FIELDS = [
    "api_url",
    "species",
    "required_score",
    "partner_limit",
    "network_type",
    "caller_identity",
    "query_symbol_count",
    "mapped_symbol_count",
    "interaction_count",
    "status",
    "error",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fetch_effector_target_network_string_network(
    seed_path: Path,
    mapping_output: Path,
    interactions_output: Path,
    manifest_output: Path,
    api_url: str,
    species: str,
    required_score: str,
    partner_limit: str,
    network_type: str,
    caller_identity: str,
    api_post=None,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    api_post = api_post or (lambda method, params: post_string_api(api_url, method, params))
    symbols = seed_symbols(read_tsv(seed_path))
    manifest = {
        "api_url": api_url,
        "species": species,
        "required_score": required_score,
        "partner_limit": partner_limit,
        "network_type": network_type,
        "caller_identity": caller_identity,
        "query_symbol_count": str(len(symbols)),
        "mapped_symbol_count": "0",
        "interaction_count": "0",
        "status": "ok",
        "error": "",
    }
    try:
        mapping_text = api_post(
            "get_string_ids",
            {
                "identifiers": "\r".join(symbols),
                "species": species,
                "echo_query": "1",
                "caller_identity": caller_identity,
            },
        )
        mapping_rows = normalize_mapping_rows(symbols, parse_tsv_text(mapping_text))
        mapped_ids = [row["string_id"] for row in mapping_rows if row["mapping_status"] == "mapped"]
        time.sleep(1)
        interaction_rows: list[dict[str, str]] = []
        if mapped_ids:
            interaction_text = api_post(
                "interaction_partners",
                {
                    "identifiers": "\r".join(mapped_ids),
                    "species": species,
                    "required_score": required_score,
                    "limit": partner_limit,
                    "network_type": network_type,
                    "caller_identity": caller_identity,
                },
            )
            interaction_rows = normalize_interaction_rows(parse_tsv_text(interaction_text))
        manifest["mapped_symbol_count"] = str(len(mapped_ids))
        manifest["interaction_count"] = str(len(interaction_rows))
    except Exception as exc:
        mapping_rows = [{"query_symbol": symbol, "string_id": "", "preferred_name": "", "annotation": "", "mapping_status": "unmapped"} for symbol in symbols]
        interaction_rows = []
        manifest["status"] = "error"
        manifest["error"] = str(exc)
    write_tsv(mapping_output, mapping_rows, MAPPING_FIELDS)
    write_tsv(interactions_output, interaction_rows, INTERACTION_FIELDS)
    write_tsv(manifest_output, [manifest], MANIFEST_FIELDS)
    return mapping_rows, interaction_rows


def seed_symbols(rows: list[dict[str, str]]) -> list[str]:
    symbols: set[str] = set()
    for row in rows:
        symbols.update(split_multi(row.get("arabidopsis_homolog_symbols", "")))
    return sorted(symbols)


def normalize_mapping_rows(symbols: list[str], api_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    by_query = {row.get("queryItem", row.get("preferredName", "")): row for row in api_rows}
    rows = []
    for symbol in symbols:
        api_row = by_query.get(symbol)
        if api_row:
            rows.append(
                {
                    "query_symbol": symbol,
                    "string_id": api_row.get("stringId", ""),
                    "preferred_name": api_row.get("preferredName", ""),
                    "annotation": api_row.get("annotation", ""),
                    "mapping_status": "mapped",
                }
            )
        else:
            rows.append(
                {
                    "query_symbol": symbol,
                    "string_id": "",
                    "preferred_name": "",
                    "annotation": "",
                    "mapping_status": "unmapped",
                }
            )
    return rows


def normalize_interaction_rows(api_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [{field: row.get(field, "") for field in INTERACTION_FIELDS} for row in api_rows]


def parse_tsv_text(text: str) -> list[dict[str, str]]:
    text = text.strip()
    if not text:
        return []
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def post_string_api(api_url: str, method: str, params: dict[str, str]) -> str:
    request_url = "/".join([api_url.rstrip("/"), "tsv", method])
    body = urllib.parse.urlencode(params).encode("utf-8")
    request = urllib.request.Request(request_url, data=body, method="POST")
    request.add_header("Content-Type", "application/x-www-form-urlencoded")
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read().decode("utf-8")


def split_multi(value: str) -> list[str]:
    return [item.strip() for item in value.replace(",", ";").split(";") if item.strip()]


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=Path, required=True)
    parser.add_argument("--mapping-output", type=Path, required=True)
    parser.add_argument("--interactions-output", type=Path, required=True)
    parser.add_argument("--manifest-output", type=Path, required=True)
    parser.add_argument("--api-url", default="https://version-12-0.string-db.org/api")
    parser.add_argument("--species", default="3702")
    parser.add_argument("--required-score", default="700")
    parser.add_argument("--partner-limit", default="20")
    parser.add_argument("--network-type", default="functional")
    parser.add_argument("--caller-identity", default="pine_pwn_effector_guided_comparative_genomics")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    mapping_rows, interaction_rows = fetch_effector_target_network_string_network(
        seed_path=args.seed,
        mapping_output=args.mapping_output,
        interactions_output=args.interactions_output,
        manifest_output=args.manifest_output,
        api_url=args.api_url,
        species=args.species,
        required_score=args.required_score,
        partner_limit=args.partner_limit,
        network_type=args.network_type,
        caller_identity=args.caller_identity,
    )
    mapped = sum(row["mapping_status"] == "mapped" for row in mapping_rows)
    print(f"Route 1 STRING mapping rows: {len(mapping_rows)}; mapped: {mapped}; interactions: {len(interaction_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
