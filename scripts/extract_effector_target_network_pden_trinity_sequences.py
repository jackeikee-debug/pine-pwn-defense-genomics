#!/usr/bin/env python3
"""Extract route 1 P. densiflora Trinity transcript sequences from the PMC OA package."""

from __future__ import annotations

import argparse
import csv
import io
import tarfile
import textwrap
import zipfile
from collections import defaultdict
from pathlib import Path


FIELDS = [
    "transcript_id",
    "sequence_found",
    "sequence_count",
    "sequence_lengths",
    "extracted_sequence_ids",
    "orthogroup_ids",
    "effector_target_network_seed_symbols",
    "matched_loci",
    "source_ids",
    "fasta_header",
    "claim_ceiling",
]

CLAIM_CEILING = (
    "Trinity transcript sequence extracted for effector_target_network indirect bridge only; not current Pd gene-ID "
    "mapping or exact differential-expression support"
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def extract_effector_target_network_pden_trinity_sequences(
    bridge_path: Path,
    oa_package_path: Path,
    package_member: str,
    zip_member: str,
    output_fasta: Path,
    manifest_path: Path,
    markdown_path: Path,
    full_reference_path: Path | None = None,
) -> list[dict[str, str]]:
    targets = collect_target_transcripts(read_tsv(bridge_path))
    found = scan_nested_trinity_fasta(
        oa_package_path,
        package_member,
        zip_member,
        set(targets),
        full_reference_path=full_reference_path,
    )
    rows = build_manifest_rows(targets, found)
    write_fasta(output_fasta, found)
    write_tsv(manifest_path, rows)
    write_markdown(markdown_path, rows, output_fasta, full_reference_path)
    return rows


def collect_target_transcripts(rows: list[dict[str, str]]) -> dict[str, dict[str, set[str]]]:
    targets: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    for row in rows:
        if row.get("bridge_status") != "indirect_arabidopsis_locus_bridge_found":
            continue
        if not row.get("source_id", "").startswith("pden"):
            continue
        for transcript_id in split_semicolon(row.get("matched_expression_ids", "")):
            targets[transcript_id]["orthogroup_ids"].add(row.get("orthogroup_id", ""))
            targets[transcript_id]["effector_target_network_seed_symbols"].update(split_semicolon(row.get("effector_target_network_seed_symbols", "")))
            targets[transcript_id]["matched_loci"].update(split_semicolon(row.get("matched_loci", "")))
            targets[transcript_id]["source_ids"].add(row.get("source_id", ""))
    return dict(sorted(targets.items()))


def scan_nested_trinity_fasta(
    oa_package_path: Path,
    package_member: str,
    zip_member: str,
    target_ids: set[str],
    full_reference_path: Path | None = None,
) -> dict[str, dict[str, str]]:
    if not target_ids and full_reference_path is None:
        return {}
    with tarfile.open(oa_package_path, mode="r:gz") as package:
        member = package.extractfile(package_member)
        if member is None:
            raise ValueError(f"Package member not found: {package_member}")
        zip_payload = member.read()
    found: dict[str, dict[str, str]] = {}
    with zipfile.ZipFile(io.BytesIO(zip_payload)) as archive:
        with archive.open(zip_member) as handle:
            if full_reference_path is None:
                parse_fasta_stream(handle, target_ids, found)
            else:
                full_reference_path.parent.mkdir(parents=True, exist_ok=True)
                with full_reference_path.open("wb") as reference_handle:
                    parse_fasta_stream(handle, target_ids, found, reference_handle)
    return found


def parse_fasta_stream(
    handle,
    target_ids: set[str],
    found: dict[str, dict[str, str]],
    full_reference_handle=None,
) -> None:
    current_header = ""
    current_id = ""
    parts: list[str] = []
    for raw_line in handle:
        if full_reference_handle is not None:
            full_reference_handle.write(raw_line)
        line = raw_line.decode("utf-8", errors="replace").strip()
        if not line:
            continue
        if line.startswith(">"):
            store_record(current_id, current_header, parts, target_ids, found)
            current_header = line[1:]
            current_id = current_header.split()[0]
            parts = []
        else:
            parts.append(line)
    store_record(current_id, current_header, parts, target_ids, found)


def store_record(
    transcript_id: str,
    header: str,
    parts: list[str],
    target_ids: set[str],
    found: dict[str, dict[str, str]],
) -> None:
    base_id = trinity_base_id(transcript_id)
    if transcript_id and base_id in target_ids:
        sequence = "".join(parts).upper()
        found.setdefault(base_id, {"records": []})["records"].append(
            {
                "sequence_id": transcript_id,
                "fasta_header": header,
                "sequence": sequence,
            }
        )


def build_manifest_rows(
    targets: dict[str, dict[str, set[str]]],
    found: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    rows = []
    for transcript_id, metadata in targets.items():
        records = found.get(transcript_id, {}).get("records", [])
        sequence_lengths = sorted(len(record["sequence"]) for record in records)
        rows.append(
            {
                "transcript_id": transcript_id,
                "sequence_found": "yes" if records else "no",
                "sequence_count": str(len(records)),
                "sequence_lengths": ";".join(str(value) for value in sequence_lengths),
                "extracted_sequence_ids": ";".join(sorted(record["sequence_id"] for record in records)),
                "orthogroup_ids": join_sorted(metadata["orthogroup_ids"]),
                "effector_target_network_seed_symbols": join_sorted(metadata["effector_target_network_seed_symbols"]),
                "matched_loci": join_sorted(metadata["matched_loci"]),
                "source_ids": join_sorted(metadata["source_ids"]),
                "fasta_header": " || ".join(record["fasta_header"] for record in records),
                "claim_ceiling": CLAIM_CEILING,
            }
        )
    return rows


def write_fasta(path: Path, found: dict[str, dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for transcript_id, payload in sorted(found.items()):
            for record in sorted(payload.get("records", []), key=lambda item: item["sequence_id"]):
                handle.write(f">{record['fasta_header']}\n")
                handle.write("\n".join(textwrap.wrap(record["sequence"], width=80)) + "\n")


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(
    path: Path,
    rows: list[dict[str, str]],
    fasta_path: Path,
    full_reference_path: Path | None = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    found_rows = [row for row in rows if row["sequence_found"] == "yes"]
    missing_rows = [row for row in rows if row["sequence_found"] != "yes"]
    extracted_record_count = sum(int(row["sequence_count"]) for row in found_rows)
    lines = [
        "# Route 1 P. densiflora Trinity Sequence Extraction",
        "",
        "This report summarizes Trinity transcript sequences extracted for route 1 indirect expression bridges.",
        "",
        f"- Target Trinity transcript IDs: {len(rows)}",
        f"- Target IDs with at least one sequence: {len(found_rows)}",
        f"- Isoform sequence records extracted: {extracted_record_count}",
        f"- Missing target sequences: {len(missing_rows)}",
        f"- Output FASTA: `{fasta_path}`",
        f"- Published full Trinity reference: `{full_reference_path}`" if full_reference_path else "",
        "",
        "Evidence boundary: extracted transcript sequences support downstream mapping only; they are not exact current Pd gene IDs or DEG calls.",
        "",
    ]
    by_orthogroup = defaultdict(int)
    for row in found_rows:
        for orthogroup_id in split_semicolon(row["orthogroup_ids"]):
            by_orthogroup[orthogroup_id] += 1
    if by_orthogroup:
        lines.extend(["## Extracted Sequences by Orthogroup", ""])
        for orthogroup_id, count in sorted(by_orthogroup.items()):
            lines.append(f"- {orthogroup_id}: {count}")
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def trinity_base_id(transcript_id: str) -> str:
    if "_i" not in transcript_id:
        return transcript_id
    prefix, suffix = transcript_id.rsplit("_i", 1)
    if suffix.isdigit():
        return prefix
    return transcript_id


def split_semicolon(value: str) -> set[str]:
    return {item.strip() for item in value.split(";") if item.strip()}


def join_sorted(values: set[str]) -> str:
    return ";".join(sorted(item for item in values if item))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bridge", type=Path, required=True)
    parser.add_argument("--oa-package", type=Path, required=True)
    parser.add_argument("--package-member", required=True)
    parser.add_argument("--zip-member", required=True)
    parser.add_argument("--output-fasta", type=Path, required=True)
    parser.add_argument("--full-reference-output", type=Path)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = extract_effector_target_network_pden_trinity_sequences(
        bridge_path=args.bridge,
        oa_package_path=args.oa_package,
        package_member=args.package_member,
        zip_member=args.zip_member,
        output_fasta=args.output_fasta,
        manifest_path=args.manifest,
        markdown_path=args.markdown,
        full_reference_path=args.full_reference_output,
    )
    found = sum(row["sequence_found"] == "yes" for row in rows)
    print(f"Route 1 P. densiflora Trinity targets: {len(rows)}; extracted: {found}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
