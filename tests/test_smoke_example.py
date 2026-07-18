import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_reduced_example_matches_expected_counts():
    expected_path = ROOT / "examples/smoke/expected/proteome_summary.tsv"
    with expected_path.open(encoding="utf-8", newline="") as handle:
        expected = {row["species_id"]: int(row["protein_count"]) for row in csv.DictReader(handle, delimiter="\t")}
    observed = {}
    for path in (ROOT / "examples/smoke/proteins").glob("*.faa"):
        observed[path.stem] = sum(line.startswith(">") for line in path.read_text().splitlines())
    assert observed == expected

