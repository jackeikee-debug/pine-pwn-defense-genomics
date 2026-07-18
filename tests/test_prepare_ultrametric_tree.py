import importlib.util
from pathlib import Path


def load_script():
    module_path = Path(__file__).resolve().parents[1] / "scripts" / "prepare_ultrametric_tree.py"
    spec = importlib.util.spec_from_file_location("prepare_ultrametric_tree", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_table(path: Path, rows: list[list[str]]) -> None:
    header = [
        "node_id",
        "parent_id",
        "node_type",
        "label",
        "age_mya",
        "child_order",
        "evidence_level",
        "citation",
        "notes",
    ]
    path.write_text(
        "\t".join(header) + "\n" + "\n".join("\t".join(row) for row in rows) + "\n",
        encoding="utf-8",
    )


def test_prepare_ultrametric_tree_writes_dated_newick_and_summary(tmp_path: Path) -> None:
    tree_builder = load_script()
    calibrations = tmp_path / "calibrations.tsv"
    manifest = tmp_path / "manifest.tsv"
    tree = tmp_path / "tree.nwk"
    summary = tmp_path / "summary.tsv"

    write_table(
        calibrations,
        [
            ["root", "", "internal", "", "100", "0", "test", "NA", "root"],
            ["ab", "root", "internal", "", "40", "1", "test", "NA", "internal"],
            ["a", "ab", "leaf", "a", "0", "1", "test", "NA", "leaf"],
            ["b", "ab", "leaf", "b", "0", "2", "test", "NA", "leaf"],
            ["c", "root", "leaf", "c", "0", "2", "test", "NA", "leaf"],
        ],
    )
    manifest.write_text(
        "species_id\tinput_fasta\torthofinder_fasta\n"
        "a\ta.faa\ta.faa\n"
        "b\tb.faa\tb.faa\n"
        "c\tc.faa\tc.faa\n",
        encoding="utf-8",
    )

    summary_rows = tree_builder.prepare_ultrametric_tree(
        calibrations_path=calibrations,
        manifest_path=manifest,
        output_path=tree,
        summary_path=summary,
    )

    assert tree.read_text(encoding="utf-8") == "((a:40,b:40):60,c:100);\n"
    assert b"\r\n" not in tree.read_bytes()
    assert summary_rows[0]["root_age_mya"] == "100"
    assert summary_rows[0]["leaf_count"] == "3"
    assert summary_rows[0]["is_ultrametric"] == "true"


def test_prepare_ultrametric_tree_rejects_missing_manifest_species(tmp_path: Path) -> None:
    tree_builder = load_script()
    calibrations = tmp_path / "calibrations.tsv"
    manifest = tmp_path / "manifest.tsv"

    write_table(
        calibrations,
        [
            ["root", "", "internal", "", "10", "0", "test", "NA", "root"],
            ["a", "root", "leaf", "a", "0", "1", "test", "NA", "leaf"],
        ],
    )
    manifest.write_text(
        "species_id\tinput_fasta\torthofinder_fasta\n"
        "a\ta.faa\ta.faa\n"
        "b\tb.faa\tb.faa\n",
        encoding="utf-8",
    )

    try:
        tree_builder.prepare_ultrametric_tree(
            calibrations_path=calibrations,
            manifest_path=manifest,
            output_path=tmp_path / "tree.nwk",
            summary_path=tmp_path / "summary.tsv",
        )
    except ValueError as error:
        assert "missing from calibration leaves: b" in str(error)
    else:
        raise AssertionError("Expected missing species to be rejected")


def test_prepare_ultrametric_tree_rejects_child_older_than_parent(tmp_path: Path) -> None:
    tree_builder = load_script()
    calibrations = tmp_path / "calibrations.tsv"
    manifest = tmp_path / "manifest.tsv"

    write_table(
        calibrations,
        [
            ["root", "", "internal", "", "10", "0", "test", "NA", "root"],
            ["a", "root", "leaf", "a", "20", "1", "test", "NA", "leaf"],
        ],
    )
    manifest.write_text(
        "species_id\tinput_fasta\torthofinder_fasta\n"
        "a\ta.faa\ta.faa\n",
        encoding="utf-8",
    )

    try:
        tree_builder.prepare_ultrametric_tree(
            calibrations_path=calibrations,
            manifest_path=manifest,
            output_path=tmp_path / "tree.nwk",
            summary_path=tmp_path / "summary.tsv",
        )
    except ValueError as error:
        assert "older than parent" in str(error)
    else:
        raise AssertionError("Expected invalid node ages to be rejected")
