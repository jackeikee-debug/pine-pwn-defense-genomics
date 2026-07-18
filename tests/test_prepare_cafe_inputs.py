import importlib.util
from pathlib import Path


def load_script():
    module_path = Path(__file__).resolve().parents[1] / "scripts" / "prepare_cafe_inputs.py"
    spec = importlib.util.spec_from_file_location("prepare_cafe_inputs", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_prepare_cafe_inputs_counts_orthogroup_members_and_filters_large_families(tmp_path: Path) -> None:
    cafe = load_script()
    orthogroups = tmp_path / "Orthogroups.tsv"
    manifest = tmp_path / "manifest.tsv"
    tree = tmp_path / "tree.txt"
    output = tmp_path / "cafe_counts.tsv"
    summary = tmp_path / "summary.tsv"
    tree_output = tmp_path / "cafe_tree.nwk"

    orthogroups.write_text(
        "Orthogroup\tpden\tptae\tpabi\n"
        "OG1\tpden|a, pden|b\tptae|a\tpabi|a\n"
        "OG2\t\tptae|b, ptae|c, ptae|d\t\n"
        "OG3\tpden|c\t\tpabi|b, pabi|c, pabi|d, pabi|e\n",
        encoding="utf-8",
    )
    manifest.write_text(
        "species_id\tinput_fasta\torthofinder_fasta\n"
        "pden\tpden.faa\tpden.faa\n"
        "ptae\tptae.faa\tptae.faa\n"
        "pabi\tpabi.faa\tpabi.faa\n",
        encoding="utf-8",
    )
    tree.write_text("((pden:1,ptae:1):1,pabi:2);\n", encoding="utf-8")

    rows, summary_rows = cafe.prepare_cafe_inputs(
        orthogroups_path=orthogroups,
        manifest_path=manifest,
        tree_path=tree,
        output_path=output,
        summary_path=summary,
        tree_output_path=tree_output,
        max_family_size=4,
        min_species_present=2,
    )

    assert [row["Family ID"] for row in rows] == ["OG1"]
    assert rows[0]["Desc"] == "(null)"
    assert rows[0]["pden"] == "2"
    assert rows[0]["ptae"] == "1"
    assert rows[0]["pabi"] == "1"
    assert summary_rows[0]["total_orthogroups"] == "3"
    assert summary_rows[0]["written_families"] == "1"
    assert summary_rows[0]["filtered_by_max_family_size"] == "1"
    assert summary_rows[0]["filtered_by_min_species_present"] == "1"
    assert output.read_text(encoding="utf-8").splitlines()[0] == "Desc\tFamily ID\tpden\tptae\tpabi"
    assert b"\r\n" not in output.read_bytes()
    assert b"\r\n" not in summary.read_bytes()
    assert tree_output.read_text(encoding="utf-8") == "((pden:1,ptae:1):1,pabi:2);\n"
