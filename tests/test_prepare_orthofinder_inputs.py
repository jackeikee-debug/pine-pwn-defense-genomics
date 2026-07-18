import importlib.util
from pathlib import Path


def load_script():
    module_path = Path(__file__).resolve().parents[1] / "scripts" / "prepare_orthofinder_inputs.py"
    spec = importlib.util.spec_from_file_location("prepare_orthofinder_inputs", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_fasta(path: Path, record_id: str) -> None:
    path.write_text(f">{record_id}\nMPEPTIDE\n", encoding="utf-8")


def test_prepare_inputs_copies_selected_species_and_writes_manifest(tmp_path: Path) -> None:
    preparer = load_script()
    source_dir = tmp_path / "proteomes"
    output_dir = tmp_path / "orthofinder"
    source_dir.mkdir()
    write_fasta(source_dir / "ptae.faa", "ptae|gene1")
    write_fasta(source_dir / "ptab.faa", "ptab|gene1")
    write_fasta(source_dir / "pabi.faa", "pabi|gene1")

    preparer.prepare_inputs(
        source_dir=source_dir,
        output_dir=output_dir,
        species_ids=["ptae", "ptab"],
        manifest_path=tmp_path / "manifest.tsv",
    )

    assert sorted(path.name for path in output_dir.glob("*.faa")) == ["ptab.faa", "ptae.faa"]
    assert (output_dir / "ptae.faa").read_text(encoding="utf-8").startswith(">ptae|gene1")
    assert (tmp_path / "manifest.tsv").read_text(encoding="utf-8").splitlines() == [
        "species_id\tinput_fasta\torthofinder_fasta",
        f"ptae\t{source_dir / 'ptae.faa'}\t{output_dir / 'ptae.faa'}",
        f"ptab\t{source_dir / 'ptab.faa'}\t{output_dir / 'ptab.faa'}",
    ]


def test_prepare_inputs_reports_missing_species(tmp_path: Path) -> None:
    preparer = load_script()
    source_dir = tmp_path / "proteomes"
    output_dir = tmp_path / "orthofinder"
    source_dir.mkdir()
    write_fasta(source_dir / "ptae.faa", "ptae|gene1")

    try:
        preparer.prepare_inputs(
            source_dir=source_dir,
            output_dir=output_dir,
            species_ids=["ptae", "pden"],
            manifest_path=tmp_path / "manifest.tsv",
        )
    except FileNotFoundError as error:
        assert "pden" in str(error)
    else:
        raise AssertionError("prepare_inputs should fail when a selected species FASTA is missing")
