import importlib.util
from pathlib import Path


def load_script():
    module_path = Path(__file__).resolve().parents[1] / "scripts" / "predict_secretome.py"
    spec = importlib.util.spec_from_file_location("predict_secretome", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_predict_secretome_flags_signal_peptide_candidates_and_filters_tm_after_signal(tmp_path: Path) -> None:
    predictor = load_script()
    proteins = tmp_path / "proteins.faa"
    output = tmp_path / "secretome.tsv"

    proteins.write_text(
        ">secreted cellulase-like protein\n"
        "MKKLLLLLLLLLLAAAASAQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ\n"
        ">membrane protein\n"
        "MKKLLLLLLLLLLAAAASAQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ"
        "LLLLLLLLLLLLLLLLLLLLQQQQQQ\n"
        ">cytosolic protein\n"
        "MAQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ\n",
        encoding="utf-8",
    )

    rows = predictor.predict_secretome(
        proteins_path=proteins,
        output_path=output,
        min_length=30,
        small_secreted_max_length=80,
        cysteine_rich_threshold=0.03,
        max_transmembrane_domains=0,
    )

    assert rows[0]["protein_id"] == "secreted"
    assert rows[0]["signal_peptide_heuristic"] == "yes"
    assert rows[0]["transmembrane_domain_heuristic_count"] == "0"
    assert rows[0]["secretome_candidate"] == "yes"
    assert rows[0]["small_secreted_candidate"] == "yes"
    assert rows[1]["protein_id"] == "membrane"
    assert rows[1]["secretome_candidate"] == "no"
    assert rows[1]["transmembrane_domain_heuristic_count"] == "1"
    assert rows[2]["protein_id"] == "cytosolic"
    assert rows[2]["signal_peptide_heuristic"] == "no"
    assert output.read_text(encoding="utf-8").splitlines()[0].startswith(
        "protein_id\tannotation_text\tlength"
    )
