import importlib.util
from pathlib import Path


def load_script():
    module_path = Path(__file__).resolve().parents[1] / "scripts" / "classify_effectors.py"
    spec = importlib.util.spec_from_file_location("classify_effectors", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_classify_effectors_uses_secretome_and_keyword_evidence(tmp_path: Path) -> None:
    classifier = load_script()
    secretome = tmp_path / "secretome.tsv"
    keywords = tmp_path / "keywords.yaml"
    output = tmp_path / "effectors.tsv"

    secretome.write_text(
        "protein_id\tannotation_text\tlength\tcysteine_fraction\tn_terminal_hydrophobic_count\t"
        "signal_peptide_heuristic\ttransmembrane_domain_heuristic_count\tsmall_secreted_candidate\t"
        "cysteine_rich_candidate\tsecretome_candidate\tevidence\n"
        "BX1\tputative cellulase protein\t240\t0.0100\t12\tyes\t0\tyes\tno\tyes\tn_terminal_signal_peptide_heuristic;small_secreted\n"
        "BX2\tunnamed protein\t180\t0.0800\t10\tyes\t0\tyes\tyes\tyes\tn_terminal_signal_peptide_heuristic;small_secreted;cysteine_rich\n"
        "BX3\tcytosolic protease\t500\t0.0100\t2\tno\t0\tno\tno\tno\tnone\n",
        encoding="utf-8",
    )
    keywords.write_text(
        "effector_classes:\n"
        "  cell_wall_modifying:\n"
        "    keywords:\n"
        "      - cellulase\n"
        "  protease:\n"
        "    keywords:\n"
        "      - protease\n",
        encoding="utf-8",
    )

    rows = classifier.classify_effectors(
        secretome_path=secretome,
        keywords_path=keywords,
        output_path=output,
    )

    assert rows == [
        {
            "protein_id": "BX1",
            "effector_class": "cell_wall_modifying",
            "evidence": "keyword:cellulase;secretome:n_terminal_signal_peptide_heuristic;small_secreted",
            "confidence": "medium",
            "annotation_text": "putative cellulase protein",
            "length": "240",
            "cysteine_fraction": "0.0100",
        },
        {
            "protein_id": "BX2",
            "effector_class": "small_cysteine_rich_secreted",
            "evidence": "secretome:n_terminal_signal_peptide_heuristic;small_secreted;cysteine_rich",
            "confidence": "low",
            "annotation_text": "unnamed protein",
            "length": "180",
            "cysteine_fraction": "0.0800",
        },
    ]
    assert output.read_text(encoding="utf-8").splitlines()[0] == (
        "protein_id\teffector_class\tevidence\tconfidence\tannotation_text\tlength\tcysteine_fraction"
    )


def test_classify_effectors_uses_external_annotations_for_unnamed_secreted_proteins(tmp_path: Path) -> None:
    classifier = load_script()
    secretome = tmp_path / "secretome.tsv"
    keywords = tmp_path / "keywords.yaml"
    annotations = tmp_path / "annotations.tsv"
    output = tmp_path / "effectors.tsv"

    secretome.write_text(
        "protein_id\tannotation_text\tlength\tcysteine_fraction\tn_terminal_hydrophobic_count\t"
        "signal_peptide_heuristic\ttransmembrane_domain_heuristic_count\tsmall_secreted_candidate\t"
        "cysteine_rich_candidate\tsecretome_candidate\tevidence\n"
        "BX1\tunnamed protein product\t240\t0.0100\t12\tyes\t0\tyes\tno\tyes\tn_terminal_signal_peptide_heuristic;small_secreted\n",
        encoding="utf-8",
    )
    keywords.write_text(
        "effector_classes:\n"
        "  cell_wall_modifying:\n"
        "    keywords:\n"
        "      - cellulase\n",
        encoding="utf-8",
    )
    annotations.write_text(
        "gene_id\tannotation_source\tannotation_text\n"
        "BX1\tSwissProt_DIAMOND\tEndoglucanase cellulase OS=Caenorhabditis elegans\n",
        encoding="utf-8",
    )

    rows = classifier.classify_effectors(
        secretome_path=secretome,
        keywords_path=keywords,
        output_path=output,
        annotations_path=annotations,
    )

    assert rows[0]["effector_class"] == "cell_wall_modifying"
    assert rows[0]["confidence"] == "medium"
    assert rows[0]["annotation_text"] == (
        "unnamed protein product | SwissProt_DIAMOND: Endoglucanase cellulase OS=Caenorhabditis elegans"
    )
    assert rows[0]["evidence"].startswith("keyword:cellulase")
