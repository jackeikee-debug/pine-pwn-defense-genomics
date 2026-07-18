import importlib.util
from pathlib import Path


def load_script():
    module_path = Path(__file__).resolve().parents[1] / "scripts" / "prepare_longest_isoforms.py"
    spec = importlib.util.spec_from_file_location("prepare_longest_isoforms", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_selects_longest_isoform_by_numeric_suffix(tmp_path):
    selector = load_script()
    input_fasta = tmp_path / "input.faa"
    output_fasta = tmp_path / "longest.faa"
    summary_tsv = tmp_path / "summary.tsv"
    log_path = tmp_path / "longest.log"
    input_fasta.write_text(
        ">ptab|Pt1G60410.1\nAAAA\n"
        ">ptab|Pt1G60410.2\nAAAAAAAA\n"
        ">ptab|Pt1G60410.3\nAAAAAA\n"
        ">ptab|Pt1G11111\nMMMMM\n"
        ">ptab|Pt1G22222.1\nKKKK\n"
        ">ptab|Pt1G22222.2\nLLLL\n",
        encoding="utf-8",
    )

    summary = selector.select_longest_isoforms(
        input_path=input_fasta,
        output_path=output_fasta,
        summary_path=summary_tsv,
        log_path=log_path,
    )

    assert output_fasta.read_text(encoding="utf-8") == (
        ">ptab|Pt1G60410.2\nAAAAAAAA\n"
        ">ptab|Pt1G11111\nMMMMM\n"
        ">ptab|Pt1G22222.1\nKKKK\n"
    )
    assert summary == {
        "input_records": 6,
        "genes_written": 3,
        "isoforms_removed": 3,
        "multi_isoform_genes": 2,
    }
    assert summary_tsv.read_text(encoding="utf-8").splitlines() == [
        "input_records\tgenes_written\tisoforms_removed\tmulti_isoform_genes",
        "6\t3\t3\t2",
    ]
    assert "Wrote 3 longest isoforms" in log_path.read_text(encoding="utf-8")


def test_gene_id_from_header_removes_only_terminal_numeric_isoform_suffix():
    selector = load_script()

    assert selector.gene_id_from_header("ptab|Pt1G60410.3") == "Pt1G60410"
    assert selector.gene_id_from_header("ptab|Pt1G60410") == "Pt1G60410"
    assert selector.gene_id_from_header("ptab|Pt1G60410.t1") == "Pt1G60410.t1"
