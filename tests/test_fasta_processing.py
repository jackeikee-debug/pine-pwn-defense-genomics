import importlib.util
import gzip
from pathlib import Path


def load_script(script_name):
    module_path = Path(__file__).resolve().parents[1] / "scripts" / script_name
    spec = importlib.util.spec_from_file_location(script_name.replace(".py", ""), module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_clean_fasta_standardizes_headers_filters_sequences_and_writes_summary(tmp_path):
    cleaner = load_script("clean_fasta_headers.py")
    input_fasta = tmp_path / "raw.faa"
    output_fasta = tmp_path / "processed.faa"
    summary_tsv = tmp_path / "summary.tsv"
    log_path = tmp_path / "clean.log"
    input_fasta.write_text(
        ">geneA.1 description text\nMKTAA*\n"
        ">geneA.1 duplicate should be removed\nMKTAAAA\n"
        ">geneB short\nMK\n"
        ">geneC internal stop\nMA*TA\n"
        ">sp|Q1|geneD extra\nMKKKK\n",
        encoding="utf-8",
    )

    summary = cleaner.clean_fasta(
        input_path=input_fasta,
        output_path=output_fasta,
        species_id="ptae",
        min_length=5,
        summary_path=summary_tsv,
        log_path=log_path,
    )

    assert output_fasta.read_text(encoding="utf-8") == (
        ">ptae|geneA.1\nMKTAA\n"
        ">ptae|geneC\nMAXTA\n"
        ">ptae|Q1\nMKKKK\n"
    )
    assert summary["proteins_written"] == 3
    assert summary["duplicated_id_count"] == 1
    assert summary["short_sequence_count"] == 1
    assert summary["terminal_stop_count"] == 1
    assert summary["internal_stop_count"] == 1
    assert "species_id\tinput_records\tproteins_written" in summary_tsv.read_text(encoding="utf-8")
    assert "ptae\t5\t3" in summary_tsv.read_text(encoding="utf-8")
    assert "Wrote 3 proteins" in log_path.read_text(encoding="utf-8")


def test_clean_fasta_reads_gzipped_input(tmp_path):
    cleaner = load_script("clean_fasta_headers.py")
    input_fasta = tmp_path / "raw.faa.gz"
    output_fasta = tmp_path / "processed.faa"
    with gzip.open(input_fasta, "wt", encoding="utf-8") as handle:
        handle.write(">gene1\nMKTAA\n")

    cleaner.clean_fasta(input_fasta, output_fasta, species_id="pabi", min_length=5)

    assert output_fasta.read_text(encoding="utf-8") == ">pabi|gene1\nMKTAA\n"


def test_summarize_proteomes_reports_lengths_and_n50(tmp_path):
    summarizer = load_script("summarize_proteomes.py")
    proteome_dir = tmp_path / "proteomes"
    proteome_dir.mkdir()
    output = tmp_path / "species_summary.tsv"
    log_path = tmp_path / "summary.log"
    (proteome_dir / "ptae.faa").write_text(
        ">ptae|gene1\nAAAA\n>ptae|gene2\nAAAAAA\n>ptae|gene3\nAAAAAAAAAA\n",
        encoding="utf-8",
    )
    (proteome_dir / "pmas.faa").write_text(">pmas|gene1\nAAA\n", encoding="utf-8")

    rows = summarizer.summarize_directory(proteome_dir, output, log_path)

    assert rows == [
        {
            "species_id": "pmas",
            "protein_count": 1,
            "median_length": 3,
            "protein_n50": 3,
            "missing_sequences": 0,
            "duplicated_ids": 0,
        },
        {
            "species_id": "ptae",
            "protein_count": 3,
            "median_length": 6,
            "protein_n50": 10,
            "missing_sequences": 0,
            "duplicated_ids": 0,
        },
    ]
    assert output.read_text(encoding="utf-8").splitlines() == [
        "species_id\tprotein_count\tmedian_length\tprotein_n50\tmissing_sequences\tduplicated_ids",
        "pmas\t1\t3\t3\t0\t0",
        "ptae\t3\t6\t10\t0\t0",
    ]
    assert "Summarized 2 proteomes" in log_path.read_text(encoding="utf-8")


def test_summarize_proteomes_can_restrict_to_explicit_inputs(tmp_path):
    summarizer = load_script("summarize_proteomes.py")
    proteome_dir = tmp_path / "proteomes"
    proteome_dir.mkdir()
    output = tmp_path / "species_summary.tsv"
    log_path = tmp_path / "summary.log"
    active_fasta = proteome_dir / "ptae.faa"
    active_fasta.write_text(">ptae|gene1\nAAAA\n", encoding="utf-8")
    (proteome_dir / "holdout.faa").write_text(">holdout|gene1\nAAAAAA\n", encoding="utf-8")

    rows = summarizer.summarize_fastas([active_fasta], output, log_path)

    assert [row["species_id"] for row in rows] == ["ptae"]
    assert output.read_text(encoding="utf-8").splitlines() == [
        "species_id\tprotein_count\tmedian_length\tprotein_n50\tmissing_sequences\tduplicated_ids",
        "ptae\t1\t4\t4\t0\t0",
    ]
    assert "Summarized 1 proteomes" in log_path.read_text(encoding="utf-8")
