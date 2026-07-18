import importlib.util
from pathlib import Path


def load_validator():
    module_path = Path(__file__).resolve().parents[1] / "scripts" / "validate_metadata.py"
    spec = importlib.util.spec_from_file_location("validate_metadata", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_tsv(path, header, rows):
    path.write_text(
        "\t".join(header) + "\n" + "\n".join("\t".join(row) for row in rows) + "\n",
        encoding="utf-8",
    )


def test_validate_species_metadata_accepts_required_columns(tmp_path):
    validator = load_validator()
    path = tmp_path / "species_metadata.tsv"
    write_tsv(
        path,
        [
            "species_id",
            "species_name",
            "common_name",
            "region",
            "group",
            "data_type",
            "protein_fasta",
            "gff3",
            "genome_fasta",
            "source",
            "accession",
            "version",
            "notes",
        ],
        [["ptae", "Pinus_taeda", "Loblolly_pine", "North_America", "tolerant", "proteome", "NA", "NA", "NA", "NA", "NA", "NA", "ok"]],
    )

    assert validator.validate_species_metadata(path) == []


def test_validate_species_metadata_reports_duplicate_ids(tmp_path):
    validator = load_validator()
    path = tmp_path / "species_metadata.tsv"
    write_tsv(
        path,
        [
            "species_id",
            "species_name",
            "common_name",
            "region",
            "group",
            "data_type",
            "protein_fasta",
            "gff3",
            "genome_fasta",
            "source",
            "accession",
            "version",
            "notes",
        ],
        [
            ["ptae", "Pinus_taeda", "Loblolly_pine", "North_America", "tolerant", "proteome", "NA", "NA", "NA", "NA", "NA", "NA", "ok"],
            ["ptae", "Pinus_taeda", "Loblolly_pine", "North_America", "tolerant", "proteome", "NA", "NA", "NA", "NA", "NA", "NA", "duplicate"],
        ],
    )

    errors = validator.validate_species_metadata(path)

    assert any("Duplicate species_id" in error for error in errors)


def test_validate_species_metadata_strict_paths_reports_missing_fasta(tmp_path):
    validator = load_validator()
    path = tmp_path / "species_metadata.tsv"
    write_tsv(
        path,
        [
            "species_id",
            "species_name",
            "common_name",
            "region",
            "group",
            "data_type",
            "protein_fasta",
            "gff3",
            "genome_fasta",
            "source",
            "accession",
            "version",
            "notes",
        ],
        [["ptae", "Pinus_taeda", "Loblolly_pine", "North_America", "tolerant", "proteome", "data/raw/proteomes/missing.faa", "NA", "NA", "NA", "NA", "NA", "missing"]],
    )

    errors = validator.validate_species_metadata(path, strict_paths=True)

    assert any("Missing protein_fasta file" in error for error in errors)


def test_validate_species_metadata_strict_paths_allows_na_placeholder(tmp_path):
    validator = load_validator()
    path = tmp_path / "species_metadata.tsv"
    write_tsv(
        path,
        [
            "species_id",
            "species_name",
            "common_name",
            "region",
            "group",
            "data_type",
            "protein_fasta",
            "gff3",
            "genome_fasta",
            "source",
            "accession",
            "version",
            "notes",
        ],
        [["ptae", "Pinus_taeda", "Loblolly_pine", "North_America", "tolerant", "proteome", "NA", "NA", "NA", "NA", "NA", "NA", "template"]],
    )

    assert validator.validate_species_metadata(path, strict_paths=True) == []


def test_validate_phenotype_matrix_reports_invalid_response_score(tmp_path):
    validator = load_validator()
    path = tmp_path / "phenotype_matrix.tsv"
    write_tsv(
        path,
        [
            "species_id",
            "species_name",
            "pwd_response",
            "response_score",
            "evidence_type",
            "evidence_strength",
            "key_reference",
            "notes",
        ],
        [["ptae", "Pinus_taeda", "tolerant", "9", "inoculation_pathology", "high", "NA", "bad score"]],
    )

    errors = validator.validate_phenotype_matrix(path)

    assert any("Invalid response_score" in error for error in errors)


def test_validate_from_config_works_without_pyyaml(tmp_path, monkeypatch):
    validator = load_validator()
    monkeypatch.setattr(validator, "yaml", None)

    config_dir = tmp_path / "config"
    config_dir.mkdir()
    species_path = config_dir / "species_metadata.tsv"
    phenotype_path = config_dir / "phenotype_matrix.tsv"
    config_path = config_dir / "config.yaml"

    write_tsv(
        species_path,
        [
            "species_id",
            "species_name",
            "common_name",
            "region",
            "group",
            "data_type",
            "protein_fasta",
            "gff3",
            "genome_fasta",
            "source",
            "accession",
            "version",
            "notes",
        ],
        [["ptae", "Pinus_taeda", "Loblolly_pine", "North_America", "tolerant", "proteome", "NA", "NA", "NA", "NA", "NA", "NA", "ok"]],
    )
    write_tsv(
        phenotype_path,
        [
            "species_id",
            "species_name",
            "pwd_response",
            "response_score",
            "evidence_type",
            "evidence_strength",
            "key_reference",
            "notes",
        ],
        [["ptae", "Pinus_taeda", "tolerant", "2", "inoculation_pathology", "high", "NA", "ok"]],
    )
    config_path.write_text(
        f"metadata:\n  species: {species_path.as_posix()}\n  phenotype: {phenotype_path.as_posix()}\n",
        encoding="utf-8",
    )

    assert validator.validate_from_config(config_path) == []


def test_validate_expansion_candidates_requires_columns_and_paths(tmp_path):
    validator = load_validator()
    fasta = tmp_path / "psme.faa"
    fasta.write_text(">gene1\nMPEPTIDE\n", encoding="utf-8")
    path = tmp_path / "expansion_candidate_species.tsv"
    write_tsv(
        path,
        [
            "species_id",
            "species_name",
            "common_name",
            "region",
            "phylogenetic_group",
            "role",
            "recommended_use",
            "protein_fasta",
            "source",
            "accession",
            "version",
            "download_status",
            "notes",
        ],
        [["psme", "Pseudotsuga_menziesii", "Douglas_fir", "North_America", "conifer", "North_American_expansion", "outgroup_or_sensitivity_analysis", str(fasta), "TreeGenes", "NA", "v1", "downloaded", "ok"]],
    )

    assert validator.validate_expansion_candidates(path, strict_paths=True) == []


def test_selected_expansion_species_excludes_distant_and_not_ready(tmp_path):
    validator = load_validator()
    path = tmp_path / "expansion_candidate_species.tsv"
    write_tsv(
        path,
        [
            "species_id",
            "species_name",
            "common_name",
            "region",
            "phylogenetic_group",
            "role",
            "recommended_use",
            "protein_fasta",
            "source",
            "accession",
            "version",
            "download_status",
            "notes",
        ],
        [
            ["psme", "Pseudotsuga_menziesii", "Douglas_fir", "North_America", "conifer", "North_American_expansion", "outgroup_or_sensitivity_analysis", "psme.faa", "TreeGenes", "NA", "v1", "downloaded", "ok"],
            ["gibi", "Ginkgo_biloba", "Ginkgo", "China_mainland", "gymnosperm", "China_mainland_expansion", "distant_gymnosperm_outgroup", "gibi.faa", "TreeGenes", "NA", "v1", "downloaded", "not default"],
            ["prad", "Pinus_radiata", "Monterey_pine", "North_America", "pine", "North_American_candidate", "not_ready_for_orthology", "NA", "NCBI", "NA", "NA", "reviewed_not_downloaded", "not ready"],
        ],
    )

    selected = validator.selected_expansion_species(path)

    assert [row["species_id"] for row in selected] == ["psme"]
