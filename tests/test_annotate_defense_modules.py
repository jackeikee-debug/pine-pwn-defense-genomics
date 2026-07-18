import importlib.util
from pathlib import Path


def load_script():
    module_path = Path(__file__).resolve().parents[1] / "scripts" / "annotate_defense_modules.py"
    spec = importlib.util.spec_from_file_location("annotate_defense_modules", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_annotate_defense_modules_matches_keywords_and_aggregates_by_species(tmp_path: Path) -> None:
    annotator = load_script()
    modules = tmp_path / "defense_modules.yaml"
    stable_genes = tmp_path / "stable_orthogroup_genes.tsv"
    annotations = tmp_path / "annotations.tsv"
    matrix = tmp_path / "matrix.tsv"
    gene_hits = tmp_path / "gene_hits.tsv"

    modules.write_text(
        "modules:\n"
        "  resin_terpenoid:\n"
        "    description: Resin defense\n"
        "    keywords:\n"
        "      - terpene synthase\n"
        "      - cytochrome P450\n"
        "  ros_detoxification:\n"
        "    description: ROS defense\n"
        "    keywords:\n"
        "      - peroxidase\n",
        encoding="utf-8",
    )
    stable_genes.write_text(
        "orthogroup_id\tspecies_id\tgene_id\n"
        "OG1\tptae\tptae|g1\n"
        "OG1\tptae\tptae|g2\n"
        "OG1\tpmas\tpmas|g3\n"
        "OG2\tptae\tptae|g4\n",
        encoding="utf-8",
    )
    annotations.write_text(
        "gene_id\tannotation_source\tannotation_text\n"
        "ptae|g1\teggNOG\tputative terpene synthase family protein\n"
        "ptae|g2\tInterPro\tclass III peroxidase domain\n"
        "pmas|g3\tSwissProt\tcytochrome P450 monooxygenase\n"
        "ptae|g4\tSwissProt\tunknown protein\n",
        encoding="utf-8",
    )

    matrix_rows, hit_rows = annotator.annotate_defense_modules(
        stable_genes_path=stable_genes,
        annotations_path=annotations,
        modules_path=modules,
        matrix_output=matrix,
        gene_hits_output=gene_hits,
    )

    assert matrix_rows == [
        {
            "orthogroup_id": "OG1",
            "module_id": "resin_terpenoid",
            "species_id": "pmas",
            "gene_count": "1",
            "gene_ids": "pmas|g3",
            "matched_keywords": "cytochrome P450",
            "annotation_sources": "SwissProt",
            "evidence_type": "keyword_annotation_match",
        },
        {
            "orthogroup_id": "OG1",
            "module_id": "resin_terpenoid",
            "species_id": "ptae",
            "gene_count": "1",
            "gene_ids": "ptae|g1",
            "matched_keywords": "terpene synthase",
            "annotation_sources": "eggNOG",
            "evidence_type": "keyword_annotation_match",
        },
        {
            "orthogroup_id": "OG1",
            "module_id": "ros_detoxification",
            "species_id": "ptae",
            "gene_count": "1",
            "gene_ids": "ptae|g2",
            "matched_keywords": "peroxidase",
            "annotation_sources": "InterPro",
            "evidence_type": "keyword_annotation_match",
        },
    ]
    assert hit_rows[0]["annotation_text"] == "putative terpene synthase family protein"
    assert matrix.read_text(encoding="utf-8").splitlines()[0].startswith(
        "orthogroup_id\tmodule_id\tspecies_id\tgene_count"
    )
    assert gene_hits.read_text(encoding="utf-8").splitlines()[0] == (
        "orthogroup_id\tspecies_id\tgene_id\tmodule_id\tmatched_keyword\t"
        "annotation_source\tannotation_text\tevidence_type"
    )


def test_annotate_defense_modules_does_not_match_uppercase_acronyms_inside_lowercase_words(
    tmp_path: Path,
) -> None:
    annotator = load_script()
    modules = tmp_path / "defense_modules.yaml"
    stable_genes = tmp_path / "stable_orthogroup_genes.tsv"
    annotations = tmp_path / "annotations.tsv"

    modules.write_text(
        "modules:\n"
        "  hormone_signaling:\n"
        "    description: Hormone signaling\n"
        "    keywords:\n"
        "      - EIN\n",
        encoding="utf-8",
    )
    stable_genes.write_text(
        "orthogroup_id\tspecies_id\tgene_id\n"
        "OG1\tptae\tptae|g1\n"
        "OG1\tptae\tptae|g2\n",
        encoding="utf-8",
    )
    annotations.write_text(
        "gene_id\tannotation_source\tannotation_text\n"
        "ptae|g1\tSwissProt\tuncharacterized protein\n"
        "ptae|g2\tSwissProt\tethylene insensitive protein EIN2\n",
        encoding="utf-8",
    )

    matrix_rows, hit_rows = annotator.annotate_defense_modules(
        stable_genes_path=stable_genes,
        annotations_path=annotations,
        modules_path=modules,
        matrix_output=tmp_path / "matrix.tsv",
        gene_hits_output=tmp_path / "gene_hits.tsv",
    )

    assert [row["gene_id"] for row in hit_rows] == ["ptae|g2"]
    assert matrix_rows[0]["gene_ids"] == "ptae|g2"
