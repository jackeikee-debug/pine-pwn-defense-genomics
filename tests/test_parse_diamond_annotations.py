import importlib.util
from pathlib import Path


def load_script():
    module_path = Path(__file__).resolve().parents[1] / "scripts" / "parse_diamond_annotations.py"
    spec = importlib.util.spec_from_file_location("parse_diamond_annotations", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_parse_diamond_annotations_keeps_best_hit_and_writes_annotation_table(tmp_path: Path) -> None:
    parser = load_script()
    diamond = tmp_path / "diamond.tsv"
    annotations = tmp_path / "annotations.tsv"
    hits = tmp_path / "hits.tsv"

    diamond.write_text(
        "ptae|g1\tsp|A0A000|TPS_PINUS\t63.2\t401\t1e-80\t280\tTerpene synthase OS=Pinus taeda\n"
        "ptae|g1\tsp|B0B000|LOW_PINUS\t51.0\t301\t1e-20\t110\tLower scoring hit OS=Pinus\n"
        "pmas|g2\tsp|C0C000|PERO_ARATH\t42.5\t220\t2e-30\t140\tPeroxidase 12 OS=Arabidopsis thaliana\n",
        encoding="utf-8",
    )

    annotation_rows, hit_rows = parser.parse_diamond_annotations(
        diamond_path=diamond,
        annotation_output=annotations,
        hits_output=hits,
        min_pident=40.0,
        min_bitscore=100.0,
        source_name="SwissProt_DIAMOND",
    )

    assert annotation_rows == [
        {
            "gene_id": "ptae|g1",
            "annotation_source": "SwissProt_DIAMOND",
            "annotation_text": "Terpene synthase OS=Pinus taeda",
        },
        {
            "gene_id": "pmas|g2",
            "annotation_source": "SwissProt_DIAMOND",
            "annotation_text": "Peroxidase 12 OS=Arabidopsis thaliana",
        },
    ]
    assert hit_rows[0]["subject_id"] == "sp|A0A000|TPS_PINUS"
    assert hit_rows[0]["bitscore"] == "280"
    assert annotations.read_text(encoding="utf-8").splitlines()[0] == (
        "gene_id\tannotation_source\tannotation_text"
    )
