import importlib.util
from pathlib import Path


def load_script():
    module_path = Path(__file__).resolve().parents[1] / "scripts" / "map_effectors_to_host_modules.py"
    spec = importlib.util.spec_from_file_location("map_effectors_to_host_modules", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_map_effectors_to_host_modules_links_medium_effectors_to_host_shortlist(tmp_path: Path) -> None:
    mapper = load_script()
    effectors = tmp_path / "effectors.tsv"
    mapping = tmp_path / "mapping.tsv"
    host = tmp_path / "host.tsv"
    output = tmp_path / "links.tsv"

    effectors.write_text(
        "protein_id\teffector_class\tevidence\tconfidence\tannotation_text\tlength\tcysteine_fraction\n"
        "BX1\tcell_wall_modifying\tkeyword:pectate lyase\tmedium\tPectate lyase\t220\t0.0200\n"
        "BX2\tcell_wall_modifying\tkeyword:cellulase\tmedium\tCellulase\t210\t0.0100\n"
        "BX3\tsmall_secreted\tsecretome:small_secreted\tlow\tunnamed\t100\t0.0100\n"
        "BX4\tprotease\tkeyword:peptidase\tmedium\tPeptidase\t300\t0.0300\n",
        encoding="utf-8",
    )
    mapping.write_text(
        "effector_class\thost_module\tlink_type\tevidence_basis\tnotes\n"
        "cell_wall_modifying\twound_periderm\tputative_functional_link\tkeyword_mapping\twall response\n"
        "protease\timmune_signaling\tputative_functional_link\tkeyword_mapping\tdefense signaling\n",
        encoding="utf-8",
    )
    host.write_text(
        "shortlist_rank\tshortlist_group\tshortlist_group_rank\tshortlist_reason\tpriority_rank\tpriority_tier\t"
        "priority_score\tmodule_id\torthogroup_id\tcopy_bias_direction\tspecies_count\tgene_count\t"
        "susceptible_species_count\ttolerant_species_count\toutgroup_species_count\tsusceptible_gene_count\t"
        "tolerant_gene_count\toutgroup_gene_count\tmax_species_gene_count\tmin_species_gene_count\t"
        "copy_number_range\tmean_pident\tmean_bitscore\tbest_bitscore\tmatched_keywords\tannotation_sources\t"
        "evidence_types\tbest_expansion_orthogroup\tpilot_core_retention\tcore_jaccard\n"
        "1\twound_periderm|susceptible_enriched\t1\treason\t10\thigh\t11\twound_periderm\tOG1\t"
        "susceptible_enriched\t4\t12\t2\t1\t1\t10\t1\t1\t8\t1\t7\t60.0\t200.0\t260.0\tlaccase\tSwissProt\tkeyword\tOGE1\t1.000\t0.950\n"
        "2\timmune_signaling|tolerant_enriched\t1\treason\t11\thigh\t10\timmune_signaling\tOG2\t"
        "tolerant_enriched\t4\t8\t1\t2\t1\t1\t6\t1\t5\t1\t4\t55.0\t180.0\t240.0\tMYB\tSwissProt\tkeyword\tOGE2\t0.980\t0.900\n",
        encoding="utf-8",
    )

    rows = mapper.map_effectors_to_host_modules(
        effectors_path=effectors,
        mapping_path=mapping,
        host_shortlist_path=host,
        output_path=output,
        min_confidence="medium",
    )

    assert len(rows) == 2
    assert rows[0]["effector_class"] == "cell_wall_modifying"
    assert rows[0]["host_module"] == "wound_periderm"
    assert rows[0]["effector_count"] == "2"
    assert rows[0]["host_candidate_count"] == "1"
    assert rows[0]["effector_ids"] == "BX1;BX2"
    assert rows[0]["host_orthogroups"] == "OG1"
    assert rows[0]["link_confidence"] == "candidate"
    assert rows[1]["effector_class"] == "protease"
    assert rows[1]["host_module"] == "immune_signaling"
    assert output.read_text(encoding="utf-8").splitlines()[0].startswith(
        "effector_class\thost_module\tlink_type"
    )
