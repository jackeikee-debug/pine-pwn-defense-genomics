import csv
import glob
import sys


def stable_pine_proteomes(_wildcards):
    with open("config/species_metadata.tsv", encoding="utf-8", newline="") as handle:
        return [
            row["protein_fasta"]
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("protein_fasta") not in {"", "NA"}
            and row.get("group") != "outgroup"
        ]


rule build_manuscript_rnaseq_qc_audit:
    input:
        provenance="config/rnaseq_quantification_runs.tsv",
        pmas=lambda wildcards: sorted(glob.glob("data/external/galaxy/effector_target_network_pmas_1dpi_fastqc/*.fastqc.txt")),
        pden=lambda wildcards: sorted(glob.glob("data/external/galaxy/effector_target_network_pden_full_fastqc/*.fastqc.txt"))
    output:
        runs="results/tables/manuscript_rnaseq_fastqc_audit.tsv",
        summary="results/tables/manuscript_rnaseq_qc_summary.tsv",
        report="results/reports/manuscript_rnaseq_qc_audit.md"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/manuscript_rnaseq_qc_audit.log"
    benchmark:
        "results/benchmarks/manuscript_rnaseq_qc_audit.tsv"
    shell:
        "python scripts/build_rnaseq_qc_audit.py "
        "--pmas-fastqc-dir data/external/galaxy/effector_target_network_pmas_1dpi_fastqc "
        "--pden-fastqc-dir data/external/galaxy/effector_target_network_pden_full_fastqc "
        "--provenance {input.provenance} "
        "--run-output {output.runs} "
        "--summary-output {output.summary} "
        "--report {output.report} "
        "--log {log}"


rule build_candidate_coevolution_genes:
    input:
        links="results/tables/effector_host_module_links.tsv",
        host_shortlist="results/tables/defense_module_candidate_shortlist.tsv",
        effectors="results/tables/nematode_candidate_effectors.tsv"
    output:
        "results/tables/candidate_coevolution_genes.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/candidate_coevolution_genes.log"
    shell:
        "python scripts/build_candidate_coevolution_genes.py "
        "--links {input.links} "
        "--host-shortlist {input.host_shortlist} "
        "--effectors {input.effectors} "
        "--output {output} > {log} 2>&1"


rule contrast_candidate_regions:
    input:
        candidates="results/tables/candidate_coevolution_genes.tsv",
        stable_orthogroups="results/tables/stable_orthogroups.tsv",
        species_metadata=config["metadata"]["species"]
    output:
        "results/tables/east_asia_vs_north_america_candidate_contrast.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/east_asia_vs_north_america_candidate_contrast.log"
    shell:
        "python scripts/contrast_candidate_regions.py "
        "--candidates {input.candidates} "
        "--stable-orthogroups {input.stable_orthogroups} "
        "--species-metadata {input.species_metadata} "
        "--output {output} "
        "--ratio-threshold 1.5 "
        "--pseudocount 0.5 > {log} 2>&1"


rule select_top_regional_contrasts:
    input:
        "results/tables/east_asia_vs_north_america_candidate_contrast.tsv"
    output:
        associations="results/tables/regional_top_contrasts.tsv",
        orthogroups="results/tables/regional_top_orthogroups.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/regional_top_candidates.log"
    shell:
        "python scripts/select_top_regional_contrasts.py "
        "--contrast {input} "
        "--association-output {output.associations} "
        "--orthogroup-output {output.orthogroups} "
        "--max-per-direction-module-effector 5 > {log} 2>&1"


rule build_regional_family_review:
    input:
        top_orthogroups="results/tables/regional_top_orthogroups.tsv",
        stable_genes="results/tables/stable_orthogroup_genes.tsv",
        species_metadata=config["metadata"]["species"],
        proteomes=expand("data/processed/proteomes_longest/{species_id}.faa", species_id=PROTEOME_SPECIES)
    output:
        review_dir=directory("results/reports/regional_family_review"),
        summary="results/tables/regional_family_review_summary.tsv",
        manifest="results/tables/regional_family_review_manifest.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/regional_family_review.log"
    shell:
        "python scripts/build_regional_family_review.py "
        "--top-orthogroups {input.top_orthogroups} "
        "--stable-genes {input.stable_genes} "
        "--species-metadata {input.species_metadata} "
        "--proteome-dir data/processed/proteomes_longest "
        "--output-dir {output.review_dir} "
        "--summary {output.summary} "
        "--manifest {output.manifest} > {log} 2>&1"


rule run_regional_family_validation:
    input:
        review_summary="results/tables/regional_family_review_summary.tsv",
        review_manifest="results/tables/regional_family_review_manifest.tsv"
    output:
        alignment_dir=directory("results/reports/regional_family_validation/alignments"),
        tree_dir=directory("results/reports/regional_family_validation/trees"),
        summary="results/tables/regional_family_validation_summary.tsv"
    conda:
        "../envs/phylo.yaml"
    log:
        "results/logs/regional_family_validation.log"
    shell:
        "python scripts/run_family_validation.py "
        "--review-summary {input.review_summary} "
        "--alignment-dir {output.alignment_dir} "
        "--tree-dir {output.tree_dir} "
        "--summary {output.summary} > {log} 2>&1"


rule rank_candidate_evidence:
    input:
        regional_top="results/tables/regional_top_orthogroups.tsv",
        validation_summary="results/tables/regional_family_validation_summary.tsv",
        cafe_intersections="results/tables/cafe5_top_candidate_intersections.tsv"
    output:
        ranking="results/tables/candidate_evidence_ranking.tsv",
        strong="results/tables/strong_candidate_orthogroups.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/candidate_evidence_ranking.log"
    shell:
        "python scripts/rank_candidate_evidence.py "
        "--regional-top {input.regional_top} "
        "--validation-summary {input.validation_summary} "
        "--cafe-intersections {input.cafe_intersections} "
        "--output {output.ranking} "
        "--strong-output {output.strong} > {log} 2>&1"


rule build_model_supported_candidate_review:
    input:
        ranking="results/tables/candidate_evidence_ranking.tsv"
    output:
        table="results/tables/model_supported_candidate_review.tsv",
        markdown="results/reports/model_supported_candidate_review.md"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/model_supported_candidate_review.log"
    shell:
        "python scripts/build_model_supported_review.py "
        "--ranking {input.ranking} "
        "--output {output.table} "
        "--markdown {output.markdown} > {log} 2>&1"


rule build_candidate_manual_checklist:
    input:
        review="results/tables/model_supported_candidate_review.tsv"
    output:
        "results/tables/model_supported_candidate_manual_checklist.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/model_supported_candidate_manual_checklist.log"
    shell:
        "python scripts/build_candidate_manual_checklist.py "
        "--review {input.review} "
        "--output {output} > {log} 2>&1"


rule build_candidate_tree_metrics:
    input:
        checklist="results/tables/model_supported_candidate_manual_checklist.tsv"
    output:
        "results/tables/model_supported_candidate_tree_metrics.tsv"
    conda:
        "../envs/base.yaml"
    params:
        python=sys.executable
    log:
        "results/logs/model_supported_candidate_tree_metrics.log"
    shell:
        "\"{params.python}\" scripts/build_candidate_tree_metrics.py "
        "--checklist {input.checklist} "
        "--output {output} > {log} 2>&1"


rule build_candidate_tree_review_report:
    input:
        checklist="results/tables/model_supported_candidate_manual_checklist.tsv",
        metrics="results/tables/model_supported_candidate_tree_metrics.tsv"
    output:
        markdown="results/reports/model_supported_candidate_tree_review.md",
        ascii_dir=directory("results/reports/model_supported_candidate_tree_review_ascii")
    conda:
        "../envs/base.yaml"
    params:
        python=sys.executable
    log:
        "results/logs/model_supported_candidate_tree_review.log"
    shell:
        "\"{params.python}\" scripts/build_candidate_tree_review_report.py "
        "--checklist {input.checklist} "
        "--metrics {input.metrics} "
        "--markdown {output.markdown} "
        "--ascii-dir {output.ascii_dir} > {log} 2>&1"


rule build_candidate_review_decisions:
    input:
        checklist="results/tables/model_supported_candidate_manual_checklist.tsv",
        metrics="results/tables/model_supported_candidate_tree_metrics.tsv"
    output:
        table="results/tables/model_supported_candidate_review_decisions.tsv",
        markdown="results/reports/model_supported_candidate_review_decisions.md"
    conda:
        "../envs/base.yaml"
    params:
        python=sys.executable
    log:
        "results/logs/model_supported_candidate_review_decisions.log"
    shell:
        "\"{params.python}\" scripts/build_candidate_review_decisions.py "
        "--checklist {input.checklist} "
        "--metrics {input.metrics} "
        "--output {output.table} "
        "--markdown {output.markdown} > {log} 2>&1"


rule build_manuscript_candidate_summary:
    input:
        review="results/tables/model_supported_candidate_review.tsv",
        decisions="results/tables/model_supported_candidate_review_decisions.tsv"
    output:
        table="results/tables/manuscript_candidate_summary.tsv",
        markdown="results/reports/manuscript_candidate_summary.md"
    conda:
        "../envs/base.yaml"
    params:
        python=sys.executable
    log:
        "results/logs/manuscript_candidate_summary.log"
    shell:
        "\"{params.python}\" scripts/build_manuscript_candidate_summary.py "
        "--review {input.review} "
        "--decisions {input.decisions} "
        "--output {output.table} "
        "--markdown {output.markdown} > {log} 2>&1"


rule build_mechanism_axis_summary:
    input:
        candidates="results/tables/manuscript_candidate_summary.tsv",
        links="results/tables/effector_host_module_links.tsv",
        axes=config["metadata"]["mechanism_axes"]
    output:
        axis="results/tables/mechanism_axis_summary.tsv",
        gaps="results/tables/mechanism_evidence_gaps.tsv",
        markdown="results/reports/mechanism_strategy.md"
    conda:
        "../envs/base.yaml"
    params:
        python=sys.executable
    log:
        "results/logs/mechanism_axis_summary.log"
    shell:
        "\"{params.python}\" scripts/build_mechanism_axis_summary.py "
        "--candidates {input.candidates} "
        "--links {input.links} "
        "--axes {input.axes} "
        "--axis-output {output.axis} "
        "--gaps-output {output.gaps} "
        "--markdown {output.markdown} > {log} 2>&1"


rule build_public_expression_evidence_plan:
    input:
        gaps="results/tables/mechanism_evidence_gaps.tsv",
        sources=config["metadata"]["public_expression_evidence_sources"]
    output:
        table="results/tables/public_expression_evidence_plan.tsv",
        markdown="results/reports/public_expression_evidence_plan.md"
    conda:
        "../envs/base.yaml"
    params:
        python=sys.executable
    log:
        "results/logs/public_expression_evidence_plan.log"
    shell:
        "\"{params.python}\" scripts/build_public_expression_evidence_plan.py "
        "--gaps {input.gaps} "
        "--sources {input.sources} "
        "--output {output.table} "
        "--markdown {output.markdown} > {log} 2>&1"


rule fetch_public_expression_texts:
    input:
        sources=config["metadata"]["public_expression_evidence_sources"]
    output:
        text_dir=directory(config["metadata"]["public_expression_text_dir"]),
        manifest="results/tables/public_expression_text_download_manifest.tsv"
    conda:
        "../envs/base.yaml"
    params:
        python=sys.executable
    log:
        "results/logs/public_expression_text_download.log"
    shell:
        "\"{params.python}\" scripts/fetch_public_expression_texts.py "
        "--sources {input.sources} "
        "--output-dir {output.text_dir} "
        "--manifest {output.manifest} > {log} 2>&1"


rule extract_public_expression_literature_evidence:
    input:
        candidates="results/tables/manuscript_candidate_summary.tsv",
        evidence_plan="results/tables/public_expression_evidence_plan.tsv",
        text_dir=config["metadata"]["public_expression_text_dir"]
    output:
        table="results/tables/public_expression_literature_evidence.tsv",
        markdown="results/reports/public_expression_literature_evidence.md"
    conda:
        "../envs/base.yaml"
    params:
        python=sys.executable
    log:
        "results/logs/public_expression_literature_evidence.log"
    shell:
        "\"{params.python}\" scripts/extract_public_expression_literature_evidence.py "
        "--candidates {input.candidates} "
        "--evidence-plan {input.evidence_plan} "
        "--text-dir {input.text_dir} "
        "--output {output.table} "
        "--markdown {output.markdown} > {log} 2>&1"


rule fetch_public_expression_supplements:
    input:
        sources=config["metadata"]["public_expression_supplement_sources"]
    output:
        supplement_dir=directory(config["metadata"]["public_expression_supplement_dir"]),
        manifest="results/tables/public_expression_supplement_download_manifest.tsv"
    conda:
        "../envs/base.yaml"
    params:
        python=sys.executable
    log:
        "results/logs/public_expression_supplement_download.log"
    shell:
        "\"{params.python}\" scripts/fetch_public_expression_supplements.py "
        "--sources {input.sources} "
        "--output-dir {output.supplement_dir} "
        "--manifest {output.manifest} > {log} 2>&1"


rule inventory_public_expression_supplements:
    input:
        manifest="results/tables/public_expression_supplement_download_manifest.tsv"
    output:
        table="results/tables/public_expression_supplement_inventory.tsv",
        markdown="results/reports/public_expression_supplement_inventory.md"
    conda:
        "../envs/base.yaml"
    params:
        python=sys.executable
    log:
        "results/logs/public_expression_supplement_inventory.log"
    shell:
        "\"{params.python}\" scripts/inventory_public_expression_supplements.py "
        "--manifest {input.manifest} "
        "--output {output.table} "
        "--markdown {output.markdown} > {log} 2>&1"


rule build_effector_target_network_expression_support:
    input:
        mechanism_axis="results/tables/mechanism_axis_summary.tsv",
        candidates="results/tables/manuscript_candidate_summary.tsv",
        literature="results/tables/public_expression_literature_evidence.tsv",
        stable_genes="results/tables/stable_orthogroup_genes.tsv",
        annotations=config["metadata"]["protein_function_annotations"],
        supplement_manifest="results/tables/public_expression_supplement_download_manifest.tsv"
    output:
        table="results/tables/effector_target_network_expression_support.tsv",
        markdown="results/reports/effector_target_network_expression_support.md"
    conda:
        "../envs/base.yaml"
    params:
        python=sys.executable
    log:
        "results/logs/effector_target_network_expression_support.log"
    shell:
        "\"{params.python}\" scripts/build_effector_target_network_expression_support.py "
        "--mechanism-axis {input.mechanism_axis} "
        "--candidates {input.candidates} "
        "--literature {input.literature} "
        "--stable-genes {input.stable_genes} "
        "--annotations {input.annotations} "
        "--supplement-manifest {input.supplement_manifest} "
        "--output {output.table} "
        "--markdown {output.markdown} > {log} 2>&1"


rule extract_effector_target_network_supplement_gene_hits:
    input:
        effector_target_network="results/tables/effector_target_network_expression_support.tsv",
        manifest="results/tables/public_expression_supplement_download_manifest.tsv"
    output:
        table="results/tables/effector_target_network_supplement_gene_hits.tsv",
        markdown="results/reports/effector_target_network_supplement_gene_hits.md"
    conda:
        "../envs/base.yaml"
    params:
        python=sys.executable
    log:
        "results/logs/effector_target_network_supplement_gene_hits.log"
    shell:
        "\"{params.python}\" scripts/extract_effector_target_network_supplement_gene_hits.py "
        "--effector_target_network {input.effector_target_network} "
        "--manifest {input.manifest} "
        "--output {output.table} "
        "--markdown {output.markdown} > {log} 2>&1"


rule build_effector_target_network_identifier_crosswalk_audit:
    input:
        effector_target_network="results/tables/effector_target_network_expression_support.tsv",
        hits="results/tables/effector_target_network_supplement_gene_hits.tsv",
        inventory="results/tables/public_expression_supplement_inventory.tsv"
    output:
        table="results/tables/effector_target_network_identifier_crosswalk_audit.tsv",
        markdown="results/reports/effector_target_network_identifier_crosswalk_audit.md"
    conda:
        "../envs/base.yaml"
    params:
        python=sys.executable
    log:
        "results/logs/effector_target_network_identifier_crosswalk_audit.log"
    shell:
        "\"{params.python}\" scripts/build_effector_target_network_identifier_crosswalk_audit.py "
        "--effector_target_network {input.effector_target_network} "
        "--hits {input.hits} "
        "--inventory {input.inventory} "
        "--output {output.table} "
        "--markdown {output.markdown} > {log} 2>&1"


rule build_effector_target_network_functional_bridge_audit:
    input:
        network_seed="results/tables/effector_target_network_effector_target_network_seed.tsv",
        locus_map="config/effector_target_network_arabidopsis_locus_map.tsv",
        supplement_manifest="results/tables/public_expression_supplement_download_manifest.tsv"
    output:
        table="results/tables/effector_target_network_functional_bridge_audit.tsv",
        markdown="results/reports/effector_target_network_functional_bridge_audit.md"
    conda:
        "../envs/base.yaml"
    benchmark:
        "results/benchmarks/effector_target_network_functional_bridge_audit.txt"
    params:
        python=sys.executable
    log:
        "results/logs/effector_target_network_functional_bridge_audit.log"
    shell:
        "\"{params.python}\" scripts/build_effector_target_network_functional_bridge_audit.py "
        "--network-seed {input.network_seed} "
        "--locus-map {input.locus_map} "
        "--supplement-manifest {input.supplement_manifest} "
        "--output {output.table} "
        "--markdown {output.markdown} > {log} 2>&1"


rule extract_effector_target_network_pden_trinity_sequences:
    input:
        bridge="results/tables/effector_target_network_functional_bridge_audit.tsv",
        oa_package=config["metadata"]["pden_trinity_oa_package"]
    output:
        fasta="data/interim/effector_target_network_pden_trinity_bridge_transcripts.fasta",
        full_reference="data/interim/pden_published_trinity_reference.fasta",
        manifest="results/tables/effector_target_network_pden_trinity_sequence_manifest.tsv",
        markdown="results/reports/effector_target_network_pden_trinity_sequence_extraction.md"
    conda:
        "../envs/base.yaml"
    benchmark:
        "results/benchmarks/effector_target_network_pden_trinity_sequence_extraction.txt"
    params:
        python=sys.executable,
        package_member=config["metadata"]["pden_trinity_package_member"],
        zip_member=config["metadata"]["pden_trinity_zip_member"]
    log:
        "results/logs/effector_target_network_pden_trinity_sequence_extraction.log"
    shell:
        "\"{params.python}\" scripts/extract_effector_target_network_pden_trinity_sequences.py "
        "--bridge {input.bridge} "
        "--oa-package {input.oa_package} "
        "--package-member {params.package_member} "
        "--zip-member {params.zip_member} "
        "--output-fasta {output.fasta} "
        "--full-reference-output {output.full_reference} "
        "--manifest {output.manifest} "
        "--markdown {output.markdown} > {log} 2>&1"


rule build_pden_trinity_tx2gene:
    input:
        fasta="data/interim/pden_published_trinity_reference.fasta"
    output:
        table="data/interim/pden_published_trinity_tx2gene.tsv",
        markdown="results/reports/pden_published_trinity_tx2gene.md"
    conda:
        "../envs/base.yaml"
    benchmark:
        "results/benchmarks/pden_published_trinity_tx2gene.txt"
    params:
        python=sys.executable
    log:
        "results/logs/pden_published_trinity_tx2gene.log"
    shell:
        "\"{params.python}\" scripts/build_trinity_tx2gene.py "
        "--fasta {input.fasta} "
        "--output {output.table} "
        "--markdown {output.markdown} > {log} 2>&1"


rule summarize_effector_target_network_salmon_pilot:
    input:
        candidates="results/tables/effector_target_network_pden_candidate_expression_audit.tsv",
        transcript_quant="data/external/galaxy/effector_target_network_pden_pilot/SRR8061568.salmon.quant.sf",
        gene_quant="data/external/galaxy/effector_target_network_pden_pilot/SRR8061568.salmon.genes.sf"
    output:
        table="results/tables/effector_target_network_pden_salmon_pilot_candidate_abundance.tsv",
        markdown="results/reports/effector_target_network_pden_salmon_pilot_candidate_abundance.md"
    conda:
        "../envs/base.yaml"
    benchmark:
        "results/benchmarks/effector_target_network_pden_salmon_pilot_candidate_abundance.txt"
    params:
        python=sys.executable,
        sample_id=config["galaxy_effector_target_network_pilot"]["run_accession"],
        condition=config["galaxy_effector_target_network_pilot"]["condition"]
    log:
        "results/logs/effector_target_network_pden_salmon_pilot_candidate_abundance.log"
    shell:
        "\"{params.python}\" scripts/summarize_effector_target_network_salmon_pilot.py "
        "--candidates {input.candidates} "
        "--transcript-quant {input.transcript_quant} "
        "--gene-quant {input.gene_quant} "
        "--output {output.table} "
        "--markdown {output.markdown} "
        "--sample-id {params.sample_id} "
        "--condition {params.condition} > {log} 2>&1"


rule build_pden_ha_diamond_db:
    input:
        proteome=config["effector_target_network_expression_bridge"]["pden_current_proteome"]
    output:
        db="data/interim/pden_ha_protein.dmnd"
    conda:
        "../envs/annotation.yaml"
    benchmark:
        "results/benchmarks/pden_ha_diamond_db.txt"
    params:
        db_prefix="data/interim/pden_ha_protein"
    log:
        "results/logs/pden_ha_diamond_db.log"
    shell:
        "diamond makedb --in {input.proteome} -d {params.db_prefix} > {log} 2>&1"


rule blastx_effector_target_network_pden_trinity_to_current_pd:
    input:
        db="data/interim/pden_ha_protein.dmnd",
        fasta="data/interim/effector_target_network_pden_trinity_bridge_transcripts.fasta"
    output:
        "results/tables/effector_target_network_pden_trinity_vs_current_pd_blastx.tsv"
    conda:
        "../envs/annotation.yaml"
    benchmark:
        "results/benchmarks/effector_target_network_pden_trinity_vs_current_pd_blastx.txt"
    threads:
        4
    params:
        db_prefix="data/interim/pden_ha_protein",
        evalue=config["effector_target_network_expression_bridge"]["pden_blastx_evalue"],
        max_target_seqs=config["effector_target_network_expression_bridge"]["pden_blastx_max_target_seqs"]
    log:
        "results/logs/effector_target_network_pden_trinity_vs_current_pd_blastx.log"
    shell:
        "diamond blastx "
        "-d {params.db_prefix} "
        "-q {input.fasta} "
        "-o {output} "
        "--outfmt 6 qseqid sseqid pident length evalue bitscore qlen slen qstart qend sstart send "
        "--max-target-seqs {params.max_target_seqs} "
        "--evalue {params.evalue} "
        "--threads {threads} > {log} 2>&1"


rule build_effector_target_network_pden_trinity_to_pd_mapping:
    input:
        blastx="results/tables/effector_target_network_pden_trinity_vs_current_pd_blastx.tsv",
        sequence_manifest="results/tables/effector_target_network_pden_trinity_sequence_manifest.tsv",
        network_seed="results/tables/effector_target_network_effector_target_network_seed.tsv"
    output:
        table="results/tables/effector_target_network_pden_trinity_to_pd_mapping.tsv",
        markdown="results/reports/effector_target_network_pden_trinity_to_pd_mapping.md"
    conda:
        "../envs/base.yaml"
    benchmark:
        "results/benchmarks/effector_target_network_pden_trinity_to_pd_mapping.txt"
    params:
        python=sys.executable,
        min_bitscore=config["effector_target_network_expression_bridge"]["pden_blastx_min_bitscore"],
        min_pident=config["effector_target_network_expression_bridge"]["pden_blastx_min_pident"]
    log:
        "results/logs/effector_target_network_pden_trinity_to_pd_mapping.log"
    shell:
        "\"{params.python}\" scripts/build_effector_target_network_pden_trinity_to_pd_mapping.py "
        "--blastx {input.blastx} "
        "--sequence-manifest {input.sequence_manifest} "
        "--network-seed {input.network_seed} "
        "--output {output.table} "
        "--markdown {output.markdown} "
        "--min-bitscore {params.min_bitscore} "
        "--min-pident {params.min_pident} > {log} 2>&1"


rule build_effector_target_network_pden_expression_source_audit:
    input:
        oa_package=config["metadata"]["pden_trinity_oa_package"],
        mapping="results/tables/effector_target_network_pden_trinity_to_pd_mapping.tsv"
    output:
        source_table="results/tables/effector_target_network_pden_expression_source_audit.tsv",
        candidate_table="results/tables/effector_target_network_pden_candidate_expression_audit.tsv",
        markdown="results/reports/effector_target_network_pden_expression_source_audit.md"
    conda:
        "../envs/base.yaml"
    benchmark:
        "results/benchmarks/effector_target_network_pden_expression_source_audit.txt"
    params:
        python=sys.executable
    log:
        "results/logs/effector_target_network_pden_expression_source_audit.log"
    shell:
        "\"{params.python}\" scripts/build_effector_target_network_pden_expression_source_audit.py "
        "--oa-package {input.oa_package} "
        "--mapping {input.mapping} "
        "--source-output {output.source_table} "
        "--candidate-output {output.candidate_table} "
        "--markdown {output.markdown} > {log} 2>&1"


rule build_effector_target_network_pden_rnaseq_manifest:
    output:
        manifest="results/tables/effector_target_network_pden_rnaseq_run_manifest.tsv",
        audit="results/tables/effector_target_network_pden_rnaseq_metadata_audit.tsv",
        markdown="results/reports/effector_target_network_pden_rnaseq_metadata.md"
    conda:
        "../envs/base.yaml"
    benchmark:
        "results/benchmarks/effector_target_network_pden_rnaseq_metadata.txt"
    params:
        python=sys.executable,
        accession=config["effector_target_network_pden_rnaseq"]["accession"],
        api_url=config["effector_target_network_pden_rnaseq"]["ena_api_url"],
        tissue=config["effector_target_network_pden_rnaseq"]["tissue"],
        collection_time=config["effector_target_network_pden_rnaseq"]["collection_time"],
        expected_runs=config["effector_target_network_pden_rnaseq"]["expected_runs"],
        expected_replicates=config["effector_target_network_pden_rnaseq"]["expected_replicates_per_condition"],
        expected_bioproject=config["effector_target_network_pden_rnaseq"]["bioproject_accession"],
        timeout_seconds=config["effector_target_network_pden_rnaseq"]["timeout_seconds"]
    log:
        "results/logs/effector_target_network_pden_rnaseq_metadata.log"
    shell:
        "\"{params.python}\" scripts/build_effector_target_network_pden_rnaseq_manifest.py "
        "--accession {params.accession} "
        "--api-url {params.api_url} "
        "--tissue {params.tissue} "
        "--collection-time {params.collection_time} "
        "--expected-runs {params.expected_runs} "
        "--expected-replicates {params.expected_replicates} "
        "--expected-bioproject {params.expected_bioproject} "
        "--expected-secondary-study {params.accession} "
        "--timeout-seconds {params.timeout_seconds} "
        "--output {output.manifest} "
        "--audit-output {output.audit} "
        "--markdown {output.markdown} "
        "--strict > {log} 2>&1"


rule build_effector_target_network_network_seed:
    input:
        effector_target_network="results/tables/effector_target_network_expression_support.tsv",
        links="results/tables/effector_host_module_links.tsv"
    output:
        table="results/tables/effector_target_network_effector_target_network_seed.tsv",
        markdown="results/reports/effector_target_network_effector_target_network_seed.md"
    conda:
        "../envs/base.yaml"
    params:
        python=sys.executable
    log:
        "results/logs/effector_target_network_effector_target_network_seed.log"
    shell:
        "\"{params.python}\" scripts/build_effector_target_network_network_seed.py "
        "--effector_target_network {input.effector_target_network} "
        "--links {input.links} "
        "--output {output.table} "
        "--markdown {output.markdown} > {log} 2>&1"


rule fetch_effector_target_network_string_network:
    input:
        seed="results/tables/effector_target_network_effector_target_network_seed.tsv"
    output:
        mapping="results/tables/effector_target_network_string_seed_mapping.tsv",
        interactions="results/tables/effector_target_network_string_interactions.tsv",
        manifest="results/tables/effector_target_network_string_download_manifest.tsv"
    conda:
        "../envs/base.yaml"
    params:
        python=sys.executable,
        api_url=config["effector_target_network_network"]["string_api_url"],
        species=config["effector_target_network_network"]["string_species"],
        required_score=config["effector_target_network_network"]["string_required_score"],
        partner_limit=config["effector_target_network_network"]["string_partner_limit"],
        network_type=config["effector_target_network_network"]["string_network_type"],
        caller_identity=config["effector_target_network_network"]["string_caller_identity"]
    log:
        "results/logs/effector_target_network_string_download.log"
    shell:
        "\"{params.python}\" scripts/fetch_effector_target_network_string_network.py "
        "--seed {input.seed} "
        "--mapping-output {output.mapping} "
        "--interactions-output {output.interactions} "
        "--manifest-output {output.manifest} "
        "--api-url {params.api_url} "
        "--species {params.species} "
        "--required-score {params.required_score} "
        "--partner-limit {params.partner_limit} "
        "--network-type {params.network_type} "
        "--caller-identity {params.caller_identity} > {log} 2>&1"


rule build_effector_target_network_string_centrality:
    input:
        seed="results/tables/effector_target_network_effector_target_network_seed.tsv",
        mapping="results/tables/effector_target_network_string_seed_mapping.tsv",
        interactions="results/tables/effector_target_network_string_interactions.tsv"
    output:
        table="results/tables/effector_target_network_string_network_centrality.tsv",
        markdown="results/reports/effector_target_network_string_network_centrality.md"
    conda:
        "../envs/base.yaml"
    params:
        python=sys.executable
    log:
        "results/logs/effector_target_network_string_network_centrality.log"
    shell:
        "\"{params.python}\" scripts/build_effector_target_network_string_centrality.py "
        "--seed {input.seed} "
        "--mapping {input.mapping} "
        "--interactions {input.interactions} "
        "--output {output.table} "
        "--markdown {output.markdown} > {log} 2>&1"


rule build_effector_target_network_string_evidence_channels:
    input:
        centrality="results/tables/effector_target_network_string_network_centrality.tsv",
        interactions="results/tables/effector_target_network_string_interactions.tsv"
    output:
        table="results/tables/effector_target_network_string_evidence_channels.tsv",
        markdown="results/reports/effector_target_network_string_evidence_channels.md"
    conda:
        "../envs/base.yaml"
    params:
        python=sys.executable,
        channel_threshold=config["effector_target_network_network"]["string_channel_score_threshold"]
    log:
        "results/logs/effector_target_network_string_evidence_channels.log"
    shell:
        "\"{params.python}\" scripts/build_effector_target_network_string_evidence_channels.py "
        "--centrality {input.centrality} "
        "--interactions {input.interactions} "
        "--output {output.table} "
        "--markdown {output.markdown} "
        "--channel-threshold {params.channel_threshold} > {log} 2>&1"


rule build_effector_target_network_ros_neighbor_modules:
    input:
        channels="results/tables/effector_target_network_string_evidence_channels.tsv",
        interactions="results/tables/effector_target_network_string_interactions.tsv"
    output:
        table="results/tables/effector_target_network_ros_neighbor_modules.tsv",
        summary="results/tables/effector_target_network_ros_neighbor_module_summary.tsv",
        markdown="results/reports/effector_target_network_ros_neighbor_modules.md"
    conda:
        "../envs/base.yaml"
    params:
        python=sys.executable,
        channel_threshold=config["effector_target_network_network"]["string_channel_score_threshold"]
    log:
        "results/logs/effector_target_network_ros_neighbor_modules.log"
    shell:
        "\"{params.python}\" scripts/build_effector_target_network_ros_neighbor_modules.py "
        "--channels {input.channels} "
        "--interactions {input.interactions} "
        "--output {output.table} "
        "--summary {output.summary} "
        "--markdown {output.markdown} "
        "--channel-threshold {params.channel_threshold} > {log} 2>&1"


rule build_effector_target_network_ros_expression_overlay:
    input:
        modules="results/tables/effector_target_network_ros_neighbor_modules.tsv",
        manifest="results/tables/public_expression_text_download_manifest.tsv",
        text_dir=config["metadata"]["public_expression_text_dir"]
    output:
        table="results/tables/effector_target_network_ros_expression_overlay.tsv",
        markdown="results/reports/effector_target_network_ros_expression_overlay.md"
    conda:
        "../envs/base.yaml"
    params:
        python=sys.executable
    log:
        "results/logs/effector_target_network_ros_expression_overlay.log"
    shell:
        "\"{params.python}\" scripts/build_effector_target_network_ros_expression_overlay.py "
        "--modules {input.modules} "
        "--manifest {input.manifest} "
        "--text-dir {input.text_dir} "
        "--output {output.table} "
        "--markdown {output.markdown} > {log} 2>&1"


rule build_effector_target_network_ros_supplement_overlay:
    input:
        modules="results/tables/effector_target_network_ros_neighbor_modules.tsv",
        manifest="results/tables/public_expression_supplement_download_manifest.tsv"
    output:
        table="results/tables/effector_target_network_ros_supplement_overlay.tsv",
        markdown="results/reports/effector_target_network_ros_supplement_overlay.md"
    conda:
        "../envs/base.yaml"
    benchmark:
        "results/benchmarks/effector_target_network_ros_supplement_overlay.txt"
    params:
        python=sys.executable
    log:
        "results/logs/effector_target_network_ros_supplement_overlay.log"
    shell:
        "\"{params.python}\" scripts/build_effector_target_network_ros_supplement_overlay.py "
        "--modules {input.modules} "
        "--manifest {input.manifest} "
        "--output {output.table} "
        "--markdown {output.markdown} > {log} 2>&1"


rule build_effector_target_network_integrated_effector_target_hypotheses:
    input:
        secretion="results/tables/effector_target_network_effector_secretion_filter.tsv",
        interpro="results/tables/effector_target_network_effector_interproscan_priority.tsv",
        structures="results/tables/effector_target_network_alphafold_server_best_models.tsv",
        mechanisms="results/tables/effector_target_network_mechanism_hypotheses.tsv"
    output:
        effectors="results/tables/effector_target_network_integrated_effector_evidence.tsv",
        pairs="results/tables/effector_target_network_effector_host_target_hypotheses.tsv",
        report="results/reports/effector_target_network_integrated_effector_target_hypotheses.md"
    conda:
        "../envs/base.yaml"
    benchmark:
        "results/benchmarks/effector_target_network_integrated_effector_target_hypotheses.txt"
    params:
        python=sys.executable
    log:
        "results/logs/effector_target_network_integrated_effector_target_hypotheses.log"
    shell:
        "\"{params.python}\" scripts/build_effector_target_network_integrated_effector_target_hypotheses.py "
        "--secretion {input.secretion} "
        "--interpro {input.interpro} "
        "--structures {input.structures} "
        "--mechanisms {input.mechanisms} "
        "--effector-output {output.effectors} "
        "--pair-output {output.pairs} "
        "--report {output.report} > {log} 2>&1"


rule audit_effector_target_network_target_compatibility:
    input:
        pairs="results/tables/effector_target_network_effector_host_target_hypotheses.tsv",
        interpro="results/tables/effector_target_network_effector_interproscan_priority.tsv",
        constraints="config/effector_target_network_protease_family_constraints.tsv",
        compartments="config/effector_target_network_host_target_compartments.tsv"
    output:
        audit="results/tables/effector_target_network_effector_host_compatibility_audit.tsv",
        shortlist="results/tables/effector_target_network_direct_target_shortlist.tsv",
        report="results/reports/effector_target_network_target_compatibility_audit.md"
    conda:
        "../envs/base.yaml"
    benchmark:
        "results/benchmarks/effector_target_network_target_compatibility_audit.txt"
    params:
        python=sys.executable
    log:
        "results/logs/effector_target_network_target_compatibility_audit.log"
    shell:
        "\"{params.python}\" scripts/audit_effector_target_network_target_compatibility.py "
        "--pairs {input.pairs} "
        "--interpro {input.interpro} "
        "--constraints {input.constraints} "
        "--compartments {input.compartments} "
        "--audit-output {output.audit} "
        "--shortlist-output {output.shortlist} "
        "--report {output.report} > {log} 2>&1"


rule build_apoplastic_cell_wall_apoplastic_target_hypotheses:
    input:
        compatibility="results/tables/effector_target_network_effector_host_compatibility_audit.tsv",
        candidates="results/tables/pmas_integrated_mechanism_candidates_pthun_sensitivity.tsv",
        families="config/apoplastic_cell_wall_apoplastic_target_families.tsv"
    output:
        table="results/tables/apoplastic_cell_wall_apoplastic_target_hypotheses.tsv",
        shortlist="results/tables/apoplastic_cell_wall_apoplastic_target_shortlist.tsv",
        report="results/reports/apoplastic_cell_wall_apoplastic_target_hypotheses.md"
    conda:
        "../envs/base.yaml"
    benchmark:
        "results/benchmarks/apoplastic_cell_wall_apoplastic_target_hypotheses.txt"
    params:
        python=sys.executable
    log:
        "results/logs/apoplastic_cell_wall_apoplastic_target_hypotheses.log"
    shell:
        "\"{params.python}\" scripts/build_apoplastic_cell_wall_apoplastic_target_hypotheses.py "
        "--compatibility-audit {input.compatibility} "
        "--candidates {input.candidates} "
        "--families {input.families} "
        "--output {output.table} "
        "--shortlist {output.shortlist} "
        "--report {output.report} > {log} 2>&1"


rule extract_apoplastic_cell_wall_apoplastic_target_proteins:
    input:
        shortlist="results/tables/apoplastic_cell_wall_apoplastic_target_shortlist.tsv",
        genes="results/tables/stable_orthogroup_genes.tsv",
        species=config["metadata"]["species"],
        expression="results/tables/pmas_mechanism_candidate_shortlist.tsv"
    output:
        fasta="data/interim/apoplastic_cell_wall_apoplastic_target_proteins.faa",
        manifest="results/tables/apoplastic_cell_wall_apoplastic_target_protein_manifest.tsv"
    conda:
        "../envs/base.yaml"
    benchmark:
        "results/benchmarks/apoplastic_cell_wall_apoplastic_target_protein_extraction.txt"
    params:
        python=sys.executable
    log:
        "results/logs/apoplastic_cell_wall_apoplastic_target_protein_extraction.log"
    shell:
        "\"{params.python}\" scripts/extract_apoplastic_cell_wall_apoplastic_target_proteins.py "
        "--shortlist {input.shortlist} "
        "--genes {input.genes} "
        "--species {input.species} "
        "--expression {input.expression} "
        "--fasta {output.fasta} "
        "--manifest {output.manifest} > {log} 2>&1"


rule build_apoplastic_cell_wall_sequence_validation:
    input:
        manifest="results/tables/apoplastic_cell_wall_apoplastic_target_protein_manifest.tsv",
        shortlist="results/tables/apoplastic_cell_wall_apoplastic_target_shortlist.tsv",
        annotations=config["metadata"]["protein_function_annotations"],
        deepsig="data/interim/apoplastic_cell_wall_apoplastic_target_deepsig.gff3",
        tmhmm="data/interim/apoplastic_cell_wall_apoplastic_target_tmhmm.tsv"
    output:
        table="results/tables/apoplastic_cell_wall_apoplastic_target_sequence_validation.tsv",
        summary="results/tables/apoplastic_cell_wall_apoplastic_target_sequence_summary.tsv",
        report="results/reports/apoplastic_cell_wall_apoplastic_target_sequence_validation.md"
    conda:
        "../envs/base.yaml"
    benchmark:
        "results/benchmarks/apoplastic_cell_wall_apoplastic_target_sequence_validation.txt"
    params:
        python=sys.executable
    log:
        "results/logs/apoplastic_cell_wall_apoplastic_target_sequence_validation.log"
    shell:
        "\"{params.python}\" scripts/build_apoplastic_cell_wall_sequence_validation.py "
        "--manifest {input.manifest} "
        "--shortlist {input.shortlist} "
        "--annotations {input.annotations} "
        "--deepsig {input.deepsig} "
        "--tmhmm {input.tmhmm} "
        "--output {output.table} "
        "--summary {output.summary} "
        "--report {output.report} > {log} 2>&1"


rule build_apoplastic_cell_wall_sequence_resolved_shortlist:
    input:
        targets="results/tables/apoplastic_cell_wall_apoplastic_target_shortlist.tsv",
        summary="results/tables/apoplastic_cell_wall_apoplastic_target_sequence_summary.tsv"
    output:
        table="results/tables/apoplastic_cell_wall_sequence_resolved_target_shortlist.tsv",
        report="results/reports/apoplastic_cell_wall_sequence_resolved_target_shortlist.md"
    conda:
        "../envs/base.yaml"
    benchmark:
        "results/benchmarks/apoplastic_cell_wall_sequence_resolved_target_shortlist.txt"
    params:
        python=sys.executable
    log:
        "results/logs/apoplastic_cell_wall_sequence_resolved_target_shortlist.log"
    shell:
        "\"{params.python}\" scripts/build_apoplastic_cell_wall_sequence_resolved_shortlist.py "
        "--targets {input.targets} "
        "--summary {input.summary} "
        "--output {output.table} "
        "--report {output.report} > {log} 2>&1"


rule build_apoplastic_cell_wall_tier1_experimental_design:
    input:
        fasta="data/interim/apoplastic_cell_wall_tier1_m8_target_triplet.faa",
        interpro="data/interim/apoplastic_cell_wall_tier1_interproscan.tsv",
        shortlist="results/tables/apoplastic_cell_wall_sequence_resolved_target_shortlist.tsv",
        validation="results/tables/apoplastic_cell_wall_apoplastic_target_sequence_validation.tsv",
        effector_deepsig="data/interim/effector_target_network_deepsig_representatives.gff3",
        nematode_hits="results/tables/nematode_secretome_swissprot_hits.tsv",
        pine_hits="results/tables/swissprot_diamond_hits.tsv",
        alphafold="results/tables/effector_target_network_alphafold_server_best_models.tsv",
        interaction="results/tables/pmas_interaction_feature_annotations.tsv"
    output:
        audit="results/tables/apoplastic_cell_wall_tier1_domain_completeness_audit.tsv",
        assay="results/tables/apoplastic_cell_wall_tier1_cleavage_assay_design.tsv",
        report="results/reports/apoplastic_cell_wall_tier1_domain_completeness_audit.md"
    conda:
        "../envs/base.yaml"
    benchmark:
        "results/benchmarks/apoplastic_cell_wall_tier1_experimental_design.txt"
    params:
        python=sys.executable
    log:
        "results/logs/apoplastic_cell_wall_tier1_experimental_design.log"
    shell:
        "\"{params.python}\" scripts/build_apoplastic_cell_wall_tier1_experimental_design.py "
        "--fasta {input.fasta} "
        "--interpro {input.interpro} "
        "--shortlist {input.shortlist} "
        "--validation {input.validation} "
        "--effector-deepsig {input.effector_deepsig} "
        "--nematode-hits {input.nematode_hits} "
        "--pine-hits {input.pine_hits} "
        "--alphafold {input.alphafold} "
        "--interaction {input.interaction} "
        "--audit {output.audit} "
        "--assay {output.assay} "
        "--report {output.report} > {log} 2>&1"


rule build_apoplastic_cell_wall_tier1_constructs:
    input:
        fasta="data/interim/apoplastic_cell_wall_tier1_m8_target_triplet.faa",
        audit="results/tables/apoplastic_cell_wall_tier1_domain_completeness_audit.tsv"
    output:
        fasta="data/interim/apoplastic_cell_wall_tier1_recombinant_constructs.faa",
        manifest="results/tables/apoplastic_cell_wall_tier1_recombinant_construct_manifest.tsv"
    conda:
        "../envs/base.yaml"
    benchmark:
        "results/benchmarks/apoplastic_cell_wall_tier1_constructs.txt"
    params:
        python=sys.executable
    log:
        "results/logs/apoplastic_cell_wall_tier1_constructs.log"
    shell:
        "\"{params.python}\" scripts/build_apoplastic_cell_wall_tier1_constructs.py "
        "--fasta {input.fasta} "
        "--audit {input.audit} "
        "--output-fasta {output.fasta} "
        "--manifest {output.manifest} > {log} 2>&1"


rule build_candidate_orthogroup_reference:
    input:
        shortlist="results/tables/pmas_mechanism_candidate_shortlist.tsv",
        genes="results/tables/stable_orthogroup_genes.tsv",
        metadata="config/species_metadata.tsv",
        proteomes=stable_pine_proteomes
    output:
        fasta="data/interim/pine_candidate_orthogroup_reference.faa",
        manifest="results/tables/pine_candidate_orthogroup_reference_manifest.tsv"
    conda:
        "../envs/base.yaml"
    benchmark:
        "results/benchmarks/candidate_orthogroup_reference.txt"
    params:
        python=sys.executable
    log:
        "results/logs/candidate_orthogroup_reference.log"
    shell:
        "\"{params.python}\" scripts/build_candidate_orthogroup_reference.py "
        "--shortlist {input.shortlist} "
        "--genes {input.genes} "
        "--metadata {input.metadata} "
        "--fasta-output {output.fasta} "
        "--manifest-output {output.manifest} > {log} 2>&1"


rule build_candidate_orthogroup_diamond_db:
    input:
        "data/interim/pine_candidate_orthogroup_reference.faa"
    output:
        "data/interim/pine_candidate_orthogroup_reference.dmnd"
    conda:
        "../envs/annotation.yaml"
    benchmark:
        "results/benchmarks/candidate_orthogroup_diamond_db.txt"
    log:
        "results/logs/candidate_orthogroup_diamond_db.log"
    shell:
        "diamond makedb --in {input} --db data/interim/pine_candidate_orthogroup_reference > {log} 2>&1"


rule align_pstrobus_to_candidate_orthogroups:
    input:
        query="data/raw/transcriptomes/Pinus_strobus_GIIE01000000_NCBI_TSA_stem_PWN.fasta.gz",
        database="data/interim/pine_candidate_orthogroup_reference.dmnd"
    output:
        "data/interim/pstrobus_vs_pine_candidate_orthogroups.tsv"
    conda:
        "../envs/annotation.yaml"
    benchmark:
        "results/benchmarks/pstrobus_candidate_orthogroup_blastx.txt"
    threads: 16
    log:
        "results/logs/pstrobus_candidate_orthogroup_blastx.log"
    shell:
        "diamond blastx --query {input.query} --db {input.database} --out {output} "
        "--outfmt 6 qseqid sseqid pident length qlen qstart qend slen evalue bitscore "
        "--evalue 1e-5 --max-target-seqs 100 --threads {threads} > {log} 2>&1"


rule build_pstrobus_orthogroup_consensus_mapping:
    input:
        blast="data/interim/pstrobus_vs_pine_candidate_orthogroups.tsv",
        manifest="results/tables/pine_candidate_orthogroup_reference_manifest.tsv"
    output:
        "results/tables/pstrobus_orthogroup_consensus_mapping.tsv"
    conda:
        "../envs/base.yaml"
    benchmark:
        "results/benchmarks/pstrobus_orthogroup_consensus_mapping.txt"
    params:
        python=sys.executable
    log:
        "results/logs/pstrobus_orthogroup_consensus_mapping.log"
    shell:
        "\"{params.python}\" scripts/build_orthogroup_consensus_mapping.py "
        "--blast {input.blast} --manifest {input.manifest} --output {output} > {log} 2>&1"


rule build_cross_species_expression_cross_species_validation:
    input:
        shortlist="results/tables/pmas_mechanism_candidate_shortlist.tsv",
        orthogroups="results/tables/pmas_interaction_orthogroup_summary.tsv",
        pden_mapping="results/tables/cross_species_expression_pden_transcript_mapping.tsv",
        pstrobus_mapping="results/tables/pstrobus_orthogroup_consensus_mapping.tsv",
        pden_tx2gene="data/interim/pden_published_trinity_tx2gene.tsv",
        pstrobus_tx2gene="results/tables/effector_target_network_pstrobus_GIIE_tx2gene.tsv",
        pden_primary="data/external/galaxy/effector_target_network_pden_full_deseq2/effector_target_network_bxyl_vs_bthai_result.tabular",
        pden_bxyl_water="data/external/galaxy/effector_target_network_pden_full_deseq2/effector_target_network_bxyl_vs_water_result.tabular",
        pden_bthai_water="data/external/galaxy/effector_target_network_pden_full_deseq2/effector_target_network_bthai_vs_water_result.tabular",
        pstrobus_two_week="results/tables/effector_target_network_pstrobus_GIIE_DESeq2_2w_vs_0w.tsv",
        pstrobus_four_week="results/tables/effector_target_network_pstrobus_GIIE_DESeq2_4w_vs_0w.tsv",
        stable_genes="results/tables/stable_orthogroup_genes.tsv",
        hybrid_two_week="results/tables/apoplastic_cell_wall_prig_xtae_DESeq2_2w_vs_0w.tsv",
        hybrid_four_week="results/tables/apoplastic_cell_wall_prig_xtae_DESeq2_4w_vs_0w.tsv",
        pthun_two_week="results/tables/apoplastic_cell_wall_pthun_DESeq2_2w_vs_0w.tsv",
        pthun_four_week="results/tables/apoplastic_cell_wall_pthun_DESeq2_4w_vs_0w.tsv"
    output:
        "results/tables/pmas_mechanism_cross_species_validation.tsv"
    conda:
        "../envs/base.yaml"
    benchmark:
        "results/benchmarks/cross_species_expression_cross_species_validation.txt"
    params:
        python=sys.executable
    log:
        "results/logs/cross_species_expression_cross_species_validation.log"
    shell:
        "\"{params.python}\" scripts/build_cross_species_expression_cross_species_validation.py "
        "--shortlist {input.shortlist} --orthogroups {input.orthogroups} "
        "--pden-mapping {input.pden_mapping} --pstrobus-mapping {input.pstrobus_mapping} "
        "--pden-tx2gene {input.pden_tx2gene} --pstrobus-tx2gene {input.pstrobus_tx2gene} "
        "--pden-primary {input.pden_primary} --pden-bxyl-water {input.pden_bxyl_water} "
        "--pden-bthai-water {input.pden_bthai_water} "
        "--pstrobus-two-week {input.pstrobus_two_week} --pstrobus-four-week {input.pstrobus_four_week} "
        "--stable-genes {input.stable_genes} "
        "--prig-xtae-two-week {input.hybrid_two_week} --prig-xtae-four-week {input.hybrid_four_week} "
        "--pthun-two-week {input.pthun_two_week} --pthun-four-week {input.pthun_four_week} "
        "--output {output} > {log} 2>&1"


rule build_host_mechanism_evidence_integration:
    input:
        shortlist="results/tables/pmas_mechanism_candidate_shortlist.tsv",
        axes="results/tables/mechanism_axis_summary.tsv",
        validation="results/tables/pmas_mechanism_cross_species_validation.tsv",
        centrality="results/tables/pmas_mechanism_string_centrality.tsv",
        interactions="results/tables/pmas_mechanism_string_interactions.tsv",
        effector_target_network="results/tables/effector_target_network_effector_host_target_hypotheses.tsv",
        apoplastic_cell_wall="results/tables/apoplastic_cell_wall_apoplastic_target_hypotheses.tsv",
        effectors="results/tables/effector_target_network_integrated_effector_evidence.tsv"
    output:
        priority="results/tables/mechanism_host_candidate_priority.tsv",
        audit="results/tables/mechanism_host_evidence_audit.tsv",
        convergence="results/tables/mechanism_route_convergence.tsv",
        report="results/reports/mechanism_host_candidate_priority.md"
    conda:
        "../envs/base.yaml"
    benchmark:
        "results/benchmarks/host_mechanism_evidence_integration.txt"
    params:
        python=sys.executable
    log:
        "results/logs/host_mechanism_evidence_integration.log"
    shell:
        "\"{params.python}\" scripts/build_host_mechanism_evidence_integration.py "
        "--shortlist {input.shortlist} "
        "--axes {input.axes} "
        "--validation {input.validation} "
        "--centrality {input.centrality} "
        "--interactions {input.interactions} "
        "--effector_target_network-hypotheses {input.effector_target_network} "
        "--apoplastic_cell_wall-hypotheses {input.apoplastic_cell_wall} "
        "--effectors {input.effectors} "
        "--priority-output {output.priority} "
        "--audit-output {output.audit} "
        "--convergence-output {output.convergence} "
        "--report-output {output.report} > {log} 2>&1"


rule build_og0005853_pr4_dossier:
    input:
        stable_genes="results/tables/stable_orthogroup_genes.tsv",
        metadata="config/species_metadata.tsv",
        fasta="results/orthofinder/pilot_run/Results_Jul03/Orthogroup_Sequences/OG0005853.fa",
        tree="results/orthofinder/pilot_run/Results_Jul03/Gene_Trees/OG0005853_tree.txt",
        swissprot="results/tables/swissprot_diamond_hits.tsv",
        defense_hits="results/tables/defense_module_gene_hits.tsv",
        pden_mapping="results/tables/cross_species_expression_pden_transcript_mapping.tsv",
        pstrobus_mapping="results/tables/pstrobus_orthogroup_consensus_mapping.tsv",
        pden_tx2gene="data/interim/pden_published_trinity_tx2gene.tsv",
        pstrobus_tx2gene="results/tables/effector_target_network_pstrobus_GIIE_tx2gene.tsv",
        pden_primary="data/external/galaxy/effector_target_network_pden_full_deseq2/effector_target_network_bxyl_vs_bthai_result.tabular",
        pden_bxyl_water="data/external/galaxy/effector_target_network_pden_full_deseq2/effector_target_network_bxyl_vs_water_result.tabular",
        pden_bthai_water="data/external/galaxy/effector_target_network_pden_full_deseq2/effector_target_network_bthai_vs_water_result.tabular",
        pstrobus_two_week="results/tables/effector_target_network_pstrobus_GIIE_DESeq2_2w_vs_0w.tsv",
        pstrobus_four_week="results/tables/effector_target_network_pstrobus_GIIE_DESeq2_4w_vs_0w.tsv"
    output:
        members="results/tables/og0005853_member_annotation.tsv",
        expression="results/tables/og0005853_expression_evidence.tsv",
        tree_audit="results/tables/og0005853_tree_audit.tsv",
        report="results/reports/og0005853_pr4_mechanism_dossier.md"
    conda:
        "../envs/base.yaml"
    benchmark:
        "results/benchmarks/og0005853_pr4_dossier.txt"
    params:
        python=sys.executable
    log:
        "results/logs/og0005853_pr4_dossier.log"
    shell:
        "\"{params.python}\" scripts/build_og0005853_pr4_dossier.py "
        "--stable-genes {input.stable_genes} --metadata {input.metadata} "
        "--fasta {input.fasta} --tree {input.tree} "
        "--swissprot {input.swissprot} --defense-hits {input.defense_hits} "
        "--pden-mapping {input.pden_mapping} --pstrobus-mapping {input.pstrobus_mapping} "
        "--pden-tx2gene {input.pden_tx2gene} --pstrobus-tx2gene {input.pstrobus_tx2gene} "
        "--pden-primary {input.pden_primary} --pden-bxyl-water {input.pden_bxyl_water} "
        "--pden-bthai-water {input.pden_bthai_water} "
        "--pstrobus-two-week {input.pstrobus_two_week} "
        "--pstrobus-four-week {input.pstrobus_four_week} "
        "--member-output {output.members} --expression-output {output.expression} "
        "--tree-audit-output {output.tree_audit} --report-output {output.report} "
        "> {log} 2>&1"


rule build_manuscript_mechanism_summary:
    input:
        apoplastic_cell_wall_hypotheses="results/tables/apoplastic_cell_wall_apoplastic_target_hypotheses.tsv",
        host_priorities="results/tables/mechanism_host_candidate_priority.tsv",
        sequence_shortlist="results/tables/apoplastic_cell_wall_sequence_resolved_target_shortlist.tsv",
        og0005853_audit="results/tables/og0005853_tree_audit.tsv",
        og0005853_expression="results/tables/og0005853_expression_evidence.tsv"
    output:
        apoplastic_cell_wall_audit="results/tables/apoplastic_cell_wall_multinode_independence_audit.tsv",
        apoplastic_cell_wall_summary="results/tables/apoplastic_cell_wall_multinode_independence_summary.tsv",
        apoplastic_cell_wall_report="results/reports/apoplastic_cell_wall_multinode_independence_audit.md",
        candidates="results/tables/manuscript_mechanism_candidate_summary.tsv",
        candidate_report="results/reports/manuscript_mechanism_candidate_summary.md"
    conda:
        "../envs/base.yaml"
    benchmark:
        "results/benchmarks/manuscript_mechanism_summary.txt"
    params:
        python=sys.executable
    log:
        "results/logs/manuscript_mechanism_summary.log"
    shell:
        "\"{params.python}\" scripts/build_manuscript_mechanism_summary.py "
        "--apoplastic_cell_wall-hypotheses {input.apoplastic_cell_wall_hypotheses} "
        "--host-priorities {input.host_priorities} "
        "--sequence-shortlist {input.sequence_shortlist} "
        "--og0005853-audit {input.og0005853_audit} "
        "--og0005853-expression {input.og0005853_expression} "
        "--apoplastic_cell_wall-audit-output {output.apoplastic_cell_wall_audit} "
        "--apoplastic_cell_wall-summary-output {output.apoplastic_cell_wall_summary} "
        "--apoplastic_cell_wall-report-output {output.apoplastic_cell_wall_report} "
        "--candidate-output {output.candidates} "
        "--candidate-report-output {output.candidate_report} > {log} 2>&1"
