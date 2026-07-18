FIGURES = [
    "results/figures/Fig1_project_design.pdf",
    "results/figures/Fig2_secretome_audit.pdf",
    "results/figures/Fig3_integrated_candidate_evidence.pdf",
    "results/figures/Fig5_apoplastic_validation_path.pdf",
    "results/figures/Fig2_defense_module_heatmap.pdf",
    "results/figures/Fig3_effector_classes.pdf",
    "results/figures/Fig5_og0005853_family.pdf",
    "results/figures/Fig6_integrated_cross_species_evidence_landscape.pdf",
    "results/figures/Fig7_effector_informed_functional_prior_response_map.pdf",
    "results/figures/Fig7_multilayer_candidate_circos.pdf",
    "results/figures/FigS_cafe5_candidate_turnover_sensitivity.pdf",
]


rule build_candidate_proteome_sensitivity:
    input:
        script="scripts/build_candidate_proteome_sensitivity.py",
        stable="results/tables/stable_orthogroups.tsv",
        candidates="results/tables/manuscript_mechanism_candidate_summary.tsv"
    output:
        "results/tables/manuscript_candidate_proteome_sensitivity.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/candidate_proteome_sensitivity.log"
    benchmark:
        "results/benchmarks/candidate_proteome_sensitivity.tsv"
    shell:
        "python {input.script} --stable-orthogroups {input.stable} "
        "--candidates {input.candidates} --output {output} > {log} 2>&1"


rule build_secretome_audit_figure_data:
    input:
        script="scripts/build_secretome_audit_figure_data.py",
        integrated="results/tables/nematode_secretome_integrated_evidence.tsv"
    output:
        stages="results/tables/manuscript_secretome_audit_stages.tsv",
        combinations="results/tables/manuscript_secretome_audit_upset.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/secretome_audit_figure_data.log"
    benchmark:
        "results/benchmarks/secretome_audit_figure_data.tsv"
    shell:
        "python {input.script} --integrated {input.integrated} "
        "--stages-output {output.stages} --combinations-output {output.combinations} "
        "> {log} 2>&1"


rule plot_secretome_audit:
    input:
        script="scripts/plot_secretome_audit.R",
        stages="results/tables/manuscript_secretome_audit_stages.tsv",
        combinations="results/tables/manuscript_secretome_audit_upset.tsv"
    output:
        pdf="results/figures/Fig2_secretome_audit.pdf",
        png="results/figures/Fig2_secretome_audit.png"
    conda:
        "../envs/r.yaml"
    log:
        "results/logs/secretome_audit_figure.log"
    benchmark:
        "results/benchmarks/secretome_audit_figure.tsv"
    shell:
        "Rscript {input.script} {input.stages} {input.combinations} "
        "{output.pdf} {output.png} > {log} 2>&1"


rule plot_integrated_candidate_evidence:
    input:
        script="scripts/plot_integrated_candidate_evidence.R",
        bubble="results/tables/manuscript_attack_defense_bubble.tsv",
        member_audit="results/tables/manuscript_expression_member_audit.tsv",
        candidates="results/tables/manuscript_mechanism_candidate_summary.tsv"
    output:
        pdf="results/figures/Fig3_integrated_candidate_evidence.pdf",
        png="results/figures/Fig3_integrated_candidate_evidence.png"
    conda:
        "../envs/r.yaml"
    log:
        "results/logs/integrated_candidate_evidence.log"
    benchmark:
        "results/benchmarks/integrated_candidate_evidence.tsv"
    shell:
        "Rscript {input.script} {input.bubble} {input.member_audit} {input.candidates} "
        "{output.pdf} {output.png} > {log} 2>&1"


rule build_attack_defense_figure_data:
    input:
        script="scripts/build_attack_defense_figure_data.py",
        candidates="results/tables/manuscript_multilayer_circos_candidates.tsv",
        member_audit="results/tables/manuscript_expression_member_audit.tsv",
        mechanism_edges="results/tables/manuscript_mechanism_chord_edges.tsv",
        effector_context="results/tables/manuscript_circos_effector_class_context.tsv",
        sectors="results/tables/manuscript_mechanism_chord_sectors.tsv"
    output:
        bubble="results/tables/manuscript_attack_defense_bubble.tsv",
        nodes="results/tables/manuscript_attack_defense_nodes.tsv",
        edges="results/tables/manuscript_attack_defense_edges.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/attack_defense_figure_data.log"
    benchmark:
        "results/benchmarks/attack_defense_figure_data.tsv"
    shell:
        "python {input.script} "
        "--candidates {input.candidates} --member-audit {input.member_audit} "
        "--mechanism-edges {input.mechanism_edges} "
        "--effector-context {input.effector_context} --sectors {input.sectors} "
        "--bubble-output {output.bubble} --nodes-output {output.nodes} "
        "--edges-output {output.edges} > {log} 2>&1"


rule plot_pwn_pine_attack_defense:
    input:
        script="scripts/plot_pwn_pine_attack_defense.R",
        bubble="results/tables/manuscript_attack_defense_bubble.tsv",
        nodes="results/tables/manuscript_attack_defense_nodes.tsv",
        edges="results/tables/manuscript_attack_defense_edges.tsv"
    output:
        pdf="results/figures/Fig7_effector_informed_functional_prior_response_map.pdf",
        png="results/figures/Fig7_effector_informed_functional_prior_response_map.png"
    conda:
        "../envs/r.yaml"
    log:
        "results/logs/functional_prior_response_map.log"
    benchmark:
        "results/benchmarks/functional_prior_response_map.tsv"
    shell:
        "Rscript {input.script} {input.bubble} {input.nodes} {input.edges} "
        "{output.pdf} {output.png} > {log} 2>&1"


rule build_evidence_sector_circos_data:
    input:
        script="scripts/build_evidence_sector_circos_data.py",
        candidates="results/tables/manuscript_multilayer_circos_candidates.tsv",
        effector_context="results/tables/manuscript_circos_effector_class_context.tsv",
        attack_edges="results/tables/manuscript_attack_defense_edges.tsv",
        cross_species="results/tables/pmas_mechanism_cross_species_validation.tsv",
        member_audit="results/tables/manuscript_expression_member_audit.tsv"
    output:
        sectors="results/tables/manuscript_evidence_circos_sectors_v1.tsv",
        tracks="results/tables/manuscript_evidence_circos_tracks_long_v1.tsv",
        links="results/tables/manuscript_evidence_circos_links_v1.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/evidence_sector_circos_data.log"
    benchmark:
        "results/benchmarks/evidence_sector_circos_data.tsv"
    shell:
        "python {input.script} --candidates {input.candidates} "
        "--effector-context {input.effector_context} --attack-edges {input.attack_edges} "
        "--cross-species-expression {input.cross_species} "
        "--member-audit {input.member_audit} "
        "--sectors-output {output.sectors} --tracks-output {output.tracks} "
        "--links-output {output.links} > {log} 2>&1"


rule build_evidence_circos_gene_sample_mosaic:
    input:
        script="scripts/build_evidence_circos_gene_sample_mosaic.py",
        candidates="results/tables/manuscript_multilayer_circos_candidates.tsv",
        member_audit="results/tables/manuscript_expression_member_audit.tsv",
        stable_genes="results/tables/stable_orthogroup_genes.tsv",
        pden_mapping="results/tables/cross_species_expression_pden_transcript_mapping.tsv",
        pstrobus_mapping="results/tables/pstrobus_orthogroup_consensus_mapping.tsv",
        pmas_tx2gene="results/tables/effector_target_network_pmas_GigaDB102688_tx2gene.tsv",
        pden_tx2gene="data/interim/pden_published_trinity_tx2gene.tsv",
        pstrobus_tx2gene="results/tables/effector_target_network_pstrobus_GIIE_tx2gene.tsv",
        pmas_samples="config/effector_target_network_pmas_1dpi_primary_sample_sheet.tsv",
        pden_samples="config/effector_target_network_pden_deseq2_sample_sheet.tsv",
        pstrobus_samples="config/effector_target_network_pstrobus_deseq2_sample_sheet.tsv",
        pthun_samples="config/apoplastic_cell_wall_pthun_deseq2_sample_sheet.tsv",
        hybrid_samples="config/apoplastic_cell_wall_prig_xtae_deseq2_sample_sheet.tsv",
        pstrobus_manifest="results/tables/effector_target_network_pstrobus_quant_manifest.tsv",
        pden_counts="data/external/galaxy/effector_target_network_pden_full_salmon/effector_target_network_pden_tximport_gene_level_values.tabular",
        pden_de_1="data/external/galaxy/effector_target_network_pden_full_deseq2/effector_target_network_bxyl_vs_bthai_result.tabular",
        pden_de_2="data/external/galaxy/effector_target_network_pden_full_deseq2/effector_target_network_bxyl_vs_water_result.tabular",
        pden_de_3="data/external/galaxy/effector_target_network_pden_full_deseq2/effector_target_network_bthai_vs_water_result.tabular",
        pstrobus_de_2w="results/tables/effector_target_network_pstrobus_GIIE_DESeq2_2w_vs_0w.tsv",
        pstrobus_de_4w="results/tables/effector_target_network_pstrobus_GIIE_DESeq2_4w_vs_0w.tsv",
        pthun_de_2w="results/tables/apoplastic_cell_wall_pthun_DESeq2_2w_vs_0w.tsv",
        pthun_de_4w="results/tables/apoplastic_cell_wall_pthun_DESeq2_4w_vs_0w.tsv",
        hybrid_de_2w="results/tables/apoplastic_cell_wall_prig_xtae_DESeq2_2w_vs_0w.tsv",
        hybrid_de_4w="results/tables/apoplastic_cell_wall_prig_xtae_DESeq2_4w_vs_0w.tsv"
    output:
        "results/tables/manuscript_evidence_circos_gene_sample_mosaic_v1.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/evidence_circos_gene_sample_mosaic.log"
    benchmark:
        "results/benchmarks/evidence_circos_gene_sample_mosaic.tsv"
    shell:
        "python {input.script} --project-root . --output {output} > {log} 2>&1"


rule plot_evidence_sector_circos_preview:
    input:
        script="scripts/plot_evidence_sector_circos.R",
        sectors=ancient("results/tables/manuscript_evidence_circos_sectors_v1.tsv"),
        tracks=ancient("results/tables/manuscript_evidence_circos_tracks_long_v1.tsv"),
        links=ancient("results/tables/manuscript_evidence_circos_links_v1.tsv"),
        mosaic=ancient("results/tables/manuscript_evidence_circos_gene_sample_mosaic_v1.tsv")
    output:
        pdf="results/figures/FigS_full_cross_species_evidence_landscape.pdf",
        png="results/figures/FigS_full_cross_species_evidence_landscape.png"
    conda:
        "../envs/r.yaml"
    log:
        "results/logs/full_cross_species_evidence_landscape.log"
    benchmark:
        "results/benchmarks/full_cross_species_evidence_landscape.tsv"
    shell:
        "Rscript {input.script} {input.sectors} {input.tracks} {input.links} {input.mosaic} "
        "{output.pdf} {output.png} > {log} 2>&1"


rule build_main_circos_tables:
    input:
        script="scripts/build_main_circos_tables.py",
        sectors=ancient("results/tables/manuscript_evidence_circos_sectors_v1.tsv"),
        tracks=ancient("results/tables/manuscript_evidence_circos_tracks_long_v1.tsv"),
        links=ancient("results/tables/manuscript_evidence_circos_links_v1.tsv"),
        mosaic=ancient("results/tables/manuscript_evidence_circos_gene_sample_mosaic_v1.tsv")
    output:
        sectors="results/tables/manuscript_main_circos_sectors_v1.tsv",
        tracks="results/tables/manuscript_main_circos_tracks_long_v1.tsv",
        links="results/tables/manuscript_main_circos_links_v1.tsv",
        mosaic="results/tables/manuscript_main_circos_gene_sample_mosaic_v1.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/main_circos_tables.log"
    benchmark:
        "results/benchmarks/main_circos_tables.tsv"
    shell:
        "python {input.script} --sectors {input.sectors} --tracks {input.tracks} "
        "--links {input.links} --mosaic {input.mosaic} "
        "--sectors-output {output.sectors} --tracks-output {output.tracks} "
        "--links-output {output.links} --mosaic-output {output.mosaic} > {log} 2>&1"


rule plot_main_integrated_circos:
    input:
        script="scripts/plot_evidence_sector_circos.R",
        sectors="results/tables/manuscript_main_circos_sectors_v1.tsv",
        tracks="results/tables/manuscript_main_circos_tracks_long_v1.tsv",
        links="results/tables/manuscript_main_circos_links_v1.tsv",
        mosaic="results/tables/manuscript_main_circos_gene_sample_mosaic_v1.tsv"
    output:
        pdf="results/figures/Fig6_integrated_cross_species_evidence_landscape.pdf",
        png="results/figures/Fig6_integrated_cross_species_evidence_landscape.png"
    conda:
        "../envs/r.yaml"
    log:
        "results/logs/main_integrated_circos.log"
    benchmark:
        "results/benchmarks/main_integrated_circos.tsv"
    shell:
        "Rscript {input.script} {input.sectors} {input.tracks} {input.links} {input.mosaic} "
        "{output.pdf} {output.png} main > {log} 2>&1"


rule build_manuscript_revision_audits:
    input:
        script="scripts/build_manuscript_revision_audits.py",
        pmas_samples="config/effector_target_network_pmas_1dpi_primary_sample_sheet.tsv",
        pden_samples="config/effector_target_network_pden_deseq2_sample_sheet.tsv",
        pstrobus_samples="config/effector_target_network_pstrobus_deseq2_sample_sheet.tsv",
        stable_genes="results/tables/stable_orthogroup_genes.tsv",
        defense_matrix="results/tables/orthogroup_defense_matrix.tsv",
        module_links="results/tables/effector_host_module_links.tsv",
        priorities="results/tables/mechanism_host_candidate_priority.tsv",
        candidates="results/tables/manuscript_mechanism_candidate_summary.tsv",
        circos_candidates="results/tables/manuscript_multilayer_circos_candidates.tsv",
        pmas_resistant="results/tables/effector_target_network_pmas_GigaDB102688_DESeq2_resistant_pwn_vs_water.tsv",
        pmas_susceptible="results/tables/effector_target_network_pmas_GigaDB102688_DESeq2_susceptible_pwn_vs_water.tsv",
        pmas_interaction="results/tables/effector_target_network_pmas_GigaDB102688_DESeq2_genotype_by_inoculum.tsv",
        cross_species="results/tables/pmas_mechanism_cross_species_validation.tsv",
        apoplastic_audit="results/tables/apoplastic_cell_wall_multinode_independence_audit.tsv"
    output:
        datasets="results/tables/manuscript_expression_dataset_summary.tsv",
        attrition="results/tables/manuscript_candidate_attrition.tsv",
        rules="results/tables/manuscript_priority_rule_matrix.tsv",
        members="results/tables/manuscript_expression_member_audit.tsv",
        biochemical="results/tables/manuscript_biochemical_test_candidates.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/manuscript_revision_audits.log"
    benchmark:
        "results/benchmarks/manuscript_revision_audits.tsv"
    shell:
        "python {input.script} --root . --output-dir results/tables > {log} 2>&1"


rule build_mechanism_chord_inputs:
    input:
        candidates="results/tables/manuscript_mechanism_candidate_summary.tsv",
        links="results/tables/effector_host_module_links.tsv",
        effectors="results/tables/nematode_candidate_effectors.tsv"
    output:
        edges="results/tables/manuscript_mechanism_chord_edges.tsv",
        sectors="results/tables/manuscript_mechanism_chord_sectors.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/mechanism_chord_inputs.log"
    benchmark:
        "results/benchmarks/mechanism_chord_inputs.tsv"
    shell:
        "python scripts/build_mechanism_chord_inputs.py "
        "--candidates {input.candidates} --links {input.links} --effectors {input.effectors} "
        "--edges-output {output.edges} --sectors-output {output.sectors} > {log} 2>&1"


rule plot_mechanism_chord:
    input:
        edges="results/tables/manuscript_mechanism_chord_edges.tsv",
        sectors="results/tables/manuscript_mechanism_chord_sectors.tsv"
    output:
        pdf="results/figures/Fig7_mechanism_chord.pdf",
        png="results/figures/Fig7_mechanism_chord.png"
    conda:
        "../envs/r.yaml"
    log:
        "results/logs/mechanism_chord.log"
    benchmark:
        "results/benchmarks/mechanism_chord.tsv"
    shell:
        "Rscript scripts/plot_mechanism_chord.R {input.edges} {input.sectors} "
        "{output.pdf} {output.png} > {log} 2>&1"


rule build_multilayer_candidate_circos:
    input:
        script="scripts/build_multilayer_candidate_circos.py",
        candidates="results/tables/manuscript_mechanism_candidate_summary.tsv",
        stable_genes="results/tables/stable_orthogroup_genes.tsv",
        regional="results/tables/east_asia_vs_north_america_candidate_contrast.tsv",
        resistant="results/tables/effector_target_network_pmas_GigaDB102688_DESeq2_resistant_pwn_vs_water.tsv",
        susceptible="results/tables/effector_target_network_pmas_GigaDB102688_DESeq2_susceptible_pwn_vs_water.tsv",
        interaction="results/tables/effector_target_network_pmas_GigaDB102688_DESeq2_genotype_by_inoculum.tsv",
        interaction_summary="results/tables/pmas_interaction_orthogroup_summary.tsv",
        cross_species="results/tables/pmas_mechanism_cross_species_validation.tsv",
        effectors="results/tables/nematode_candidate_effectors.tsv",
        secretome_audit="results/tables/nematode_full_secretome_evidence.tsv",
        edges="results/tables/manuscript_mechanism_chord_edges.tsv"
    output:
        candidates="results/tables/manuscript_multilayer_circos_candidates.tsv",
        effector_context="results/tables/manuscript_circos_effector_class_context.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/multilayer_candidate_circos_inputs.log"
    benchmark:
        "results/benchmarks/multilayer_candidate_circos_inputs.tsv"
    shell:
        "python {input.script} "
        "--candidates {input.candidates} --stable-genes {input.stable_genes} "
        "--regional-contrasts {input.regional} --resistant {input.resistant} "
        "--susceptible {input.susceptible} --interaction {input.interaction} "
        "--interaction-summary {input.interaction_summary} "
        "--cross-species-validation {input.cross_species} "
        "--effectors {input.effectors} --secretome-audit {input.secretome_audit} "
        "--edges {input.edges} --output {output.candidates} "
        "--effector-context-output {output.effector_context} > {log} 2>&1"


rule plot_multilayer_candidate_circos:
    input:
        script="scripts/plot_multilayer_candidate_circos.R",
        freeze_script="scripts/freeze_circos_figure_inputs.py",
        freeze_manifest="results/tables/manuscript_circos_figure_data_freeze_v1.tsv",
        candidates="results/tables/manuscript_multilayer_circos_candidates.tsv",
        effector_context="results/tables/manuscript_circos_effector_class_context.tsv",
        edges="results/tables/manuscript_mechanism_chord_edges.tsv",
        sectors="results/tables/manuscript_mechanism_chord_sectors.tsv"
    output:
        pdf="results/figures/Fig7_multilayer_candidate_circos.pdf",
        png="results/figures/Fig7_multilayer_candidate_circos.png"
    conda:
        "../envs/r.yaml"
    log:
        "results/logs/multilayer_candidate_circos.log"
    benchmark:
        "results/benchmarks/multilayer_candidate_circos.tsv"
    shell:
        "python {input.freeze_script} --mode verify --freeze-version circos_v1 "
        "--candidate-tracks {input.candidates} --effector-context {input.effector_context} "
        "--edges {input.edges} --sectors {input.sectors} --manifest {input.freeze_manifest} "
        "&& Rscript {input.script} "
        "{input.candidates} {input.effector_context} {input.edges} {input.sectors} "
        "{output.pdf} {output.png} > {log} 2>&1"


rule plot_effector_module_candidate_network:
    input:
        script="scripts/plot_effector_module_candidate_network.R",
        freeze_script="scripts/freeze_circos_figure_inputs.py",
        freeze_manifest="results/tables/manuscript_circos_figure_data_freeze_v1.tsv",
        candidates="results/tables/manuscript_multilayer_circos_candidates.tsv",
        effector_context="results/tables/manuscript_circos_effector_class_context.tsv",
        edges="results/tables/manuscript_mechanism_chord_edges.tsv",
        sectors="results/tables/manuscript_mechanism_chord_sectors.tsv"
    output:
        pdf="results/figures/Fig7_effector_module_candidate_network.pdf",
        png="results/figures/Fig7_effector_module_candidate_network.png"
    conda:
        "../envs/r.yaml"
    log:
        "results/logs/effector_module_candidate_network.log"
    benchmark:
        "results/benchmarks/effector_module_candidate_network.tsv"
    shell:
        "python {input.freeze_script} --mode verify --freeze-version circos_v1 "
        "--candidate-tracks {input.candidates} --effector-context {input.effector_context} "
        "--edges {input.edges} --sectors {input.sectors} --manifest {input.freeze_manifest} "
        "&& Rscript {input.script} "
        "{input.candidates} {input.effector_context} {input.edges} {input.sectors} "
        "{output.pdf} {output.png} > {log} 2>&1"


rule build_cafe5_candidate_turnover_figure:
    input:
        script="scripts/build_cafe5_candidate_turnover_figure.py",
        candidates="results/tables/manuscript_mechanism_candidate_summary.tsv",
        families="results/tables/cafe5_full_family_results.tsv",
        branches="results/tables/cafe5_full_branch_results.tsv",
        model="results/tables/cafe5_full_model_summary.tsv",
        diagnostics="results/tables/cafe5_diagnostic_model_comparison.tsv",
        convergence="results/tables/cafe5_diagnostic_convergence.tsv"
    output:
        "results/tables/cafe5_candidate_turnover_figure.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/cafe5_candidate_turnover_figure.log"
    benchmark:
        "results/benchmarks/cafe5_candidate_turnover_figure.tsv"
    shell:
        "python {input.script} "
        "--candidates {input.candidates} --family-results {input.families} "
        "--branch-results {input.branches} --model-summary {input.model} "
        "--diagnostics {input.diagnostics} --convergence {input.convergence} "
        "--output {output} > {log} 2>&1"


rule plot_cafe5_candidate_turnover_sensitivity:
    input:
        script="scripts/plot_cafe5_candidate_turnover_sensitivity.R",
        table="results/tables/cafe5_candidate_turnover_figure.tsv"
    output:
        pdf="results/figures/FigS_cafe5_candidate_turnover_sensitivity.pdf",
        png="results/figures/FigS_cafe5_candidate_turnover_sensitivity.png"
    conda:
        "../envs/r.yaml"
    log:
        "results/logs/cafe5_candidate_turnover_sensitivity.log"
    benchmark:
        "results/benchmarks/cafe5_candidate_turnover_sensitivity.tsv"
    shell:
        "Rscript {input.script} {input.table} "
        "{output.pdf} {output.png} > {log} 2>&1"

rule plot_manuscript_core_figures:
    input:
        species="results/tables/species_summary_longest_isoforms.tsv",
        defense="results/tables/defense_module_species_summary.tsv",
        orthology="results/tables/of_pilot_vs_expansion_summary.tsv",
        candidates="results/tables/manuscript_mechanism_candidate_summary.tsv",
        effectors="results/tables/nematode_candidate_effectors.tsv",
        links="results/tables/effector_host_module_links.tsv",
        og_members="results/tables/og0005853_member_annotation.tsv",
        og_expression="results/tables/og0005853_expression_evidence.tsv",
        og_tree="results/orthofinder/pilot_run/Results_Jul03/Gene_Trees/OG0005853_tree.txt"
    output:
        fig1_pdf="results/figures/Fig1_project_design.pdf",
        fig1_png="results/figures/Fig1_project_design.png",
        fig2_pdf="results/figures/Fig2_defense_module_heatmap.pdf",
        fig2_png="results/figures/Fig2_defense_module_heatmap.png",
        fig3_pdf="results/figures/Fig3_effector_classes.pdf",
        fig3_png="results/figures/Fig3_effector_classes.png",
        fig5_pdf="results/figures/Fig5_og0005853_family.pdf",
        fig5_png="results/figures/Fig5_og0005853_family.png"
    conda:
        "../envs/r.yaml"
    log:
        "results/logs/manuscript_core_figures.log"
    benchmark:
        "results/benchmarks/manuscript_core_figures.tsv"
    shell:
        "Rscript scripts/plot_manuscript_core_figures.R "
        "{input.species} {input.defense} {input.orthology} {input.candidates} "
        "{input.effectors} {input.links} {input.og_members} {input.og_expression} {input.og_tree} "
        "{output.fig1_pdf} {output.fig1_png} {output.fig2_pdf} {output.fig2_png} "
        "{output.fig3_pdf} {output.fig3_png} {output.fig5_pdf} {output.fig5_png} "
        "> {log} 2>&1"


rule plot_manuscript_evidence_figures:
    input:
        candidates="results/tables/manuscript_mechanism_candidate_summary.tsv",
        audit="results/tables/apoplastic_cell_wall_multinode_independence_audit.tsv"
    output:
        fig4_pdf="results/figures/FigS_candidate_evidence_matrix.pdf",
        fig4_png="results/figures/FigS_candidate_evidence_matrix.png",
        fig6_pdf="results/figures/Fig5_apoplastic_validation_path.pdf",
        fig6_png="results/figures/Fig5_apoplastic_validation_path.png"
    conda:
        "../envs/r.yaml"
    log:
        "results/logs/manuscript_evidence_figures.log"
    benchmark:
        "results/benchmarks/manuscript_evidence_figures.tsv"
    shell:
        "Rscript scripts/plot_manuscript_evidence_figures.R "
        "{input.candidates} {input.audit} "
        "{output.fig4_pdf} {output.fig4_png} {output.fig6_pdf} {output.fig6_png} "
        "> {log} 2>&1"


rule summarize_candidate_coevolution_for_fig5:
    input:
        "results/tables/candidate_coevolution_genes.tsv"
    output:
        summary="results/tables/candidate_coevolution_summary.tsv",
        top="results/tables/candidate_coevolution_top_candidates.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/candidate_coevolution_summary.log"
    shell:
        "python scripts/summarize_candidate_coevolution.py "
        "--candidates {input} "
        "--summary-output {output.summary} "
        "--top-output {output.top} "
        "--max-per-group 5 > {log} 2>&1"


rule plot_fig5_candidate_coevolution:
    input:
        summary="results/tables/candidate_coevolution_summary.tsv",
        top="results/tables/candidate_coevolution_top_candidates.tsv"
    output:
        "results/figures/Fig5_candidate_coevolution.pdf"
    conda:
        "../envs/r.yaml"
    log:
        "results/logs/Fig5_candidate_coevolution.log"
    shell:
        "Rscript scripts/plot_fig5_candidate_coevolution.R "
        "{input.summary} {input.top} {output} > {log} 2>&1"
