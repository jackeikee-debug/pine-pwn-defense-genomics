rule annotate_defense_modules:
    input:
        stable_genes="results/tables/stable_orthogroup_genes.tsv",
        annotations=config["metadata"]["protein_function_annotations"],
        modules=config["metadata"]["defense_modules"]
    output:
        matrix="results/tables/orthogroup_defense_matrix.tsv",
        gene_hits="results/tables/defense_module_gene_hits.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/defense_module_annotation.log"
    shell:
        "python scripts/annotate_defense_modules.py "
        "--stable-genes {input.stable_genes} "
        "--annotations {input.annotations} "
        "--modules {input.modules} "
        "--matrix-output {output.matrix} "
        "--gene-hits-output {output.gene_hits} > {log} 2>&1"


rule summarize_defense_module_hits:
    input:
        matrix="results/tables/orthogroup_defense_matrix.tsv",
        phenotype=config["metadata"]["phenotype"]
    output:
        species="results/tables/defense_module_species_summary.tsv",
        modules="results/tables/defense_module_summary.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/defense_module_summary.log"
    shell:
        "python scripts/summarize_defense_module_hits.py "
        "--matrix {input.matrix} "
        "--phenotype {input.phenotype} "
        "--species-output {output.species} "
        "--module-output {output.modules} > {log} 2>&1"


rule prioritize_defense_module_candidates:
    input:
        matrix="results/tables/orthogroup_defense_matrix.tsv",
        gene_hits="results/tables/defense_module_gene_hits.tsv",
        diamond_hits="results/tables/swissprot_diamond_hits.tsv",
        stable_orthogroups="results/tables/stable_orthogroups.tsv",
        phenotype=config["metadata"]["phenotype"]
    output:
        "results/tables/defense_module_candidate_priority.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/defense_module_candidate_priority.log"
    shell:
        "python scripts/prioritize_defense_module_candidates.py "
        "--matrix {input.matrix} "
        "--gene-hits {input.gene_hits} "
        "--diamond-hits {input.diamond_hits} "
        "--stable-orthogroups {input.stable_orthogroups} "
        "--phenotype {input.phenotype} "
        "--output {output} > {log} 2>&1"


rule shortlist_defense_module_candidates:
    input:
        "results/tables/defense_module_candidate_priority.tsv"
    output:
        "results/tables/defense_module_candidate_shortlist.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/defense_module_candidate_shortlist.log"
    params:
        max_per_module_direction=config["candidate_orthogroups"]["shortlist_max_per_module_direction"],
        min_mean_bitscore=config["candidate_orthogroups"]["shortlist_min_mean_bitscore"],
        min_core_retention=config["candidate_orthogroups"]["shortlist_min_core_retention"],
        min_core_jaccard=config["candidate_orthogroups"]["shortlist_min_core_jaccard"],
        allowed_tiers=config["candidate_orthogroups"]["shortlist_allowed_tiers"]
    shell:
        "python scripts/shortlist_defense_module_candidates.py "
        "--priority {input} "
        "--output {output} "
        "--max-per-module-direction {params.max_per_module_direction} "
        "--min-mean-bitscore {params.min_mean_bitscore} "
        "--min-core-retention {params.min_core_retention} "
        "--min-core-jaccard {params.min_core_jaccard} "
        "--allowed-tiers {params.allowed_tiers} > {log} 2>&1"
