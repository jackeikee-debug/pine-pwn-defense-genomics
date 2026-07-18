rule validate_metadata:
    input:
        config="config/config.yaml",
        species="config/species_metadata.tsv",
        phenotype="config/phenotype_matrix.tsv"
    output:
        touch("results/reports/metadata.validated.txt")
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/validate_metadata.log"
    shell:
        "python scripts/validate_metadata.py --config {input.config} > {log} 2>&1"


rule clean_proteome_headers:
    input:
        lambda wildcards: PROTEOME_FASTA_BY_SPECIES[wildcards.species_id]
    output:
        fasta="data/processed/proteomes/{species_id}.faa",
        summary="results/tables/proteome_cleaning/{species_id}.summary.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/prepare_proteomes/{species_id}.log"
    params:
        min_length=config["orthology"]["min_protein_length"]
    shell:
        "python scripts/clean_fasta_headers.py "
        "--input {input} "
        "--output {output.fasta} "
        "--species-id {wildcards.species_id} "
        "--min-length {params.min_length} "
        "--summary {output.summary} "
        "--log {log}"


rule prepare_longest_isoforms:
    input:
        "data/processed/proteomes/{species_id}.faa"
    output:
        fasta="data/processed/proteomes_longest/{species_id}.faa",
        summary="results/tables/proteome_longest_isoforms/{species_id}.summary.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/prepare_longest_isoforms/{species_id}.log"
    shell:
        "python scripts/prepare_longest_isoforms.py "
        "--input {input} "
        "--output {output.fasta} "
        "--summary {output.summary} "
        "--log {log}"


rule summarize_proteomes:
    input:
        expand("data/processed/proteomes/{species_id}.faa", species_id=PROTEOME_SPECIES)
    output:
        "results/tables/species_summary.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/species_summary.log"
    shell:
        "python scripts/summarize_proteomes.py --inputs {input} --output {output} --log {log}"


rule summarize_longest_isoform_proteomes:
    input:
        expand("data/processed/proteomes_longest/{species_id}.faa", species_id=PROTEOME_SPECIES)
    output:
        "results/tables/species_summary_longest_isoforms.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/species_summary_longest_isoforms.log"
    shell:
        "python scripts/summarize_proteomes.py --inputs {input} --output {output} --log {log}"


rule build_host_proteome_quality_audit:
    input:
        metadata="config/species_metadata.tsv",
        summary="results/tables/species_summary_longest_isoforms.tsv",
        busco="results/tables/busco_proteome_summary.tsv"
    output:
        table="results/tables/manuscript_proteome_quality_audit.tsv",
        report="results/reports/manuscript_proteome_quality_audit.md"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/manuscript_proteome_quality_audit.log"
    benchmark:
        "results/benchmarks/manuscript_proteome_quality_audit.tsv"
    shell:
        "python scripts/build_proteome_quality_audit.py "
        "--metadata {input.metadata} "
        "--summary {input.summary} "
        "--busco {input.busco} "
        "--output {output.table} "
        "--report {output.report} > {log} 2>&1"


rule clean_expansion_proteome_headers:
    input:
        lambda wildcards: EXPANSION_FASTA_BY_SPECIES[wildcards.species_id]
    output:
        fasta="data/processed/exp_proteomes/{species_id}.faa",
        summary="results/tables/exp_clean/{species_id}.summary.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/exp_prepare/{species_id}.log"
    params:
        min_length=config["orthology"]["min_protein_length"]
    shell:
        "python scripts/clean_fasta_headers.py "
        "--input {input} "
        "--output {output.fasta} "
        "--species-id {wildcards.species_id} "
        "--min-length {params.min_length} "
        "--summary {output.summary} "
        "--log {log}"


rule prepare_expansion_longest_isoforms:
    input:
        "data/processed/exp_proteomes/{species_id}.faa"
    output:
        fasta="data/processed/exp_longest/{species_id}.faa",
        summary="results/tables/exp_longest/{species_id}.summary.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/exp_longest/{species_id}.log"
    shell:
        "python scripts/prepare_longest_isoforms.py "
        "--input {input} "
        "--output {output.fasta} "
        "--summary {output.summary} "
        "--log {log}"


rule summarize_expansion_longest_isoform_proteomes:
    input:
        expand("data/processed/exp_longest/{species_id}.faa", species_id=EXPANSION_SPECIES)
    output:
        "results/tables/exp_species_summary.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/exp_species_summary.log"
    shell:
        "python scripts/summarize_proteomes.py --inputs {input} --output {output} --log {log}"
