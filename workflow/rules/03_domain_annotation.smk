import os


rule download_swissprot_fasta:
    output:
        "data/external/swissprot/uniprot_sprot.fasta.gz"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/download_swissprot.log"
    params:
        url=config["annotation"]["swissprot_url"]
    shell:
        "python -c \"from pathlib import Path; import urllib.request; "
        "out=Path(r'{output}'); out.parent.mkdir(parents=True, exist_ok=True); "
        "urllib.request.urlretrieve(r'{params.url}', out); "
        "Path(r'{log}').parent.mkdir(parents=True, exist_ok=True); "
        "Path(r'{log}').write_text('Downloaded SwissProt FASTA from {params.url}\\n', encoding='utf-8')\""


rule build_swissprot_diamond_db:
    input:
        "data/external/swissprot/uniprot_sprot.fasta.gz"
    output:
        "data/external/swissprot/uniprot_sprot.dmnd"
    conda:
        "../envs/annotation.yaml"
    log:
        "results/logs/build_swissprot_diamond_db.log"
    params:
        db_prefix=lambda wildcards, output: os.path.splitext(output[0])[0]
    shell:
        "diamond makedb --in {input} -d {params.db_prefix} > {log} 2>&1"


rule select_swissprot_annotation_query_batch:
    input:
        stable_genes="results/tables/stable_orthogroup_genes.tsv",
        proteomes=expand("data/processed/proteomes_longest/{species_id}.faa", species_id=PROTEOME_SPECIES)
    output:
        fasta="data/interim/swissprot_query_batch.faa",
        manifest="results/tables/swissprot_query_batch_manifest.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/select_swissprot_query_batch.log"
    params:
        proteome_dir=lambda wildcards, input: os.path.dirname(input.proteomes[0]),
        max_total=config["annotation"]["swissprot_query_max_total"],
        max_per_species=config["annotation"]["swissprot_query_max_per_species"]
    shell:
        "python scripts/select_annotation_query_batch.py "
        "--stable-genes {input.stable_genes} "
        "--proteome-dir {params.proteome_dir} "
        "--output-fasta {output.fasta} "
        "--manifest {output.manifest} "
        "--max-total {params.max_total} "
        "--max-per-species {params.max_per_species} > {log} 2>&1"


rule run_swissprot_diamond_batch:
    input:
        query="data/interim/swissprot_query_batch.faa",
        db="data/external/swissprot/uniprot_sprot.dmnd"
    output:
        "data/interim/swissprot_diamond_batch.tsv"
    threads:
        config["resources"]["threads"]
    conda:
        "../envs/annotation.yaml"
    log:
        "results/logs/run_swissprot_diamond_batch.log"
    params:
        db_prefix=lambda wildcards, input: os.path.splitext(input.db)[0]
    shell:
        "diamond blastp "
        "--query {input.query} "
        "--db {params.db_prefix} "
        "--out {output} "
        "--outfmt 6 qseqid sseqid pident length evalue bitscore stitle "
        "--max-target-seqs 1 "
        "--threads {threads} > {log} 2>&1"


rule parse_swissprot_diamond_batch:
    input:
        "data/interim/swissprot_diamond_batch.tsv"
    output:
        annotations=config["metadata"]["protein_function_annotations"],
        hits="results/tables/swissprot_diamond_hits.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/parse_swissprot_diamond_batch.log"
    params:
        min_pident=config["annotation"]["swissprot_min_pident"],
        min_bitscore=config["annotation"]["swissprot_min_bitscore"]
    shell:
        "python scripts/parse_diamond_annotations.py "
        "--diamond {input} "
        "--annotation-output {output.annotations} "
        "--hits-output {output.hits} "
        "--min-pident {params.min_pident} "
        "--min-bitscore {params.min_bitscore} > {log} 2>&1"
