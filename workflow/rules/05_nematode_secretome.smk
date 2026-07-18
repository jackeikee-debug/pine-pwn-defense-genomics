import os


rule predict_nematode_secretome:
    output:
        "results/tables/nematode_secretome_candidates.tsv"
    conda:
        "../envs/annotation.yaml"
    log:
        "results/logs/nematode_secretome_candidates.log"
    params:
        proteins="data/raw/nematode/Bursaphelenchus_xylophilus.BXYJv5.protein.faa.gz",
        min_length=config["secretome"]["min_protein_length"],
        small_secreted_max_length=config["secretome"]["small_secreted_max_length"],
        cysteine_rich_threshold=config["secretome"]["cysteine_rich_threshold"],
        max_transmembrane_domains=config["secretome"]["max_transmembrane_domains"]
    shell:
        "python scripts/predict_secretome.py "
        "--proteins {params.proteins} "
        "--output {output} "
        "--min-length {params.min_length} "
        "--small-secreted-max-length {params.small_secreted_max_length} "
        "--cysteine-rich-threshold {params.cysteine_rich_threshold} "
        "--max-transmembrane-domains {params.max_transmembrane_domains} > {log} 2>&1"


rule select_nematode_secretome_swissprot_query:
    input:
        secretome="results/tables/nematode_secretome_candidates.tsv"
    output:
        fasta="data/interim/nematode_secretome_swissprot_query.faa",
        manifest="results/tables/nematode_secretome_swissprot_query_manifest.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/nematode_secretome_swissprot_query.log"
    params:
        proteins="data/raw/nematode/Bursaphelenchus_xylophilus.BXYJv5.protein.faa.gz"
    shell:
        "python scripts/select_fasta_by_table.py "
        "--fasta {params.proteins} "
        "--table {input.secretome} "
        "--id-column protein_id "
        "--where-column secretome_candidate "
        "--where-value yes "
        "--output-fasta {output.fasta} "
        "--manifest {output.manifest} > {log} 2>&1"


rule run_nematode_secretome_swissprot_diamond:
    input:
        query="data/interim/nematode_secretome_swissprot_query.faa",
        db="data/external/swissprot/uniprot_sprot.dmnd"
    output:
        "data/interim/nematode_secretome_swissprot_diamond.tsv"
    threads:
        config["resources"]["threads"]
    conda:
        "../envs/annotation.yaml"
    log:
        "results/logs/nematode_secretome_swissprot_diamond.log"
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


rule parse_nematode_secretome_swissprot_annotations:
    input:
        "data/interim/nematode_secretome_swissprot_diamond.tsv"
    output:
        annotations="data/interim/nematode_secretome_swissprot_annotations.tsv",
        hits="results/tables/nematode_secretome_swissprot_hits.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/nematode_secretome_swissprot_parse.log"
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


rule classify_nematode_candidate_effectors:
    input:
        secretome="results/tables/nematode_secretome_candidates.tsv",
        keywords=config["metadata"]["effector_keywords"],
        annotations="data/interim/nematode_secretome_swissprot_annotations.tsv"
    output:
        "results/tables/nematode_candidate_effectors.tsv"
    conda:
        "../envs/annotation.yaml"
    log:
        "results/logs/nematode_candidate_effectors.log"
    shell:
        "python scripts/classify_effectors.py "
        "--secretome {input.secretome} "
        "--keywords {input.keywords} "
        "--annotations {input.annotations} "
        "--output {output} > {log} 2>&1"


rule build_full_nematode_secretome_evidence_audit:
    input:
        heuristic="results/tables/nematode_secretome_candidates.tsv",
        deepsig="data/interim/nematode_full_proteome_deepsig.gff3",
        tmhmm="data/interim/nematode_full_proteome_tmhmm.tsv"
    output:
        table="results/tables/nematode_full_secretome_sequence_baseline.tsv",
        report="results/reports/nematode_full_secretome_sequence_baseline.md"
    conda:
        "../envs/secretome.yaml"
    log:
        "results/logs/nematode_full_secretome_sequence_baseline.log"
    benchmark:
        "results/benchmarks/nematode_full_secretome_sequence_baseline.tsv"
    params:
        min_signal_score=0.80,
        signal_overlap_tolerance=5
    shell:
        "python scripts/build_secretome_evidence_audit.py "
        "--heuristic {input.heuristic} "
        "--deepsig {input.deepsig} "
        "--tmhmm {input.tmhmm} "
        "--output {output.table} "
        "--report {output.report} "
        "--min-signal-score {params.min_signal_score} "
        "--signal-overlap-tolerance {params.signal_overlap_tolerance} "
        "> {log} 2>&1"


rule extract_standard_nematode_secretome_fasta:
    input:
        source="data/raw/nematode/Bursaphelenchus_xylophilus.BXYJv5.protein.faa.gz",
        audit="results/tables/nematode_full_secretome_sequence_baseline.tsv"
    output:
        fasta="data/interim/nematode_standard_secreted_soluble.faa",
        manifest="results/tables/nematode_standard_secretome_sequence_manifest.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/nematode_standard_secretome_extract.log"
    benchmark:
        "results/benchmarks/nematode_standard_secretome_extract.tsv"
    shell:
        "python scripts/extract_standard_secretome_fasta.py "
        "--source-fasta {input.source} "
        "--audit-table {input.audit} "
        "--output-fasta {output.fasta} "
        "--manifest {output.manifest} "
        "--log {log}"


rule download_secretome_functional_databases:
    input:
        config="config/secretome_functional_databases.tsv"
    output:
        merops="data/external/functional_databases/merops_pepunit_2023-02-22.faa",
        dbcan="data/external/functional_databases/dbCAN_db_v5-2_9-13-2025.hmm",
        manifest="results/tables/secretome_functional_database_manifest.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/secretome_functional_database_download.log"
    benchmark:
        "results/benchmarks/secretome_functional_database_download.tsv"
    shell:
        "python scripts/download_secretome_functional_databases.py "
        "--config {input.config} "
        "--output-dir data/external/functional_databases "
        "--manifest {output.manifest} "
        "--log {log}"


rule normalize_merops_peptidase_unit_database:
    input:
        "data/external/functional_databases/merops_pepunit_2023-02-22.faa"
    output:
        fasta="data/interim/secretome_functional_searches/merops_pepunit_2023-02-22.normalized.faa",
        summary="results/tables/merops_pepunit_normalization.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/merops_pepunit_normalization.log"
    benchmark:
        "results/benchmarks/merops_pepunit_normalization.tsv"
    shell:
        "python scripts/normalize_database_fasta.py "
        "--input {input} --output {output.fasta} --summary {output.summary} "
        "--input-encoding cp1252 --log {log}"


rule build_merops_peptidase_unit_diamond_database:
    input:
        "data/interim/secretome_functional_searches/merops_pepunit_2023-02-22.normalized.faa"
    output:
        "data/interim/secretome_functional_searches/merops_pepunit_2023-02-22.dmnd"
    conda:
        "../envs/annotation.yaml"
    log:
        "results/logs/merops_diamond_makedb.log"
    benchmark:
        "results/benchmarks/merops_diamond_makedb.tsv"
    shell:
        "diamond makedb --in {input} "
        "--db data/interim/secretome_functional_searches/merops_pepunit_2023-02-22 "
        "> {log} 2>&1"


rule press_dbcan_hmm_database:
    input:
        "data/external/functional_databases/dbCAN_db_v5-2_9-13-2025.hmm"
    output:
        h3f="data/external/functional_databases/dbCAN_db_v5-2_9-13-2025.hmm.h3f",
        h3i="data/external/functional_databases/dbCAN_db_v5-2_9-13-2025.hmm.h3i",
        h3m="data/external/functional_databases/dbCAN_db_v5-2_9-13-2025.hmm.h3m",
        h3p="data/external/functional_databases/dbCAN_db_v5-2_9-13-2025.hmm.h3p"
    conda:
        "../envs/annotation.yaml"
    log:
        "results/logs/dbcan_hmmpress.log"
    benchmark:
        "results/benchmarks/dbcan_hmmpress.tsv"
    shell:
        "hmmpress -f {input} > {log} 2>&1"


rule search_standard_secretome_merops:
    input:
        query="data/interim/nematode_standard_secreted_soluble.faa",
        database="data/interim/secretome_functional_searches/merops_pepunit_2023-02-22.dmnd"
    output:
        "data/interim/secretome_functional_searches/merops_standard_secretome.tsv"
    threads:
        config["resources"]["threads"]
    conda:
        "../envs/annotation.yaml"
    log:
        "results/logs/merops_standard_secretome_diamond.log"
    benchmark:
        "results/benchmarks/merops_standard_secretome_diamond.tsv"
    shell:
        "diamond blastp --query {input.query} --db {input.database} --out {output} "
        "--outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore qcovhsp scovhsp stitle "
        "--max-target-seqs 25 --evalue 1e-5 --sensitive --threads {threads} "
        "> {log} 2>&1"


rule search_standard_secretome_dbcan:
    input:
        query="data/interim/nematode_standard_secreted_soluble.faa",
        database="data/external/functional_databases/dbCAN_db_v5-2_9-13-2025.hmm",
        h3f="data/external/functional_databases/dbCAN_db_v5-2_9-13-2025.hmm.h3f",
        h3i="data/external/functional_databases/dbCAN_db_v5-2_9-13-2025.hmm.h3i",
        h3m="data/external/functional_databases/dbCAN_db_v5-2_9-13-2025.hmm.h3m",
        h3p="data/external/functional_databases/dbCAN_db_v5-2_9-13-2025.hmm.h3p"
    output:
        "data/interim/secretome_functional_searches/dbcan_standard_secretome.domtblout"
    threads:
        config["resources"]["threads"]
    conda:
        "../envs/annotation.yaml"
    log:
        "results/logs/dbcan_standard_secretome_hmmscan.log"
    benchmark:
        "results/benchmarks/dbcan_standard_secretome_hmmscan.tsv"
    shell:
        "hmmscan --cpu {threads} --domtblout {output} --noali "
        "{input.database} {input.query} > {log} 2>&1"


rule integrate_standard_secretome_functional_evidence:
    input:
        manifest="results/tables/nematode_standard_secretome_sequence_manifest.tsv",
        interpro="data/interim/nematode_standard_secretome_interproscan.tsv",
        merops="data/interim/secretome_functional_searches/merops_standard_secretome.tsv",
        dbcan="data/interim/secretome_functional_searches/dbcan_standard_secretome.domtblout"
    output:
        interpro_hits="results/tables/nematode_secretome_interpro_hits.tsv",
        functional="results/tables/nematode_secretome_functional_evidence.tsv",
        report="results/reports/nematode_secretome_functional_evidence.md"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/nematode_secretome_functional_evidence.log"
    benchmark:
        "results/benchmarks/nematode_secretome_functional_evidence.tsv"
    shell:
        "python scripts/integrate_secretome_functional_evidence.py "
        "--manifest {input.manifest} --interpro {input.interpro} "
        "--merops {input.merops} --dbcan {input.dbcan} "
        "--interpro-hits {output.interpro_hits} --output {output.functional} "
        "--report {output.report} --log {log}"


rule download_published_pwn_secretome_sources:
    input:
        "config/published_secretome_downloads.tsv"
    output:
        shinya_table="data/external/published_pwn_secretomes/PMC3689755_supplementary/pone.0067377.s005.xls",
        shinya_fasta="data/external/published_pwn_secretomes/PMC3689755_supplementary/pone.0067377.s011.txt",
        cardoso_table="data/external/published_pwn_secretomes/PMC5150578_supplementary/srep39007-s2.xls",
        silva_table="data/external/published_pwn_secretomes/PMC8144518_supplementary/Data_Sheet_1.xlsx",
        wbps_proteome="data/external/published_pwn_secretomes/bursaphelenchus_xylophilus.PRJEA64437.WBPS19.protein.fa.gz",
        manifest="results/tables/published_secretome_download_manifest.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/published_secretome_download.log"
    benchmark:
        "results/benchmarks/published_secretome_download.tsv"
    shell:
        "python scripts/download_published_secretome_sources.py "
        "--config {input} --manifest {output.manifest} --log {log}"


rule build_published_pwn_secretome_queries:
    input:
        sources="config/published_secretome_sources.tsv",
        wbps_proteome="data/external/published_pwn_secretomes/bursaphelenchus_xylophilus.PRJEA64437.WBPS19.protein.fa.gz",
        shinya_fasta="data/external/published_pwn_secretomes/PMC3689755_supplementary/pone.0067377.s011.txt",
        shinya_table="data/external/published_pwn_secretomes/PMC3689755_supplementary/pone.0067377.s005.xls",
        cardoso_table="data/external/published_pwn_secretomes/PMC5150578_supplementary/srep39007-s2.xls",
        silva_table="data/external/published_pwn_secretomes/PMC8144518_supplementary/Data_Sheet_1.xlsx"
    output:
        records="results/tables/published_pwn_secretome_source_records.tsv",
        query="data/interim/published_pwn_secretome_queries.faa",
        manifest="results/tables/published_pwn_secretome_query_manifest.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/published_pwn_secretome_query_build.log"
    benchmark:
        "results/benchmarks/published_pwn_secretome_query_build.tsv"
    shell:
        "python scripts/build_published_secretome_queries.py "
        "--sources {input.sources} --wbps-proteome {input.wbps_proteome} "
        "--legacy-secretome-fasta {input.shinya_fasta} "
        "--records-output {output.records} --query-fasta {output.query} "
        "--manifest-output {output.manifest} --log {log}"


rule build_bxyjv5_published_secretome_diamond_database:
    input:
        "data/raw/nematode/Bursaphelenchus_xylophilus.BXYJv5.protein.faa.gz"
    output:
        "data/interim/published_pwn_secretome_bxyjv5.dmnd"
    conda:
        "../envs/annotation.yaml"
    log:
        "results/logs/published_pwn_secretome_diamond_makedb.log"
    benchmark:
        "results/benchmarks/published_pwn_secretome_diamond_makedb.tsv"
    shell:
        "diamond makedb --in {input} --db data/interim/published_pwn_secretome_bxyjv5 "
        "> {log} 2>&1"


rule map_published_pwn_secretomes_to_bxyjv5:
    input:
        query="data/interim/published_pwn_secretome_queries.faa",
        database="data/interim/published_pwn_secretome_bxyjv5.dmnd"
    output:
        "data/interim/published_pwn_secretome_to_bxyjv5.tsv"
    threads:
        config["resources"]["threads"]
    conda:
        "../envs/annotation.yaml"
    log:
        "results/logs/published_pwn_secretome_diamond.log"
    benchmark:
        "results/benchmarks/published_pwn_secretome_diamond.tsv"
    shell:
        "diamond blastp --query {input.query} --db {input.database} --out {output} "
        "--outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore qcovhsp scovhsp "
        "--max-target-seqs 10 --evalue 1e-10 --very-sensitive --threads {threads} "
        "> {log} 2>&1"


rule integrate_published_pwn_secretome_crossvalidation:
    input:
        manifest="results/tables/published_pwn_secretome_query_manifest.tsv",
        alignments="data/interim/published_pwn_secretome_to_bxyjv5.tsv",
        records="results/tables/published_pwn_secretome_source_records.tsv",
        functional="results/tables/nematode_secretome_functional_evidence.tsv"
    output:
        mapping="results/tables/published_pwn_secretome_bxyjv5_mapping.tsv",
        source_mappings="results/tables/published_pwn_secretome_source_mappings.tsv",
        integrated="results/tables/nematode_secretome_integrated_evidence.tsv",
        report="results/reports/published_pwn_secretome_crossvalidation.md"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/published_pwn_secretome_crossvalidation.log"
    benchmark:
        "results/benchmarks/published_pwn_secretome_crossvalidation.tsv"
    shell:
        "python scripts/integrate_published_secretome_crossvalidation.py "
        "--manifest {input.manifest} --alignments {input.alignments} "
        "--source-records {input.records} --functional {input.functional} "
        "--mapping-output {output.mapping} --source-mappings-output {output.source_mappings} "
        "--output {output.integrated} --report {output.report} --log {log}"


rule build_functionally_annotated_nematode_secretome_audit:
    input:
        heuristic="results/tables/nematode_secretome_candidates.tsv",
        deepsig="data/interim/nematode_full_proteome_deepsig.gff3",
        tmhmm="data/interim/nematode_full_proteome_tmhmm.tsv",
        functional="results/tables/nematode_secretome_integrated_evidence.tsv"
    output:
        table="results/tables/nematode_full_secretome_evidence.tsv",
        report="results/reports/nematode_full_secretome_evidence.md"
    conda:
        "../envs/secretome.yaml"
    log:
        "results/logs/nematode_full_secretome_evidence.log"
    benchmark:
        "results/benchmarks/nematode_full_secretome_evidence.tsv"
    params:
        min_signal_score=0.80,
        signal_overlap_tolerance=5
    shell:
        "python scripts/build_secretome_evidence_audit.py "
        "--heuristic {input.heuristic} --deepsig {input.deepsig} "
        "--tmhmm {input.tmhmm} --functional {input.functional} "
        "--output {output.table} --report {output.report} "
        "--min-signal-score {params.min_signal_score} "
        "--signal-overlap-tolerance {params.signal_overlap_tolerance} "
        "> {log} 2>&1"
