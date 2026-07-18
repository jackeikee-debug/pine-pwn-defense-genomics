PATHVIEW_PATHWAYS = ["04626", "04016", "04075", "00940", "00480"]
PATHVIEW_CONTRASTS = [
    "genotype_by_inoculum",
    "resistant_pwn_vs_water",
    "susceptible_pwn_vs_water",
]
PATHVIEW_FIGURES = expand(
    "results/figures/pathview/ath{pathway}.{contrast}.png",
    pathway=PATHVIEW_PATHWAYS,
    contrast=PATHVIEW_CONTRASTS,
)


rule download_pathview_arabidopsis_reference:
    output:
        fasta="data/external/pathview/arath_reviewed_UP000006548.fasta",
        metadata="data/external/pathview/arath_reviewed_UP000006548.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/pathview_arabidopsis_download.log"
    benchmark:
        "results/benchmarks/pathview_arabidopsis_download.tsv"
    shell:
        "curl -L --fail --retry 3 "
        "'https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=fasta&query=%28proteome%3AUP000006548%29%20AND%20%28reviewed%3Atrue%29' "
        "-o {output.fasta} && "
        "curl -L --fail --retry 3 "
        "'https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=tsv&fields=accession%2Cid%2Cgene_names%2Cgene_oln%2Clength&query=%28proteome%3AUP000006548%29%20AND%20%28reviewed%3Atrue%29' "
        "-o {output.metadata} > {log} 2>&1"


rule build_pathview_diamond_databases:
    input:
        arath="data/external/pathview/arath_reviewed_UP000006548.fasta",
        pmas="data/processed/proteomes_longest/pmas.faa"
    output:
        arath_db="data/external/pathview/arath_reviewed_UP000006548.dmnd",
        pmas_db="data/interim/pathview/pmas_longest.dmnd"
    conda:
        "../envs/orthology.yaml"
    log:
        "results/logs/pathview_diamond_databases.log"
    benchmark:
        "results/benchmarks/pathview_diamond_databases.tsv"
    shell:
        "diamond makedb --in {input.arath} --db data/external/pathview/arath_reviewed_UP000006548 "
        "> {log} 2>&1 && "
        "diamond makedb --in {input.pmas} --db data/interim/pathview/pmas_longest "
        ">> {log} 2>&1"


rule search_pathview_reciprocal_hits:
    input:
        arath_db="data/external/pathview/arath_reviewed_UP000006548.dmnd",
        pmas_db="data/interim/pathview/pmas_longest.dmnd",
        arath="data/external/pathview/arath_reviewed_UP000006548.fasta",
        pmas="data/processed/proteomes_longest/pmas.faa"
    output:
        forward="data/interim/pathview/pmas_to_arath.tsv",
        reverse_hits="data/interim/pathview/arath_to_pmas.tsv"
    conda:
        "../envs/orthology.yaml"
    log:
        "results/logs/pathview_reciprocal_diamond.log"
    benchmark:
        "results/benchmarks/pathview_reciprocal_diamond.tsv"
    shell:
        "diamond blastp --query {input.pmas} --db {input.arath_db} --out {output.forward} "
        "--sensitive --evalue 1e-5 --max-target-seqs 25 "
        "--outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore "
        "> {log} 2>&1 && "
        "diamond blastp --query {input.arath} --db {input.pmas_db} --out {output.reverse_hits} "
        "--sensitive --evalue 1e-5 --max-target-seqs 25 "
        "--outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore "
        ">> {log} 2>&1"


rule build_pathview_reciprocal_crosswalk:
    input:
        forward="data/interim/pathview/pmas_to_arath.tsv",
        reverse_hits="data/interim/pathview/arath_to_pmas.tsv",
        metadata="data/external/pathview/arath_reviewed_UP000006548.tsv",
        pmas="data/processed/proteomes_longest/pmas.faa"
    output:
        accepted="results/tables/pmas_arath_reciprocal_best_hits.tsv",
        audit="results/tables/pmas_arath_reciprocal_best_hits_audit.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/pathview_reciprocal_crosswalk.log"
    benchmark:
        "results/benchmarks/pathview_reciprocal_crosswalk.tsv"
    shell:
        "python scripts/build_reciprocal_best_hits.py "
        "--forward-hits {input.forward} --reverse-hits {input.reverse_hits} "
        "--uniprot-metadata {input.metadata} --pmas-fasta {input.pmas} "
        "--output {output.accepted} --audit-output {output.audit} > {log} 2>&1"


rule build_pmas_pathview_inputs:
    input:
        resistant="results/tables/effector_target_network_pmas_GigaDB102688_DESeq2_resistant_pwn_vs_water.tsv",
        susceptible="results/tables/effector_target_network_pmas_GigaDB102688_DESeq2_susceptible_pwn_vs_water.tsv",
        interaction="results/tables/effector_target_network_pmas_GigaDB102688_DESeq2_genotype_by_inoculum.tsv",
        crosswalk="results/tables/pmas_arath_reciprocal_best_hits.tsv"
    output:
        "results/tables/pmas_pathview_gene_values.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/pmas_pathview_inputs.log"
    benchmark:
        "results/benchmarks/pmas_pathview_inputs.tsv"
    shell:
        "python scripts/build_pmas_pathview_inputs.py "
        "--resistant {input.resistant} --susceptible {input.susceptible} "
        "--interaction {input.interaction} --crosswalk {input.crosswalk} "
        "--output {output} > {log} 2>&1"


rule render_pmas_pathview:
    input:
        values="results/tables/pmas_pathview_gene_values.tsv",
        pathways="config/pmas_pathview_pathways.tsv"
    output:
        coverage="results/tables/pmas_pathview_coverage.tsv",
        membership="results/tables/pmas_pathview_pathway_membership.tsv",
        figures=PATHVIEW_FIGURES
    conda:
        "../envs/r.yaml"
    log:
        "results/logs/pmas_pathview_render.log"
    benchmark:
        "results/benchmarks/pmas_pathview_render.tsv"
    shell:
        "Rscript scripts/render_pmas_pathview.R "
        "--gene-values {input.values} --pathways {input.pathways} "
        "--membership-output {output.membership} --coverage-output {output.coverage} "
        "--output-dir results/figures/pathview > {log} 2>&1"
