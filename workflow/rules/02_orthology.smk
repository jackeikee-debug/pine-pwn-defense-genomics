import os
import sys


rule prepare_orthofinder_inputs:
    input:
        proteomes=expand("data/processed/proteomes_longest/{species_id}.faa", species_id=PROTEOME_SPECIES)
    output:
        input_dir=directory("results/orthofinder/pilot_input"),
        manifest="results/tables/orthofinder_pilot_input_manifest.tsv"
    conda:
        "../envs/orthology.yaml"
    log:
        "results/logs/orthofinder_pilot_input.log"
    params:
        source_dir=lambda wildcards, input: os.path.dirname(input.proteomes[0]),
        species=PROTEOME_SPECIES
    shell:
        "python scripts/prepare_orthofinder_inputs.py "
        "--source-dir {params.source_dir} "
        "--output-dir {output.input_dir} "
        "--manifest {output.manifest} "
        "--species {params.species} > {log} 2>&1"


rule run_orthofinder_pilot:
    input:
        input_dir="results/orthofinder/pilot_input",
        manifest="results/tables/orthofinder_pilot_input_manifest.tsv"
    output:
        done="results/orthofinder/pilot_run.done"
    threads:
        config["resources"]["threads"]
    conda:
        "../envs/orthology.yaml"
    log:
        "results/logs/orthofinder_pilot_run.log"
    params:
        outdir=lambda wildcards, output: os.path.splitext(output.done)[0]
    shell:
        "orthofinder -f {input.input_dir} -S diamond -t {threads} -a {threads} -o {params.outdir} > {log} 2>&1 && "
        "python -c \"from pathlib import Path; Path(r'{output.done}').write_text('done\\n', encoding='utf-8')\""


rule collect_orthofinder_expansion_longest:
    input:
        core=expand("data/processed/proteomes_longest/{species_id}.faa", species_id=PROTEOME_SPECIES),
        expansion=expand("data/processed/exp_longest/{species_id}.faa", species_id=EXPANSION_SPECIES)
    output:
        directory("results/orthofinder/exp_source")
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/of_exp_collect.log"
    shell:
        "python -c \"from pathlib import Path; import shutil; "
        "out=Path(r'{output}'); out.mkdir(parents=True, exist_ok=True); "
        "[p.unlink() for p in out.glob('*.faa')]; "
        "[shutil.copy2(Path(p), out / Path(p).name) for p in r'{input.core} {input.expansion}'.split()]; "
        "Path(r'{log}').parent.mkdir(parents=True, exist_ok=True); Path(r'{log}').write_text('Collected expansion OrthoFinder source FASTA files.\\n', encoding='utf-8')\""


rule prepare_orthofinder_expansion_inputs:
    input:
        source_dir="results/orthofinder/exp_source"
    output:
        input_dir=directory("results/orthofinder/exp_input"),
        manifest="results/tables/of_expansion_manifest.tsv"
    conda:
        "../envs/orthology.yaml"
    log:
        "results/logs/of_exp_input.log"
    params:
        species=ORTHOFINDER_EXPANSION_SPECIES
    shell:
        "python scripts/prepare_orthofinder_inputs.py "
        "--source-dir {input.source_dir} "
        "--output-dir {output.input_dir} "
        "--manifest {output.manifest} "
        "--species {params.species} > {log} 2>&1"


rule run_orthofinder_expansion:
    input:
        input_dir="results/orthofinder/exp_input",
        manifest="results/tables/of_expansion_manifest.tsv"
    output:
        done="results/orthofinder/exp_run.done"
    threads:
        config["resources"]["threads"]
    conda:
        "../envs/orthology.yaml"
    log:
        "results/logs/of_exp_run.log"
    params:
        outdir=lambda wildcards, output: os.path.splitext(output.done)[0]
    shell:
        "orthofinder -f {input.input_dir} -S diamond -t {threads} -a {threads} -o {params.outdir} > {log} 2>&1 && "
        "python -c \"from pathlib import Path; Path(r'{output.done}').write_text('done\\n', encoding='utf-8')\""


rule compare_orthofinder_pilot_expansion:
    input:
        pilot="results/orthofinder/pilot_run/Results_Jul03/Orthogroups/Orthogroups.tsv",
        expansion="results/orthofinder/exp_run/Results_Jul03/Orthogroups/Orthogroups.tsv",
        pilot_done="results/orthofinder/pilot_run.done",
        expansion_done="results/orthofinder/exp_run.done"
    output:
        stability="results/tables/of_pilot_vs_expansion_stability.tsv",
        summary="results/tables/of_pilot_vs_expansion_summary.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/of_pilot_vs_expansion.log"
    params:
        core_species=PROTEOME_SPECIES
    shell:
        "python scripts/compare_orthofinder_runs.py "
        "--pilot-orthogroups {input.pilot} "
        "--expansion-orthogroups {input.expansion} "
        "--core-species {params.core_species} "
        "--output {output.stability} "
        "--summary {output.summary} > {log} 2>&1"


rule build_stable_orthogroup_tables:
    input:
        orthogroups="results/orthofinder/pilot_run/Results_Jul03/Orthogroups/Orthogroups.tsv",
        stability="results/tables/of_pilot_vs_expansion_stability.tsv"
    output:
        orthogroups="results/tables/stable_orthogroups.tsv",
        genes="results/tables/stable_orthogroup_genes.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/stable_orthogroup_tables.log"
    params:
        core_species=PROTEOME_SPECIES
    shell:
        "python scripts/build_stable_orthogroup_tables.py "
        "--orthogroups {input.orthogroups} "
        "--stability {input.stability} "
        "--species {params.core_species} "
        "--orthogroup-output {output.orthogroups} "
        "--gene-output {output.genes} > {log} 2>&1"


rule prepare_cafe_pilot_inputs:
    input:
        orthogroups="results/orthofinder/exp_run/Results_Jul03/Orthogroups/Orthogroups.tsv",
        manifest="results/tables/of_expansion_manifest.tsv",
        tree="results/orthofinder/exp_run/Results_Jul03/Species_Tree/SpeciesTree_rooted.txt",
        done="results/orthofinder/exp_run.done"
    output:
        counts="results/cafe/pilot/counts.tsv",
        tree="results/cafe/pilot/species_tree.nwk",
        summary="results/tables/cafe_pilot_input_summary.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/cafe_pilot_inputs.log"
    shell:
        "python scripts/prepare_cafe_inputs.py "
        "--orthogroups {input.orthogroups} "
        "--manifest {input.manifest} "
        "--tree {input.tree} "
        "--output {output.counts} "
        "--tree-output {output.tree} "
        "--summary {output.summary} "
        "--max-family-size 100 "
        "--min-species-present 2 > {log} 2>&1"


rule prepare_cafe_ultrametric_tree:
    input:
        calibrations=config["metadata"]["cafe_time_calibrations"],
        manifest="results/tables/of_expansion_manifest.tsv"
    output:
        tree="results/cafe/pilot/species_tree_ultrametric.nwk",
        summary="results/tables/cafe_time_tree_summary.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/cafe_time_tree.log"
    shell:
        "python scripts/prepare_ultrametric_tree.py "
        "--calibrations {input.calibrations} "
        "--manifest {input.manifest} "
        "--output {output.tree} "
        "--summary {output.summary} > {log} 2>&1"


rule prepare_cafe_diagnostic_input:
    input:
        orthogroups="results/orthofinder/exp_run/Results_Jul03/Orthogroups/Orthogroups.tsv",
        tree="results/cafe/pilot/species_tree_ultrametric.nwk",
        panels="config/cafe_diagnostic_panels.tsv"
    output:
        counts="results/cafe/diagnostic/inputs/{panel}.counts.tsv",
        tree="results/cafe/diagnostic/inputs/{panel}.tree.nwk",
        summary="results/cafe/diagnostic/inputs/{panel}.summary.tsv"
    wildcard_constraints:
        panel="conifer_root|conifer_strict|host_root|host_strict"
    conda:
        "../envs/base.yaml"
    benchmark:
        "results/benchmarks/cafe_diagnostic_input_{panel}.txt"
    params:
        python=sys.executable
    log:
        "results/logs/cafe_diagnostic_input_{panel}.log"
    shell:
        "\"{params.python}\" scripts/prepare_cafe_diagnostic_inputs.py "
        "--orthogroups {input.orthogroups} "
        "--tree {input.tree} "
        "--panels {input.panels} "
        "--panel-id {wildcards.panel} "
        "--counts-output {output.counts} "
        "--tree-output {output.tree} "
        "--summary-output {output.summary} > {log} 2>&1"


rule run_cafe5_diag_conifer_root_base:
    input:
        counts="results/cafe/diagnostic/inputs/conifer_root.counts.tsv",
        tree="results/cafe/diagnostic/inputs/conifer_root.tree.nwk"
    output:
        "results/cafe/diagnostic/runs/conifer_root_base_r1/Base_results.txt"
    threads: 4
    benchmark:
        "results/benchmarks/cafe_diag_conifer_root_base_r1.txt"
    log:
        "results/logs/cafe_diag_conifer_root_base_r1.log"
    params:
        outdir="results/cafe/diagnostic/runs/conifer_root_base_r1"
    shell:
        "cafe5 "
        "-i {input.counts} -t {input.tree} -o {params.outdir} -c {threads} > {log} 2>&1"


rule run_cafe5_diag_conifer_root_error:
    input:
        counts="results/cafe/diagnostic/inputs/conifer_root.counts.tsv",
        tree="results/cafe/diagnostic/inputs/conifer_root.tree.nwk"
    output:
        results="results/cafe/diagnostic/runs/conifer_root_error_est_r1/Base_results.txt",
        error="results/cafe/diagnostic/runs/conifer_root_error_est_r1/Base_error_model.txt"
    threads: 4
    benchmark:
        "results/benchmarks/cafe_diag_conifer_root_error_est_r1.txt"
    log:
        "results/logs/cafe_diag_conifer_root_error_est_r1.log"
    params:
        outdir="results/cafe/diagnostic/runs/conifer_root_error_est_r1"
    shell:
        "cafe5 "
        "-i {input.counts} -t {input.tree} -e -o {params.outdir} -c {threads} > {log} 2>&1"


rule run_cafe5_diag_conifer_root_gamma2_error:
    input:
        counts="results/cafe/diagnostic/inputs/conifer_root.counts.tsv",
        tree="results/cafe/diagnostic/inputs/conifer_root.tree.nwk",
        error="results/cafe/diagnostic/runs/conifer_root_error_est_r1/Base_error_model.txt"
    output:
        "results/cafe/diagnostic/runs/conifer_root_gamma2_error_r1/Gamma_results.txt"
    threads: 4
    benchmark:
        "results/benchmarks/cafe_diag_conifer_root_gamma2_error_r1.txt"
    log:
        "results/logs/cafe_diag_conifer_root_gamma2_error_r1.log"
    params:
        outdir="results/cafe/diagnostic/runs/conifer_root_gamma2_error_r1"
    shell:
        "cafe5 "
        "-i {input.counts} -t {input.tree} -e{input.error} -k 2 "
        "-o {params.outdir} -c {threads} > {log} 2>&1"


rule run_cafe5_diag_conifer_strict_error:
    input:
        counts="results/cafe/diagnostic/inputs/conifer_strict.counts.tsv",
        tree="results/cafe/diagnostic/inputs/conifer_strict.tree.nwk"
    output:
        results="results/cafe/diagnostic/runs/conifer_strict_error_est_r1/Base_results.txt",
        error="results/cafe/diagnostic/runs/conifer_strict_error_est_r1/Base_error_model.txt"
    threads: 4
    benchmark:
        "results/benchmarks/cafe_diag_conifer_strict_error_est_r1.txt"
    log:
        "results/logs/cafe_diag_conifer_strict_error_est_r1.log"
    params:
        outdir="results/cafe/diagnostic/runs/conifer_strict_error_est_r1"
    shell:
        "cafe5 "
        "-i {input.counts} -t {input.tree} -e -o {params.outdir} -c {threads} > {log} 2>&1"


rule run_cafe5_diag_host_root_base:
    input:
        counts="results/cafe/diagnostic/inputs/host_root.counts.tsv",
        tree="results/cafe/diagnostic/inputs/host_root.tree.nwk"
    output:
        "results/cafe/diagnostic/runs/host_root_base_r1/Base_results.txt"
    threads: 4
    benchmark:
        "results/benchmarks/cafe_diag_host_root_base_r1.txt"
    log:
        "results/logs/cafe_diag_host_root_base_r1.log"
    params:
        outdir="results/cafe/diagnostic/runs/host_root_base_r1"
    shell:
        "cafe5 "
        "-i {input.counts} -t {input.tree} -o {params.outdir} -c {threads} > {log} 2>&1"


rule run_cafe5_diag_host_root_error:
    input:
        counts="results/cafe/diagnostic/inputs/host_root.counts.tsv",
        tree="results/cafe/diagnostic/inputs/host_root.tree.nwk"
    output:
        results="results/cafe/diagnostic/runs/host_root_error_est_r1/Base_results.txt",
        error="results/cafe/diagnostic/runs/host_root_error_est_r1/Base_error_model.txt"
    threads: 4
    benchmark:
        "results/benchmarks/cafe_diag_host_root_error_est_r1.txt"
    log:
        "results/logs/cafe_diag_host_root_error_est_r1.log"
    params:
        outdir="results/cafe/diagnostic/runs/host_root_error_est_r1"
    shell:
        "cafe5 "
        "-i {input.counts} -t {input.tree} -e -o {params.outdir} -c {threads} > {log} 2>&1"


rule run_cafe5_diag_host_root_gamma2_error:
    input:
        counts="results/cafe/diagnostic/inputs/host_root.counts.tsv",
        tree="results/cafe/diagnostic/inputs/host_root.tree.nwk",
        error="results/cafe/diagnostic/runs/host_root_error_est_r1/Base_error_model.txt"
    output:
        "results/cafe/diagnostic/runs/host_root_gamma2_error_r1/Gamma_results.txt"
    threads: 4
    benchmark:
        "results/benchmarks/cafe_diag_host_root_gamma2_error_r1.txt"
    log:
        "results/logs/cafe_diag_host_root_gamma2_error_r1.log"
    params:
        outdir="results/cafe/diagnostic/runs/host_root_gamma2_error_r1"
    shell:
        "cafe5 "
        "-i {input.counts} -t {input.tree} -e{input.error} -k 2 "
        "-o {params.outdir} -c {threads} > {log} 2>&1"


rule run_cafe5_diag_host_strict_error:
    input:
        counts="results/cafe/diagnostic/inputs/host_strict.counts.tsv",
        tree="results/cafe/diagnostic/inputs/host_strict.tree.nwk"
    output:
        results="results/cafe/diagnostic/runs/host_strict_error_est_r1/Base_results.txt",
        error="results/cafe/diagnostic/runs/host_strict_error_est_r1/Base_error_model.txt"
    threads: 4
    benchmark:
        "results/benchmarks/cafe_diag_host_strict_error_est_r1.txt"
    log:
        "results/logs/cafe_diag_host_strict_error_est_r1.log"
    params:
        outdir="results/cafe/diagnostic/runs/host_strict_error_est_r1"
    shell:
        "cafe5 "
        "-i {input.counts} -t {input.tree} -e -o {params.outdir} -c {threads} > {log} 2>&1"


rule summarize_cafe5_diagnostic_runs:
    input:
        manifest="config/cafe_diagnostic_runs.tsv",
        results=[
            "results/cafe/diagnostic/runs/conifer_root_base_r1/Base_results.txt",
            "results/cafe/diagnostic/runs/conifer_root_error_est_r1/Base_results.txt",
            "results/cafe/diagnostic/runs/conifer_root_gamma2_error_r1/Gamma_results.txt",
            "results/cafe/diagnostic/runs/conifer_strict_error_est_r1/Base_results.txt",
            "results/cafe/diagnostic/runs/host_root_base_r1/Base_results.txt",
            "results/cafe/diagnostic/runs/host_root_error_est_r1/Base_results.txt",
            "results/cafe/diagnostic/runs/host_root_gamma2_error_r1/Gamma_results.txt",
            "results/cafe/diagnostic/runs/host_strict_error_est_r1/Base_results.txt",
        ]
    output:
        models="results/tables/cafe5_diagnostic_model_comparison.tsv",
        convergence="results/tables/cafe5_diagnostic_convergence.tsv",
        report="results/reports/cafe5_diagnostic_model_comparison.md"
    conda:
        "../envs/base.yaml"
    benchmark:
        "results/benchmarks/cafe5_diagnostic_summary.txt"
    log:
        "results/logs/cafe5_diagnostic_summary.log"
    shell:
        "python scripts/summarize_cafe_diagnostic_runs.py --manifest {input.manifest} "
        "--output {output.models} --convergence {output.convergence} --report {output.report} "
        "> {log} 2>&1"


rule build_cafe5_diagnostic_candidate_stability:
    input:
        models="results/tables/cafe5_diagnostic_model_comparison.tsv",
        convergence="results/tables/cafe5_diagnostic_convergence.tsv",
        manifest="config/cafe_diagnostic_runs.tsv",
        candidates="results/tables/mechanism_axis_summary.tsv",
        decisions="results/tables/model_supported_candidate_review_decisions.tsv"
    output:
        stability="results/tables/cafe5_diagnostic_candidate_stability.tsv",
        report="results/reports/cafe5_strengthening_diagnostic.md"
    conda:
        "../envs/base.yaml"
    benchmark:
        "results/benchmarks/cafe5_diagnostic_candidate_stability.txt"
    log:
        "results/logs/cafe5_diagnostic_candidate_stability.log"
    shell:
        "python scripts/build_cafe_diagnostic_candidate_stability.py "
        "--model-table {input.models} --convergence-table {input.convergence} "
        "--run-manifest {input.manifest} --candidates {input.candidates} "
        "--review-decisions {input.decisions} --output {output.stability} "
        "--report {output.report} > {log} 2>&1"


rule run_cafe5_full_pilot:
    input:
        counts="results/cafe/pilot/counts.tsv",
        tree="results/cafe/pilot/species_tree_ultrametric.nwk"
    output:
        directory("results/cafe/full/cafe5_ultrametric_base")
    log:
        "results/logs/cafe5_full_pilot.log"
    threads:
        4
    shell:
        "cafe5 "
        "-i {input.counts} "
        "-t {input.tree} "
        "-o {output} "
        "-c {threads} > {log} 2>&1"


rule parse_cafe5_full_pilot:
    input:
        cafe_dir="results/cafe/full/cafe5_ultrametric_base",
        top_candidates="results/tables/regional_family_validation_summary.tsv"
    output:
        model="results/tables/cafe5_full_model_summary.tsv",
        clades="results/tables/cafe5_full_clade_results.tsv",
        families="results/tables/cafe5_full_family_results.tsv",
        branches="results/tables/cafe5_full_branch_results.tsv",
        top_intersections="results/tables/cafe5_top_candidate_intersections.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/cafe5_full_parse.log"
    shell:
        "python scripts/parse_cafe5_results.py "
        "--cafe-dir {input.cafe_dir} "
        "--top-candidates {input.top_candidates} "
        "--model-summary-output {output.model} "
        "--clade-output {output.clades} "
        "--family-output {output.families} "
        "--branch-output {output.branches} "
        "--top-intersection-output {output.top_intersections} > {log} 2>&1"
