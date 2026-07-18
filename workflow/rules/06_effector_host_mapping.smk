rule map_effectors_to_host_modules:
    input:
        effectors="results/tables/nematode_candidate_effectors.tsv",
        mapping=config["metadata"]["effector_host_module_map"],
        host_shortlist="results/tables/defense_module_candidate_shortlist.tsv"
    output:
        "results/tables/effector_host_module_links.tsv"
    conda:
        "../envs/base.yaml"
    log:
        "results/logs/effector_host_module_links.log"
    shell:
        "python scripts/map_effectors_to_host_modules.py "
        "--effectors {input.effectors} "
        "--mapping {input.mapping} "
        "--host-shortlist {input.host_shortlist} "
        "--output {output} "
        "--min-confidence medium > {log} 2>&1"
