# Reproduction guide

1. Run the reduced FASTA example and test suite to verify the Python environment.
2. Download public source data using the accessions and provenance tables in `config/`.
3. Execute the workflow in stages: metadata, proteomes, orthology, annotations, secretome, host mapping, gene-family diagnostics, and integrated figures.
4. Compare key derived tables with the companion archive using `MANIFEST.tsv` checksums.
5. Rebuild publication figures from frozen figure-data tables rather than editing upstream results during visualization.

Software environments are separated by task in `workflow/envs/`. CAFE5 is expected to be available as `cafe5` inside its Linux environment. Resource-intensive orthology, RNA-seq, and structure steps should be scheduled independently from the reduced validation suite.

