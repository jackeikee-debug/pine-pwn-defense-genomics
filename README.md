# Pine wilt disease defense genomics

This repository contains the reproducible analysis code for an effector-guided comparative genomics study of pine defense modules associated with pine wilt disease susceptibility and tolerance. It integrates pine proteomes, orthogroups, gene-family diagnostics, public infection transcriptomes, nematode secretome evidence, functional-association networks, and selected structure predictions.

## Scope and interpretation

The workflow prioritizes candidate mechanisms rather than experimentally validated effector-target interactions. Public infection-expression data support host response patterns but do not establish direct targeting or causality. CAFE5 outputs are retained as diagnostic and candidate-ranking evidence; family-level expansion or contraction claims require convergence, time-tree, and annotation review.

## Workflow

The top-level `Snakefile` organizes metadata validation, proteome preparation, orthology, annotation, defense-module screening, nematode secretome analysis, host-module integration, gene-family diagnostics, expression integration, and figures. Conda environments are defined in `workflow/envs/`, parameters in `config/config.yaml`, and analysis steps in `scripts/`.

Scientific evidence axes are named by mechanism:

- `effector_target_network`: predicted concentration of candidate effectors on host defense-network hubs.
- `apoplastic_cell_wall`: extracellular, cell-wall, and oxidative defense candidates.
- `cross_species_expression`: expression support mapped across independent pine infection datasets.

## Reproduction

Create the base environment and inspect the workflow:

```bash
conda env create -f workflow/envs/base.yaml
conda activate pine-pwn-base
snakemake --list
```

Run the reduced example without downloading public data:

```bash
python scripts/summarize_proteomes.py \
  --input-dir examples/smoke/proteins \
  --output examples/smoke/results/proteome_summary.tsv
python -m pytest -q
```

Full reproduction requires the public accessions listed in `config/` and the frozen derived-data archive described in `config/release_data_crosswalk.tsv`. See `docs/reproduction_guide.md` for staged execution.

## Data availability

Raw sequencing reads, reference proteomes, and database resources remain at their original public repositories under the accessions recorded in `config/`. The manuscript-focused Zenodo archive contains only frozen derived tables, compact sequence sets, metadata, and checksums required to support and reproduce reported results. Version 0.1.1 is available at [https://doi.org/10.5281/zenodo.21495920](https://doi.org/10.5281/zenodo.21495920). Version 0.1.1 corrects submission metadata in the reviewer-readable supplementary workbook; scientific data and analytical results are unchanged from version 0.1.0. Large public source files and redundant intermediates are intentionally excluded.

## Repository structure

- `config/`: species, accession, phenotype, mechanism, and parameter tables.
- `workflow/`: Snakemake rules and software environments.
- `scripts/`: command-line analysis and visualization programs.
- `tests/`: synthetic and reduced-data tests.
- `examples/smoke/`: small FASTA example for installation checks.
- `docs/`: workflow, provenance, and reproduction notes.

## Citation

Use the metadata in `CITATION.cff`. The associated manuscript targets *Tree Genetics & Genomes*. Cite the frozen derived-data archive using DOI `10.5281/zenodo.21495920`.

## License

Code is released under the MIT License. Third-party datasets retain their original licenses and terms of use.
