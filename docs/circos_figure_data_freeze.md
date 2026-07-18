# Circos figure-data freeze

The manuscript Circos figure uses freeze version `circos_v1`. Its plotting inputs are locked by `results/tables/manuscript_circos_figure_data_freeze_v1.tsv`, which records one SHA-256 checksum and schema checksum per input table.

Frozen components:

- `candidate_tracks`: 12 orthogroups and A-M candidate evidence fields.
- `effector_context`: five linked effector classes and class-level evidence summaries.
- `edges`: functional-prior and candidate-membership links.
- `sectors`: stable effector, module, and candidate sector definitions.

The Snakemake plotting rule verifies all four entries before invoking R. A mismatch stops rendering and names the changed component.

To revise one component, review and regenerate only its source table, run the focused tests for that component, recreate the manifest, and confirm that `git diff` changes only the intended manifest row. A new manuscript data release should use a new freeze version rather than silently replacing `circos_v1`.

Create or verify the manifest with `scripts/freeze_circos_figure_inputs.py --mode create|verify` and the four explicit table arguments recorded in `workflow/rules/08_figures.smk`.
