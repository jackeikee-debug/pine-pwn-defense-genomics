# Data provenance

Public accessions and source URLs are recorded in `config/data_sources.tsv`, `config/species_metadata.tsv`, and the expression and secretome source tables. Analysis parameters are versioned with the workflow.

The companion derived-data archive is allowlist-based. Its `MANIFEST.tsv` records SHA-256 checksums and byte sizes, while `DATA_CATALOG.tsv` explains the evidence layer and manuscript role of every payload. `config/release_data_crosswalk.tsv` connects archive paths to workflow locations and producers.

Raw reads, complete public proteomes, database mirrors, failed runs, browser sessions, and redundant intermediate files are excluded from the archive.

