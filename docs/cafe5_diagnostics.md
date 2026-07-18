# CAFE5 pilot notes

This project uses CAFE5 only after OrthoFinder-derived gene-family counts have
been converted to CAFE format:

```text
Desc    Family ID    species_1    species_2    ...
```

The pilot input files are:

- `results/cafe/pilot/counts.tsv`
- `results/cafe/pilot/species_tree.nwk`
- `results/cafe/pilot/species_tree_ultrametric.nwk`
- `results/tables/cafe_pilot_input_summary.tsv`
- `results/tables/cafe_time_tree_summary.tsv`

`species_tree.nwk` is the OrthoFinder rooted species tree. It is useful for
topology inspection, but it should not be used as the formal CAFE tree.
`species_tree_ultrametric.nwk` is generated from
`config/cafe_time_calibrations.tsv` and is the CAFE pilot time tree.

On Windows, write CAFE input TSV files with Unix LF line endings. CRLF line
endings can cause CAFE5 under WSL/Linux to read the final species column as a
different name.

## WSL 24.04 smoke test

Install CAFE5 in a Linux micromamba environment:

```bash
micromamba create -y -n cafe5 -c conda-forge -c bioconda cafe=5.1.0
```

Run a small format and environment check:

```bash
cd /path/to/pine-pwn-defense-genomics
mkdir -p results/cafe/smoke
head -n 501 results/cafe/pilot/counts.tsv > results/cafe/smoke/counts_500.tsv
cafe5 \
  -i results/cafe/smoke/counts_500.tsv \
  -t results/cafe/pilot/species_tree_ultrametric.nwk \
  -o results/cafe/smoke/cafe5_ultrametric_base \
  -c 4
```

The smoke test is only a technical check. It should not be interpreted as a
formal expansion/contraction result because it uses a 500-family subset and a
provisional pilot time tree.

## Full pilot run

After the OrthoFinder expansion run was updated to include `plon`, the pilot
CAFE input contained 34,002 filtered families across 10 species. A full pilot
CAFE5 Base-model run was executed under WSL Ubuntu 24.04:

```bash
cd /path/to/pine-pwn-defense-genomics
cafe5 \
  -i results/cafe/pilot/counts.tsv \
  -t results/cafe/pilot/species_tree_ultrametric.nwk \
  -o results/cafe/full/cafe5_ultrametric_base \
  -c 4
```

Parsed outputs:

- `results/tables/cafe5_full_model_summary.tsv`
- `results/tables/cafe5_full_clade_results.tsv`
- `results/tables/cafe5_full_family_results.tsv`
- `results/tables/cafe5_full_branch_results.tsv`
- `results/tables/cafe5_top_candidate_intersections.tsv`

The Base model estimated lambda at `0.0050403229920973`, essentially at the
maximum possible lambda for the topology (`0.00504032`). Therefore, this full
run is a model-based pilot prioritization layer rather than a final formal
gene-family expansion/contraction analysis.

## Candidate evidence ranking

The CAFE5 pilot output is intersected with the regional top orthogroups and
the MAFFT/FastTree family-review flags in:

- `results/tables/candidate_evidence_ranking.tsv`
- `results/tables/strong_candidate_orthogroups.tsv`

The ranking is intentionally conservative:

- `strong_candidate`: CAFE5 family-level significant, at least one significant
  branch-level change, and `priority_low_risk` family-review status.
- `moderate_candidate`: CAFE5 family-level and branch-level support, but the
  family requires outgroup-context interpretation.
- `caution_candidate`: CAFE5 support is present, but the family is flagged as
  single-species-driven or a large single-species-driven family.
- `weak_or_unresolved`: CAFE5 support is absent or the family was not tested in
  the current pilot run.

In the current pilot, no orthogroup qualifies as `strong_candidate`; this is an
important result rather than a failed output. It means the strongest
model-supported candidates still require manual inspection before biological
interpretation.

The model-supported subset is exported for manual review in:

- `results/tables/model_supported_candidate_review.tsv`
- `results/tables/model_supported_candidate_manual_checklist.tsv`
- `results/tables/model_supported_candidate_tree_metrics.tsv`
- `results/reports/model_supported_candidate_review.md`
- `results/reports/model_supported_candidate_tree_review.md`
- `results/reports/model_supported_candidate_tree_review_ascii/`
- `results/tables/model_supported_candidate_review_decisions.tsv`
- `results/reports/model_supported_candidate_review_decisions.md`
- `results/tables/manuscript_candidate_summary.tsv`
- `results/reports/manuscript_candidate_summary.md`
- `results/tables/mechanism_axis_summary.tsv`
- `results/tables/mechanism_evidence_gaps.tsv`
- `results/reports/mechanism_strategy.md`
- `results/tables/public_expression_evidence_plan.tsv`
- `results/reports/public_expression_evidence_plan.md`
- `results/tables/public_expression_text_download_manifest.tsv`
- `results/tables/public_expression_literature_evidence.tsv`
- `results/reports/public_expression_literature_evidence.md`
- `results/tables/public_expression_supplement_download_manifest.tsv`
- `results/tables/public_expression_supplement_inventory.tsv`
- `results/reports/public_expression_supplement_inventory.md`
- `results/tables/effector_target_network_expression_support.tsv`
- `results/reports/effector_target_network_expression_support.md`
- `results/tables/effector_target_network_supplement_gene_hits.tsv`
- `results/reports/effector_target_network_supplement_gene_hits.md`

The current review packet contains 10 orthogroups: five `moderate_candidate`
families that should be inspected for outgroup context and branch support, and
five `caution_candidate` families that should be checked for single-species or
large-family drivers. The checklist keeps the interpretation conservative by
recording the nematode role as an indirect effector-guided anchor, not a CAFE5
covariate, and by requiring manual tree review before any model-supported
candidate can be promoted. The tree-metric prescreen provides topology
triage only: in the current output, no model-supported candidate forms a
strict focal-region clade, four show focal-region-biased tree-wide
distributions, and six are region-intermingled. Biopython is used to parse the
candidate Newick trees and export ASCII topology snapshots for manual review.
The review-decision table then assigns conservative use labels. In the current
output, two candidates are priority manual-review cases without regional-claim
support, one is retained as a functional-context candidate with mixed topology,
two remain unresolved supplementary-review candidates, and five caution-tier
families are deprioritized to background or supplement-only use. The manuscript
candidate summary converts those decisions into claim ceilings for writing:
`OG0003130` and `OG0003639` can be discussed as priority manual-review
candidates, `OG0002098` can be used as functional context, and the remaining
families should remain supplementary or background material. This packet is the
hand-off from automated ranking to manual family-tree review.

## Mechanism-axis reframing

The current higher-impact development path is tracked in:

- `results/tables/mechanism_axis_summary.tsv`
- `results/tables/mechanism_evidence_gaps.tsv`
- `results/reports/mechanism_strategy.md`

This layer reframes the project around three mechanism-generating hypotheses:

- Route 1, effector-target network vulnerability: PWN effectors may converge on
  host immune-signaling and ROS network nodes rather than on clean regional
  gene-family expansions. Current priority candidates are `OG0003130` and
  `OG0003639`, with `OG0002098` and unresolved immune-signaling families as
  supporting context.
- Route 2, hydraulic/xylem collapse: PWD outcome may reflect coupled
  cell-wall, xylem, wound, lignin, and ROS module disruption. Current evidence
  is supporting only and requires infection-response data before it can become
  a mechanism claim.
- Route 3, susceptibility vulnerability: the negative tree and copy-number
  evidence is used constructively to test whether susceptibility reflects
  defense-module vulnerability or network fragility rather than simple regional
  expansion.

All three routes currently require infection-expression support and target or
network support before being written as mechanisms. The public-expression plan
maps the infection-expression gaps to curated public PWD expression or
resistant/susceptible sources in
`config/public_expression_evidence_sources.tsv`. Rows for target or network
support are retained as explicit gap markers, because expression evidence alone
does not validate direct effector targets or network centrality.

The first public-expression extraction layer caches source text from public
article pages and searches the pre-reference body text for candidate-family
keywords. This produced a reproducible literature evidence index, not a
quantitative expression matrix. In the current run, four of six public pages
were cached successfully, two publisher pages returned HTTP 403, and 10
candidate-source rows contained pre-reference keyword hits. These hits support
manual follow-up for WRKY, ERF, and MYB immune-signaling candidates, but they
do not yet identify the exact orthogroup member that is differentially
expressed.

Route 1 is now summarized separately in
`results/tables/effector_target_network_expression_support.tsv`. The current route 1 packet
contains five candidate orthogroups. Four have family-level literature support
from infection-expression article text: the WRKY candidate `OG0003130`, ERF
candidates `OG0003639` and `OG0000888`, and the MYB candidate `OG0004410`.
The superoxide dismutase candidate `OG0002098` does not yet have a keyword hit
in the cached public-expression text. PMC supplementary xlsx links for the
P. densiflora and P. massoniana studies were recorded in
`config/public_expression_supplement_sources.tsv`. Direct PMC article-bin
downloads return HTML challenge pages, but the same files were recovered from
PMC OA package archives under the 2026 `deprecated/oa_package` FTP layout.
Five xlsx supplements are now cached and inventoried. A workbook-level scan
found no exact route 1 gene-ID matches and no route 1 keyword hits inside these
supplements. The main blocker for exact expression support is identifier-space
mismatch: the P. massoniana supplements use `PITA_...` expression IDs, while
the current P. massoniana proteome candidates use `gmmutg...` or `STRG...`
IDs; the P. densiflora supplements use Trinity and Arabidopsis/PLAZA mappings
rather than the `Pd...` genome annotation IDs used in the orthology panel.
The dedicated crosswalk audit in
`results/tables/effector_target_network_identifier_crosswalk_audit.tsv` records this as zero
exact candidate-gene matches for all route 1 supplement sources: pden requires
a Trinity-to-Pd bridge, pmas requires a PITA-to-gmmutg/STRG bridge, and ptab
currently has no route 1 supplement source. Until one of these bridges is
added, the expression layer should be described as family-level literature
context rather than exact candidate-gene DEG evidence.

A first indirect functional bridge audit is now exported in
`results/tables/effector_target_network_functional_bridge_audit.tsv` and
`results/reports/effector_target_network_functional_bridge_audit.md`. The audit adds a
UniProt-derived Arabidopsis seed locus map in
`config/effector_target_network_arabidopsis_locus_map.tsv` and scans downloaded expression
supplements for Trinity-to-Arabidopsis-locus mappings. This recovered eight
source-orthogroup indirect bridges across `pden_dataset2` and `pden_dataset4`,
covering all four route 1 orthogroups represented in the network seed table:
`OG0002098` (FSD1/FSD3; 3 Trinity IDs), `OG0003130` (WRKY4; 2 Trinity IDs),
`OG0003639` (RAV1/RAV2; 4 Trinity IDs), and `OG0004410` (MYB3/MYB5; 33
Trinity IDs). The two P. densiflora supplements support the same transcript
sets through BLAST-style annotation and curated Trinity-to-Arabidopsis mapping,
respectively. These are useful traceable leads for expression follow-up, but
they are not exact current Pd gene-ID matches because the bridge still stops at
Trinity transcript to Arabidopsis homolog annotation. The next necessary step
is to recover the DEG/statistics table for these Trinity IDs and, if possible,
align or map those transcripts back to the current *P. densiflora* protein
identifiers.

The corresponding Trinity sequences have now been extracted from the nested
`Trinity.fasta` inside `PMC6704138/41598_2019_48660_MOESM2_ESM.zip`. Outputs
are `data/interim/effector_target_network_pden_trinity_bridge_transcripts.fasta`,
`results/tables/effector_target_network_pden_trinity_sequence_manifest.tsv`, and
`results/reports/effector_target_network_pden_trinity_sequence_extraction.md`. All 42 target
Trinity IDs were found, yielding 126 isoform sequence records: 8 for
`OG0002098`, 5 for `OG0003130`, 11 for `OG0003639`, and 102 for `OG0004410`
when isoform records are counted. This prepares the sequence-level input for
Trinity-to-current-Pd mapping. It still does not constitute a current Pd gene
ID match or DEG call.

A DIAMOND blastx mapping from these extracted Trinity isoforms to the current
*P. densiflora* HA protein set is exported in
`results/tables/effector_target_network_pden_trinity_vs_current_pd_blastx.tsv`,
`results/tables/effector_target_network_pden_trinity_to_pd_mapping.tsv`, and
`results/reports/effector_target_network_pden_trinity_to_pd_mapping.md`. With filtering at
bitscore >= 80 and percent identity >= 35, all 42 Trinity targets had at least
one current Pd protein hit. Eight Trinity targets had a sequence-level hit to a
current Pd gene already present in the same route 1 candidate orthogroup:
3 for `OG0002098`, 1 for `OG0003130`, 3 for `OG0003639`, and 1 for
`OG0004410`. The remaining 34 targets hit current Pd proteins outside the
route 1 candidate set, often likely reflecting broader family-level homology.
This is the first direct sequence-level bridge from the P. densiflora
Trinity expression supplements back to the current Pd proteome, but it remains
separate from DEG/statistical expression evidence.

A dedicated public-source audit is now exported in
`results/tables/effector_target_network_pden_expression_source_audit.tsv`,
`results/tables/effector_target_network_pden_candidate_expression_audit.tsv`, and
`results/reports/effector_target_network_pden_expression_source_audit.md`. It inspects the five
available PMC supplementary materials for the P. densiflora study. The audit
classifies Dataset 1 as DEG figures and TF-family DETF summary context, Dataset
2 as the Trinity FASTA resource, Dataset 3 as a BLAST/sequence annotation
table, Dataset 4 as qRT-PCR validation code or primer context, and Dataset 5 as
a Trinity-to-Arabidopsis/PLAZA identifier crosswalk. None of these public
materials contains a transcript-level DEG statistics table with Trinity IDs,
log-fold changes, and FDR/q-value columns. Therefore the eight current-Pd route
1 candidate bridges are explicitly marked as
`sequence_bridge_only_public_deg_stats_unavailable`, with raw RNA-seq
reanalysis or recovery of the authors' DEG table as the required next step.

The raw-data reanalysis entry point is now captured in
`results/tables/effector_target_network_pden_rnaseq_run_manifest.tsv`,
`results/tables/effector_target_network_pden_rnaseq_metadata_audit.tsv`, and
`results/reports/effector_target_network_pden_rnaseq_metadata.md`. The article-linked study is
`SRP165817` / `PRJNA496563`, not `PRJDB19336`; the latter is retained only as
an exclusion audit because it represents a different mixed-species project.
ENA metadata resolves nine paired-end Illumina runs into three biological
replicates each for pathogenic `B. xylophilus`, non-pathogenic
`B. thailandae`, and water control. The primary planned contrast is
`B. xylophilus` versus `B. thailandae`, with both nematode-versus-water
contrasts retained for interpretation. FASTQ download and quantification are
kept as separate stages so the full transfer volume and MD5 checks can be
reviewed before raw reads are fetched. Quantification will use the published
Trinity reference and will not require a new transcriptome assembly.

The Galaxy API pilot was completed on UseGalaxy.eu 26.1 in the history
`PWN effector_target_network - P. densiflora RNA-seq - SRP165817`. The published Trinity
archive was exported reproducibly as
`data/interim/pden_published_trinity_reference.fasta` (185,501 transcript
records), and `data/interim/pden_published_trinity_tx2gene.tsv` maps those
records to the 72,864 Trinity genes reported by the study. Pilot run
`SRR8061568` was fetched directly from ENA with MD5 validation, represented as
a paired collection, checked with FastQC 0.74 and MultiQC 1.35, and quantified
with Salmon 1.10.1 using automatic library-type inference plus sequence- and
GC-bias correction. No BAM output was generated. The transcript and gene
quantification files were downloaded to
`data/external/galaxy/effector_target_network_pden_pilot/` and summarized in
`results/tables/effector_target_network_pden_salmon_pilot_candidate_abundance.tsv` and
`results/reports/effector_target_network_pden_salmon_pilot_candidate_abundance.md`. Thirty-eight
of 42 candidate Trinity genes, including seven of eight sequence-linked
current-Pd route 1 bridges, had positive reads in this water-control pilot.
This result validates identifier continuity and quantifiability only; it is
not infection differential-expression evidence. The remaining 16 FASTQ files
were submitted to a nine-sample `list:paired` collection after the pilot
passed, with full contrasts deferred until all fetch and quantification jobs
are complete.

The route 1 mechanism track now has a network-seed layer in
`results/tables/effector_target_network_effector_target_network_seed.tsv` and
`results/reports/effector_target_network_effector_target_network_seed.md`. This is the main
route 1 hand-off from copy-number prioritization into predicted
effector-target network analysis. It connects PWN effector classes to host
candidate modules and Arabidopsis/SwissProt homologs: protease and mimicry
classes connect to immune-signaling WRKY/ERF seeds (`OG0003130`, `OG0003639`),
and detoxification classes connect to the ROS-detoxification SOD seed
(`OG0002098`). A MYB candidate (`OG0004410`) is retained as secondary network
context. This layer should be used as input for plant interolog/PPI mapping and
centrality analysis, while retaining the claim boundary that these are
predicted network vulnerability seeds, not validated effector targets.

The first STRING overlay is implemented in
`results/tables/effector_target_network_string_seed_mapping.tsv`,
`results/tables/effector_target_network_string_interactions.tsv`, and
`results/tables/effector_target_network_string_network_centrality.tsv`. Using the
version-specific STRING 12.0 API for Arabidopsis (`species=3702`) with a
functional-network confidence threshold of 700 and up to 20 interaction
partners per seed, all eight Arabidopsis seed symbols mapped successfully.
Seven had STRING neighbors in the fetched local neighborhood subnetwork. The
SOD seeds `FSD1` and `FSD3` had 20 partners each, the ERF/RAV-family seeds
`RAV1`, `RAV2`, and `TEM1` had 8, 3, and 15 partners, and the MYB context
seeds `MYB3` and `MYB5` had 6 and 16 partners. `WRKY4` mapped successfully but
had no neighbor at the current threshold. These values are local neighborhood
metrics for hypothesis prioritization, not genome-wide Arabidopsis centrality
estimates.

The STRING evidence-channel decomposition is exported in
`results/tables/effector_target_network_string_evidence_channels.tsv` and
`results/reports/effector_target_network_string_evidence_channels.md`. Using a channel threshold
of 0.4, the SOD seeds provide the strongest route 1 network evidence:
`FSD1` has 11 experimental-channel and 9 database-channel partners, and `FSD3`
has 9 experimental-channel and 7 database-channel partners after deduplicating
reciprocal STRING edges. In contrast, `RAV1`, `RAV2`, `TEM1`, `MYB3`, and
`MYB5` are currently supported only by coexpression/text-mining channel
evidence at this threshold, and `WRKY4` remains mapped but without neighbors
under the current STRING settings. This refines the route 1 claim: the ROS/SOD
axis is the strongest PPI/interolog-supported subroute, whereas WRKY/ERF/MYB
remain predicted regulatory-network hypotheses requiring additional evidence
from a different PPI source, relaxed threshold sensitivity analysis, or
infection-expression overlay.

The ROS/SOD subroute is further decomposed in
`results/tables/effector_target_network_ros_neighbor_modules.tsv`,
`results/tables/effector_target_network_ros_neighbor_module_summary.tsv`, and
`results/reports/effector_target_network_ros_neighbor_modules.md`. This layer classifies the
FSD1/FSD3 STRING local-neighborhood edges into ROS-relevant modules. The
current table contains 40 FSD neighbor edges, including 28
experimental/database-supported edges. The clearest module is the SOD family
itself, with 14 edges and 14 experimental/database-supported associations
linking FSD1/FSD3 to CSD, FSD, and MSD family members. Additional ROS-relevant
modules include catalase ROS scavenging (CAT1-2, CAT2), ascorbate peroxidase
ROS scavenging (APX1, APX2), chloroplast redox/gene-expression context
(MRL7/PTAC/FLN neighbors), thioredoxin-related redox context (CITRX), and
copper homeostasis (SPL7). Several experimental/database-supported neighbors
remain symbolically unclassified and should be manually annotated before being
used in a mechanism claim. This supports a more specific route 1 subhypothesis:
the effector-guided SOD candidate may mark vulnerability in chloroplast-linked
ROS detoxification and redox-buffering network neighborhoods.

The first infection-expression overlay for this ROS subroute is exported in
`results/tables/effector_target_network_ros_expression_overlay.tsv` and
`results/reports/effector_target_network_ros_expression_overlay.md`. This layer scans cached
public expression article text for the FSD1/FSD3 neighbor-module keywords and
neighbor symbols. The result is intentionally conservative: among 42
module-source rows, only one text-context hit was found, in the
`pden_pathogenic_vs_nonpathogenic` source, where the article describes
reactive oxygen species (ROS) production upon nematode infection and its
relationship to defense and hypersensitive response. No cached source currently
provides module-specific text hits for FSD/CSD/MSD, CAT, APX, PTAC/FLN, CITRX,
or SPL7 terms before the reference section. Therefore, the route 1 ROS/SOD
subroute currently has stronger PPI/interolog support than infection-expression
support. The expression layer should be described as broad ROS infection
context until exact transcript/gene-level expression data can be bridged.

The workbook-level supplement overlay is exported in
`results/tables/effector_target_network_ros_supplement_overlay.tsv` and
`results/reports/effector_target_network_ros_supplement_overlay.md`. This scan adds GO/KEGG/WGCNA
supplement context for the same ROS neighbor modules and now records the best
available adjusted P value for workbook enrichment rows. Among 35 module-source
rows, three workbook keyword hits were detected in the *P. massoniana* GO
supplement: broad antioxidant/oxidoreductase context, copper-homeostasis
context, and the corresponding unclassified ROS-neighbor context. None of these
hits was significant after adjustment (`adjusted_p_lt_0.05` count = 0), and no
precise SOD/CAT/APX module hit was recovered from the supplements. This makes
the current expression conclusion sharper rather than stronger: supplement
workbooks provide weak broad ROS/redox context, not statistically supported
module enrichment or exact route 1 pine gene expression support.

## Strengthening diagnostic, 2026-07-13

The predeclared strengthening diagnostic is recorded in
`config/cafe_diagnostic_panels.tsv`, `config/cafe_diagnostic_runs.tsv`, and
`results/reports/cafe5_strengthening_diagnostic.md`. Eight first-pass models
were run across root-compatible and strict conifer-wide and Pinus host-focused
panels. All six Base or error-aware Base models remained at or near the
topology-specific maximum lambda. The two error-aware Gamma k=2 models moved
lambda into the interior, but CAFE5 reported high family-level calculation
failure rates for 13,177 of 14,212 conifer-wide families and all 15,540
host-focused families.

No first-pass model therefore met the combined acceptance criteria of an
interior lambda, clean numerical quality, and replicated convergence. The
predeclared stop rule was applied before additional replicates. CAFE5 remains
a supplementary sensitivity analysis under the provisional time tree and
heterogeneous public proteome annotations. Its lineage-level results must not
be interpreted as an East Asia versus North America resistance test, direct
host-parasite coevolution, or effector-target validation.

## Formal-analysis requirements

Before reporting CAFE5 results scientifically:

- Replace the OrthoFinder rooted tree with a justified time-calibrated
  ultrametric species tree, or replace the provisional values in
  `config/cafe_time_calibrations.tsv` with curated literature values.
- Run the full filtered family set, not the smoke subset.
- Inspect lambda convergence, maximum possible lambda warnings, family-size
  filters, and branch-level calls.
- Treat significant families as candidate lineage-specific expansions or
  contractions, not direct evidence of host-nematode coevolution.

Primary documentation:

- https://github.com/hahnlab/CAFE5
- https://gensoft.pasteur.fr/docs/CAFE5/5.1.0/
