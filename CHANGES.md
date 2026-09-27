# Consolidation into DEJUPipeline 0.1.0

The original five pipelines and `deju_helpers.R` were reviewed together. Version
3 contributed explicit file/sample/design checks; version 5 contributed handling
of the observed Rsubread junction-column schema. Later version numbers alone do
not establish correctness or identify a past successful run.

* Replaced hard-coded WT/human/paired/reverse-strand settings with a validated
  sample sheet, explicit comparison and sequencing settings, and categorical
  batch/subject covariates.
* Preserved the joint exon/junction usage model, with explicit modern and legacy
  edgeR engines, robust fitting, TMM and design-aware filtering.
* Fixed sample-column alignment by identifier rather than an order-dependent
  intersection of column names.
* Removed genomic insert-size filtering, which can reject ordinary spliced pairs.
* Applied primary/unique/mapping-quality/duplicate policies before both channels.
  Both channels now count aligned reads, with the paired-end policy documented.
* Replaced junction collection with chunked CIGAR counting after tests exposed
  Rsubread 2.10.5 nonSplitOnly and strand-specific junction-counting problems.
* Removed arbitrary selection from multigene assignments. Exact annotated
  junctions and uniquely contained novel junctions have explicit assignment
  rules; excluded junctions remain visible in an audit table.
* Derived both flattened exons and transcript junctions from one supplied GTF,
  preserving full IDs and chromosome names and rejecting conflicting loci.
* Fixed minus-strand junction matching by separating genomic coordinate order
  from biological donor/acceptor direction. Retained all compatible transcript
  positions instead of silently choosing the transcript with most introns.
* Converted natural-log usage coefficients to explicitly labeled log2 effects.
* Separated junction, exon, all-feature and gene-level testing families and
  distinguished FDR significance from a descriptive effect-size filter.
* Used eligible tested genes as the GO background with explicit ID/organism
  mapping, keeping gene-test and junction-derived gene selections separate.
* Added reproducible small fixtures, count/model tests, native BAM tests, QC plots,
  saved models, session metadata, input checksums and non-overwriting outputs.

The empty `DJU_pipeline.R` cannot be a runnable candidate. The large saved
workspace did not provide enough provenance to identify the previously used
source version; no such claim is made for this release.
