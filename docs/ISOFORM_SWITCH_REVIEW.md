# Review of the original isoform-switch scripts

Reviewed `iso_helpers.R`, `iso_pipeline.R` (188 lines) and `iso_pipeline_New.R`
(246 lines) from the user's Rscripts directory. Original files are unchanged.
This is a code/API review, not confirmation of a previously successful run.
IsoformSwitchAnalyzeR is not installed in the current R environment.

## What the scripts do

Both consume Salmon quant.sf files and a matching transcript GTF. They import
transcript estimated counts and TPM through tximport, create an
IsoformSwitchAnalyzeR object, filter low-expression/single-isoform genes, test
relative isoform usage using DEXSeq, annotate splicing and predicted ORF/NMD
consequences, and export tables and GO enrichment.

Isoform fraction (IF) is a transcript's share of its gene's abundance. dIF is its
change between conditions. This is distinct from gene-expression log2FC, DEJU
relative junction log2FC, and event-specific PSI. A significant isoform-usage
change does not necessarily mean the dominant isoform changed. The scripts do
not reconstruct transcripts or run Salmon; they start with its quantifications.

## Which version to retain as a starting point

`iso_pipeline_New.R` is preferable structurally: explicit gene/isoform/IF
filter parameters, explicit ID-import settings, and no misleading filtered-GTF
cache. Neither version should be copied unchanged into the new workflow.

The unusual package argument `isoformExonAnnoation` is the documented spelling.
The original `isoformExonAnno` can work through R partial argument matching;
it should be replaced by the complete name, not called a proven fatal typo.
The same caution applies to shortened extractNT/extractAA arguments.

## Confirmed code issues and method choices to correct

1. The original creates `filtered_gtf` but passes `gtf` to importRdata. Its claimed
   filtered-GTF cache is therefore unused. The helper retains transcript/exon
   records but drops CDS records; using that file for coding-consequence analysis
   would lose annotated coding regions.
2. `build_meta()` binds a condition vector to the alphabetical list.files order.
   Equal lengths do not prove correct biological sample labels. Use an explicit
   sample-to-condition-to-quantification sheet; reject duplicate/missing entries.
3. Both designs select only sampleID and condition, discarding batch, subject or
   other covariates present in metadata. ImportRdata can accept additional
   cofactors; comparisons and direction must be explicit.
4. Both call isoformSwitchTestDEXSeq with reduceToSwitchingGenes=TRUE before
   exporting isoform_features. That file is not the complete tested universe.
   Save the full tested object first, then create a separate subset for the
   expensive consequence analysis.
5. alpha and dIF_cutoff are not passed into the DEXSeq call or all downstream
   functions. User settings can therefore disagree with package defaults used
   to discard genes early. Propagate thresholds consistently.
6. GO enrichment omits the tested-gene universe. The background should be the
   mapped genes eligible for the switch test, not the whole organism database.
   Human-only ID guessing and multiVals='first' should be replaced by explicit
   mapping with ambiguous/unmapped-ID reports.
7. The gene summary groups only by gene_id. With multiple comparisons this
   merges distinct results; first(gene_switch_q_value) is not a valid combined
   statistic. Group by gene plus comparison. nrSwitches currently counts
   significant isoform rows, not independent isoform pairs or splicing events.
8. The volcano sig flag checks q only, while the significant table also requires
   dIF. Use separate FDR and effect-filter labels consistently.
9. Both require human hg38 and predict the longest ORF for every analyzed
   transcript. For annotated transcripts, preserve the annotation's CDS; predict
   missing ORFs separately. A predicted ORF is not proof of translation, and
   predicted NMD sensitivity is not an experimentally established NMD effect.
10. Optional consequence and GO dependencies are loaded eagerly, so an absent
    genome/GO package prevents even the basic isoform-usage test. Separate core
    testing, annotation, consequences and enrichment into optional stages.
11. Sequence export defaults can write outside outdir. Explicitly set its output
    path/prefix and document saved nucleotide/amino-acid files.
12. Restricting to standard chromosomes is a policy choice, not a universal
    biological requirement. Make it explicit and report discarded transcripts.
    Never strip transcript versions or suffixes without collision checks.

The reviewed upstream DEXSeq wrapper sets `gene_switch_q_value` to the minimum
isoform-adjusted p-value within that gene/comparison. Preserve it as a package
score, not an independently calibrated gene-level FDR.

## Version compatibility

The current IsoformSwitchAnalyzeR 2.12 manual documents preFilter with isoCount,
min.Count.prop, IFcutoff and min.IF.prop; it no longer lists the two
expression-cutoff arguments used by both originals. This is a concrete reason
not to assume either old script runs unchanged with a current package. Pin a
tested package version and adapt the filter interface with recorded semantics.

The comments saying countsFromAbundance='no' is required for DEXSeq are too
strong. The current package importer uses scaledTPM. Select the import/count
strategy using that package's documented workflow and account for effective
length; do not silently change the published analysis's input semantics. This
review does not establish that the original study's results were wrong.

## Proposed integration

Use a separate `IsoformSwitchAnalysis.R` module, callable both independently and
from the main workflow. Keep its Salmon quantifications, matched GTF/transcript
FASTA, comparison definitions and package dependencies explicit. The minimum
core should export all tested isoforms, IF/dIF, adjusted p-values, sample QC and
per-comparison gene summaries. Consequence prediction, GO and switch plots can
be separately requested. It should not be automatically treated as equivalent
to the StringTie assembly/quantification module.

The dedicated module is now implemented as `IsoformSwitchAnalysis.R` with
`IsoformSwitchHelpers.R`, and integrated under `isoform_switch`. See
[usage and method decisions](ISOFORM_SWITCH.md). Helper/import/annotation tests
passed; full IsoformSwitchAnalyzeR execution is pending a compatible installation.

Reference checked: [IsoformSwitchAnalyzeR manual](https://bioconductor.org/packages/release/bioc/manuals/IsoformSwitchAnalyzeR/man/IsoformSwitchAnalyzeR.pdf).
