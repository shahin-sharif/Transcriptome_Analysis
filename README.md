# Transcriptome_Analysis

## Integrated workflow (0.2.0)

`TranscriptomePipeline.py` is the new configurable entry point. It supports
DESeq2 gene-expression analysis, general-design DEJU, rMATS event/PSI analysis,
StringTie/GffCompare transcript reconstruction and quantification, Arriba fusion
candidates, locus plots and an HTML run report. Native R functions and external
caller adapters have different validation levels; see the validation document.

* [Configuration, modules and input requirements](docs/INTEGRATED_WORKFLOW.md)
* [Validation and remaining limits](docs/VALIDATION_V0.2.md)
* [Salmon/DEXSeq isoform-switch analysis](docs/ISOFORM_SWITCH.md)
* [Original isoform pipeline review](docs/ISOFORM_SWITCH_REVIEW.md)
* [Reference downloads, checksums and reuse guidance](references/README.md)

`IsoformSwitchAnalysis.R` is a dedicated Salmon/DEXSeq isoform-switch workflow,
callable independently or through `isoform_switch` in the main pipeline. It
retains all tested results and adds optional consequence prediction, plots and GO.
It is distinct from StringTie transcript reconstruction/quantification.

**Validation:** native DESeq2/DEJU count and BAM tests passed. External callers
and the complete IsoformSwitchAnalyzeR analysis still require execution testing
in an environment with those dependencies. No real experimental dataset has
been validated for this expanded release. See the validation document.

```bash
Rscript install_dependencies.R
python3 TranscriptomePipeline.py --config examples/integrated/config.json --check
python3 TranscriptomePipeline.py --config examples/integrated/config.json
```

The original standalone DEJU command remains available below. New junction
counting defaults require 8 contiguous aligned bases on each side of a splice
gap and intron lengths of 20..1,000,000 bases. The integrated config exposes
these thresholds; set them appropriately for the organism and alignment policy.

## DEJU: differential exon and junction usage

`DEJUPipeline.R` and `DEJUHelpers.R` consolidate the useful parts of the previous
DEJU pipelines into one command-line workflow. The model combines counts from
non-spliced exon reads and splice-junction reads, then tests whether a feature's
change differs from the other retained features of its gene.

This measures **relative usage**, not PSI, absolute expression, or differential
transcript abundance. A significant gene-level test can be driven by an exon;
it does not necessarily identify a significant junction.

DEJU version: **0.2.0**. Both files must remain together. Sourcing either file does not
launch an analysis. See [CHANGES.md](CHANGES.md) for the consolidation decisions
and [tests/README.md](tests/README.md) for what has been tested.

## Install

Use a recent R installation with a matching Bioconductor release. The modern
engine requires edgeR >= 4.6 and a compatible limma. The optional legacy engine
uses `diffSpliceDGE` with legacy QL fitting; the pipeline requires edgeR >= 4.0.
Do not mix Bioconductor releases in an existing analysis environment.

```bash
Rscript install_dependencies.R
```

To include optional GO enrichment packages:

```bash
Rscript install_dependencies.R --go
```

The installer uses the library configured for your R session. Use `R_LIBS_USER`
to select a separate library if desired. Installation needs internet access;
analysis itself does not download annotations or packages.

## Try the small example

From the repository directory:

```bash
Rscript DEJUPipeline.R \
  --samples examples/counts/samples.tsv \
  --counts examples/counts/counts.tsv \
  --features examples/counts/features.tsv \
  --reference WT --treatment KO \
  --covariates batch \
  --out example_results
```

These are synthetic counts for 120 genes, four features per gene and eight
samples. G001 has an increased junction-usage signal, G002 a decreased signal,
and G003 an overall expression increase without a designed usage change. They
are software fixtures, not biological observations. Always use a new output
folder: the pipeline refuses to overwrite an existing path.

## Run on BAM files

Create a tab-separated sample sheet, with one row per **biological replicate**:

```text
sample_id	condition	bam	batch
WT1	WT	bams/WT1.bam	A
WT2	WT	bams/WT2.bam	B
WT3	WT	bams/WT3.bam	A
KO1	KO	bams/KO1.bam	A
KO2	KO	bams/KO2.bam	B
KO3	KO	bams/KO3.bam	A
```

Replace the example text with real tab-separated columns. Relative BAM paths
are resolved against the sample-sheet directory. The GTF must match the genome
assembly and chromosome names used for alignment and must contain quoted
`gene_id` and `transcript_id` attributes on exon records. Versioned IDs are
preserved. No human-specific or chromosome-prefix assumptions are made.

```bash
Rscript DEJUPipeline.R \
  --samples samples.tsv \
  --gtf gencode.v48.annotation.gtf \
  --paired true \
  --strand 2 \
  --reference WT --treatment KO \
  --covariates batch \
  --threads 8 \
  --out DEJU_WT_vs_KO
```

Choose `--paired true` or `false` from the sequencing design. Choose `--strand 0`
for unstranded, `1` for forward/sense, or `2` for reverse/antisense libraries.
For paired libraries, strand refers to the transcript orientation relative to
read 1; read 2 has the opposite alignment orientation. These values are required,
not guessed. Mixed single/paired libraries are rejected. Omit `--covariates`
when there is no batch/subject adjustment. Covariates are **categorical**,
comma-separated sample-sheet columns. A paired biological design can use
`--covariates subject`; sequencing read pairs and paired subjects are different
concepts. The model is additive, without interactions.

Exactly two conditions and at least two biological replicates per condition are
required. More replication is preferable. Rank-deficient designs and designs
without residual degrees of freedom are rejected. Technical lanes must be
combined appropriately, not listed as independent biological replicates.
Positive usage effects mean **treatment minus reference**.

### Counting policy

* All retained alignments must be mapped, primary, non-supplementary, QC-passing,
  have mapping quality >= `--min-mapq` (default 10), and carry **`NH:i:1`**.
  Alignments without NH are excluded, not assumed unique. Confirm that the aligner
  writes this tag. A MAPQ value of 255 is accepted; its meaning is aligner-specific.
* Duplicate-marked alignments are retained by default. `--remove-duplicates true`
  excludes marked duplicates; it does not discover duplicates or perform UMI
  deduplication. Choose this based on the experiment and upstream processing.
* The counting unit is an **aligned read**, including for paired-end data.
  Overlapping mates can both contribute. No fragment-length limit or
  requirement for both mates to be mapped is imposed. This policy is explicit
  and differs from a fragment-counting analysis.
* Exons are merged within each gene. Rsubread counts non-spliced reads against
  these flattened exon features. Reads ambiguously overlapping multiple features
  are excluded from exon counts. Spliced reads do not enter the exon channel.
* Junctions are counted in chunks from CIGAR `N` operations, after the same BAM
  filters. Each read contributes once to each intron it spans; a read crossing
  several introns contributes to several junction features. Strand is inferred
  from flags, mate identity and the declared library orientation.
* Junction coordinates are the **1-based exonic bases flanking the intron**:
  `left < right` on both strands. Donor/acceptor direction is separate.
* Annotated junctions are assigned when the GTF provides one compatible gene.
  Novel junctions must lie entirely within exactly one compatible gene span.
  Ambiguous and unassigned junctions are excluded and reported. Strand is used
  for stranded libraries. Unstranded data cannot resolve opposite-strand genes
  sharing a junction. Fusion/trans-splicing analysis is outside this workflow.

Direct CIGAR counting avoids two observed behaviors in Rsubread 2.10.5:
`nonSplitOnly=TRUE` suppresses junction collection, and junction totals can be
independent of `strandSpecific`. The automated BAM fixtures check the chosen
counting behavior instead of relying on those defaults.

Temporary filtered BAMs are stored under the output directory and removed when
the counting function exits. Plan for sufficient disk space. Interrupted process
termination may leave temporary files; an incomplete run is never marked complete.

### Statistics and thresholds

1. Align sample columns by sample IDs and validate raw integer counts.
2. Apply design-aware `filterByExpr` (default minimum count 10, total count 15).
3. Retain genes with at least two expressed features and at least one expressed
   junction. A retained gene need not retain an exon after expression filtering.
4. Recalculate library sizes, apply TMM normalization, and estimate robust
   negative-binomial dispersions and quasi-likelihood models.
5. Run the explicitly selected usage engine: `modern` (default, `diffSplice`) or
   `legacy` (`diffSpliceDGE`). No silent engine fallback occurs. Modern and legacy
   p-values/effects need not agree exactly.
6. Adjust junction p-values by BH across **all tested junctions**. Separately
   report all-feature BH, exon-only BH, gene Simes BH, and gene F-test BH.

The default significance cutoff is `--fdr 0.05`. `usage_log2FC` is converted from
edgeR's underlying natural-log usage coefficient. It is a relative feature-usage
effect, not the gene-expression fold change. `--min-abs-log2 0` is the default
optional effect-size filter. Changing it only changes `passes_effect_filter`;
it is **not** a formal test against a nonzero effect-size threshold and does not
change p-values or FDR. `significant` always refers to the junction FDR alone.

`has_significant_junction` identifies membership in the junction-discovery list;
it is not another gene-level FDR result. Gene tests include both exon and junction
features. Select and report your primary testing family before interpreting
results; testing several families does not provide one combined FDR guarantee.

## Outputs

| File | Meaning |
|---|---|
| `COMPLETE.txt` | Written only when all requested steps succeed; concise run summary |
| `INCOMPLETE.txt` | Remains if a run fails; partial outputs are not a finished analysis |
| `counts.tsv.gz`, `features.tsv.gz` | Raw input feature counts and annotation |
| `junction_assignment.tsv.gz` | All observed junctions, assignment decisions and counts (BAM mode) |
| `featurecounts_assignment.tsv` | Exon-counting assignment summary **after** common BAM filtering |
| `filtering_results.tsv.gz` | Retained/excluded features and reasons |
| `junctions_results.tsv.gz` | Usage effects, p-values, all-feature FDR, junction FDR and significance |
| `exons_results.tsv.gz` | Exon effects and separate exon FDR |
| `features_results.tsv.gz` | All tested features with all-feature FDR |
| `genes_results.tsv.gz` | Gene Simes/F tests and significant-junction membership |
| `junction_transcript_positions.tsv.gz` | All compatible transcript positions, when GTF and coordinates are available |
| `flattened_exons.saf.tsv.gz`, `annotated_junctions.tsv.gz` | Derived annotation used by the run |
| `normalization.tsv`, `design.tsv` | Library sizes, TMM factors and exact design matrix |
| `diagnostics.pdf` | MDS, dispersion and junction volcano plots |
| `analysis.rds` | Settings, samples and fitted models for further R analysis |
| `settings.R`, `sessionInfo.txt`, `input_checksums.tsv` | Settings, package versions and MD5 checksums of small inputs/GTF |
| `bam_manifest.tsv` | BAM paths, sizes and modification times; BAMs are not fully hashed |

Transcript positions run from 5′ to 3′, using `(junction_rank - 0.5) / n_junctions`.
A junction may appear in several transcripts and therefore several rows. **Do
not count those rows as independent junction discoveries.** Novel junctions have
missing transcript positions; no arbitrary transcript is assigned.

## Optional GO enrichment

Supply an explicit mapping of the pipeline's full `gene_id` values to `ENTREZID`,
plus the matching organism database:

```bash
# Append to the BAM or count-table command:
  --gene-map gene_to_entrez.tsv --go-orgdb org.Hs.eg.db
```

The mapping file needs `gene_id` and `ENTREZID` columns. No automatic Ensembl
version stripping or species guessing occurs. One-to-many mappings are allowed
and deduplicated at ENTREZID level. Check mapping coverage before interpretation.
The background is the mapped set of **eligible, tested genes**, restricted further
by the annotation available to GO. Biological Process enrichment runs separately
for gene-Simes discoveries and genes containing significant junctions. Tables
include BH-adjusted GO p-values; the two analyses must not be treated as
independent confirmation. Universe/selection files document the mapping used.
GO uses FDR-significant junctions, not the optional effect-size-filtered list.
An empty discovery list yields an empty GO table.

## Count-table input and reuse

`--counts` requires a TSV whose first column is the unique feature ID and whose
remaining columns are exact sample IDs. Values must be raw, finite, nonnegative
integers. `--features` requires `feature_id`, `gene_id`, and `feature_type`
(`exon` or `junction`); optional `chr`, `left`, `right`, `strand` enable transcript
position annotation with `--gtf`. Files may be gzip-compressed. Counts and features
must contain the same IDs; sample order is checked and aligned explicitly.
BAM filtering options are rejected in this mode because they cannot change
precomputed counts. The caller is responsible for consistent annotation and
counting policy in supplied matrices.

```bash
Rscript DEJUPipeline.R --help
Rscript tests/test_pipeline.R
Rscript tests/test_pipeline.R --bam
```

## Provenance and limits

The five original pipelines share a core method but differ in validation and
junction annotation handling; none can be certified as the version previously
used successfully. A read-only marker scan of the saved R workspace found count
and model-related objects, but no reliable embedded pipeline identity or run
provenance. The original files and workspace have not been changed or uploaded.

Synthetic tests check implementation behavior; they do not establish biological
validity on the user's full experiment. Review library preparation, NH tagging,
strandedness, duplicate handling, batches, replicate identity, annotation version
and diagnostics before drawing conclusions. This first release focuses on a
single two-condition comparison. It is not an exact reproduction of every default
in the original scripts or the published DEJU workflow.

Method references: [DEJU paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12288301/),
[authors' workflow](https://github.com/TamPham271299/DEJU),
[edgeR manual](https://bioconductor.org/packages/release/bioc/manuals/edgeR/man/edgeR.pdf),
[Rsubread manual](https://bioconductor.org/packages/release/bioc/manuals/Rsubread/man/Rsubread.pdf),
[Rsamtools manual](https://bioconductor.org/packages/release/bioc/manuals/Rsamtools/man/Rsamtools.pdf).
