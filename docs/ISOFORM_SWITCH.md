# Isoform-switch analysis

`IsoformSwitchAnalysis.R` and `IsoformSwitchHelpers.R` consolidate the useful
Salmon → tximport → IsoformSwitchAnalyzeR → DEXSeq workflow in the original
`iso_pipeline.R` and `iso_pipeline_New.R`. Keep the two new files together.
Neither file runs an analysis merely because it is sourced.

This module tests differential **within-gene isoform usage**. IF is an isoform's
fraction of its gene's expression; dIF is treatment IF minus reference IF.
A significant usage change is not necessarily a reversal of the dominant isoform,
an event-level PSI change, or differential gene expression. It does not reconstruct
transcripts or run Salmon; the separate `isoforms` module handles assembly.

## Run the small example

Use a coherent current R/Bioconductor installation. The module requires the
count-based `IsoformSwitchAnalyzeR::preFilter` interface documented in version
2.12; the older TPM-based interface is rejected rather than silently reinterpreted.
Exact installed package versions and session information are saved with each run.

```bash
Rscript install_dependencies.R --isoform-switch
Rscript IsoformSwitchAnalysis.R examples/isoform_switch/config.json --check
Rscript IsoformSwitchAnalysis.R examples/isoform_switch/config.json
```

The example has 60 synthetic two-isoform genes and six samples, including five
planted usage changes. Its Salmon-shaped quant.sf files are generated test inputs,
not the output of actual read mapping. The main package analysis has **not yet
been executed in the local development environment**; see the validation page.
Helper tests really import these files with tximport and check the GTF/FASTA.

## Minimum real inputs

* One Salmon `quant.sf` per independent biological sample, all from the same index.
* A TSV with `sample_id`, `condition`, `quant`. Relative quant paths resolve against
  the TSV directory. Additional covariates require explicit types in the config.
* The matching GTF with transcript, exon, gene and (for coding consequences) CDS
  annotations. Preserve original IDs and coding records.
* Explicit pairwise comparisons and a new output directory.

At least two biological replicates per condition and a nonsaturated identifiable
model are required. Three or more replicates generally provide a stronger basis
for estimating biological variability. Extra technical quantifications must be
combined appropriately, not mislabeled as biological replicates.

The config paths resolve relative to the config file. Example real configuration:

```json
{
  "samples": "samples.tsv",
  "gtf": "reference/annotation.gtf",
  "transcript_fasta": "reference/transcripts.fa",
  "out": "isoform_results",
  "comparisons": [{"name": "KO_vs_WT", "reference": "WT", "treatment": "KO"}],
  "covariates": {"batch": "factor"},
  "alpha": 0.05,
  "delta_if": 0.1,
  "strip_pipe": true,
  "filter": {"iso_count": 10, "count_proportion": 0.7,
             "if_cutoff": 0.01, "if_proportion": 0.5},
  "consequences": true,
  "predict_novel_orfs": false,
  "plots": 10
}
```

`transcript_fasta` is optional for core testing but required for consequences.
It must contain **spliced transcript sequences**, not a whole genome. IDs and
exon-spliced lengths are checked. `strip_pipe: true` takes the text before the
first `|` in Salmon/FASTA identifiers, useful for GENCODE headers. Collisions fail;
dot/version suffixes are never removed. Do not mix different annotation releases.

## Statistical and filtering choices

The import uses `tximport(..., txOut=TRUE, countsFromAbundance="scaledTPM")`,
following the current IsoformSwitchAnalyzeR importer strategy. These are
library-scaled abundance-derived counts, not raw TPM supplied to DEXSeq.
DEXSeq receives rounded estimated counts inside the upstream package. Effective
lengths, counts and TPM are exported. This is an explicit change from the old
scripts' `countsFromAbundance="no"`; it is not proof the published study was wrong.
Salmon inferential replicate uncertainty is not modeled in this workflow.

Count/IF filters use the upstream current API. `count_proportion` and
`if_proportion` control the proportions of samples satisfying the respective
criteria; these are not the old gene/isoform TPM thresholds. Single-isoform genes
are removed. The retained-ID audit records import and testing inclusion. All
comparisons retain a common pool through `keepIsoformInAllConditions=TRUE`.

Covariates are explicitly typed `factor` or `numeric`; batch/subject columns are
retained. The wrapper checks each two-group design for confounding and residual
degrees of freedom. It rejects sparse numeric covariates that the upstream DEXSeq
wrapper would silently recast as categorical. Automatic unwanted-factor discovery
is disabled. This wrapper uses the package's pairwise usage model; it does not
accept arbitrary interaction contrasts like the native DGE/DEJU modules do.
Package IF/dIF effect estimates may incorporate its confounder adjustment; the
raw TPM exports and PCA remain unadjusted descriptive QC.

All tested results and the full core RDS are saved before consequence filtering.
Significant means package isoform q-value **< alpha** and **|dIF| > delta_if**,
matching package strict thresholds. The effect filter is descriptive, not a
separate composite-null effect-size test. Results retain package p/q values;
multiple comparisons are separate analysis families, not a global FDR claim.
Gene summaries are per gene **and comparison**, count significant isoform rows,
and preserve the package gene score; they do not label rows as independent
switch pairs. In the reviewed DEXSeq wrapper, `gene_switch_q_value` is the
minimum isoform-adjusted p-value within a gene/comparison. It is **not** a
separately calibrated gene-level FDR. Use isoform-level tests for significance,
and do not report the gene rollup as an independent gene-level test.

## Optional consequences, sequence export, plots and GO

`consequences: true` adds structure-based alternative-splicing annotations,
preserves annotated CDS/ORFs, extracts nucleotide/amino-acid sequences, and
compares NMD status, ORF sequence similarity and intron retention among candidate
switches. Sequence exports stay inside the output directory. Alternative-splicing
annotation includes intron-retention information; a duplicate IR analysis is not
needed. All tested core results remain separate from the consequence subset.

`predict_novel_orfs: true` additionally predicts ORFs only for transcripts marked
not-yet-annotated by the package; it does not replace known coding regions or
assign an ORF to every annotated noncoding transcript. This requires usable
annotated ORFs in the matched reference. Predicted coding/NMD consequences are
hypotheses, not proof of translation or RNA decay. Domains and signal peptides
are not inferred without their corresponding external evidence.

`plots` controls the number of top switch-gene plots per comparison; 0 disables
them. With no significant isoforms, downstream switch-specific stages are
recorded as having no significant switches, rather than crashing on empty data.
An optional-stage error preserves the saved core output but leaves the overall
run INCOMPLETE.

To request GO, add:

```json
"go": {"orgdb": "org.Hs.eg.db", "mapping": "gene_to_entrez.tsv", "ontology": "BP"}
```

The mapping TSV must have exact `gene_id` and `entrez_id` columns. Select the
organism database explicitly. Ambiguous one-to-many gene mappings are excluded
and reported instead of selecting an arbitrary first match. Each comparison's
mapped tested genes form its universe. All returned GO terms are exported with
adjusted p-values; apply the chosen threshold when interpreting them. This
optional stage requires clusterProfiler, AnnotationDbi and the named OrgDb.

## Outputs

Core outputs include input checksums/config/package versions, sample/design and
comparison tables, transcript-to-gene mapping, imported count/TPM/effective-length
matrices, filtering audit, sample QC/PCA, the package statistical table (including p-values),
all tested and significant isoform tables,
per-comparison gene summaries, volcano PDFs, and the all-tested RDS. Optional
stages add annotated switch RDS, ORF/splicing/consequence tables, nucleotide and
amino-acid FASTAs, switch plots, GO results and mapping/universe audits.
`stage_status.tsv`, `sessionInfo.txt` and COMPLETE/INCOMPLETE markers report what
actually finished. Output directories are never silently reused.

## Integrated invocation

Use the same sample sheet for both workflows (include `quant` and, where needed,
`bam`), and reference the standalone switch config:

```json
"modules": ["dge", "dju", "isoform_switch"],
"isoform_switch": {"config": "isoform_switch.json"}
```

These are members of the main configuration, not a complete JSON document.
The orchestrator overrides the switch output to `RUN/isoform_switch`, resolves
its paths, checks shared sample-sheet/GTF identity, runs dependency/design
preflight, and includes its tables/PDFs/status in the main HTML report. Standalone
and integrated execution call exactly the same R module.

Reference: [IsoformSwitchAnalyzeR manual](https://bioconductor.org/packages/release/bioc/manuals/IsoformSwitchAnalyzeR/man/IsoformSwitchAnalyzeR.pdf).
