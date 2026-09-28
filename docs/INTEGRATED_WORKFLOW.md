# Integrated transcriptome workflow, v0.2.0

Run `python3 TranscriptomePipeline.py --config analysis.json`. Python >=3.9 uses
only its standard library. Native modules use Rscript and Bioconductor packages
installed by `Rscript install_dependencies.R`. The original DEJUPipeline.R entry
point is retained. Do not overwrite existing output directories.

`--plan` prints resolved configuration without running or validating executable
availability. `--check` checks configuration, requested executable availability,
native model design/dependencies and supplied count tables without generating
analysis outputs. Neither is a completed biological analysis.

## Modules and evidence

| Module | Function | Required evidence/tools |
|---|---|---|
| `dge` | DESeq2 gene-level differential expression | Raw gene counts, or BAM + GTF; DESeq2 |
| `dju` | Differential exon/junction usage | Raw feature counts + annotation, or BAM + GTF; edgeR/limma |
| `events` | SE/A5SS/A3SS/MXE/RI classification, per-sample PSI and delta PSI | BAM + GTF, read length, explicit two-group comparisons; rMATS-turbo |
| `isoforms` | Assembly, common transcript reference, re-quantification and candidate classification | Coordinate-sorted BAM + GTF; StringTie + GffCompare |
| `isoform_switch` | Salmon/DEXSeq differential isoform usage and optional consequences | Salmon quant.sf + matched GTF, explicit comparisons; IsoformSwitchAnalyzeR; transcript FASTA for consequences |
| `fusions` | Fusion candidate detection | Original STAR BAM with chimeric evidence, genome FASTA, GTF and matching blacklist; Arriba |
| `loci` | Per-sample raw coverage and junction-arc plots with annotated exons | Coordinate-sorted BAM + GTF and explicit loci; Rsamtools/GenomicAlignments |

An ordinary gene-count matrix supports DGE, not reconstruction, fusion detection,
isoform switching or event PSI. Junction-relative logFC is not delta PSI. Missing
analyses are not fabricated or silently skipped. The report records requested,
completed, failed and unrequested modules.

## Start with the count-only example

`examples/integrated/config.json` uses synthetic data with two conditions, a
continuous time variable, a categorical batch and a condition-by-time interaction.
It includes baseline, later-time and interaction contrasts. The example output
path resolves relative to that configuration. Change it for each new run.

Minimal real-count configuration:

```json
{
  "samples": "samples.tsv",
  "gene_counts": "gene_counts.tsv",
  "out": "results",
  "modules": ["dge"],
  "design": {
    "formula": "~ condition",
    "factors": {"condition": ["WT", "KO"]},
    "numeric": []
  },
  "contrasts": [{"name": "KO_vs_WT", "weights": {"conditionKO": 1}}],
  "fdr": 0.05
}
```

Paths in JSON resolve relative to the JSON file. BAM paths in the sample TSV
resolve relative to that TSV. Sample columns must contain `sample_id` and every
variable in the design. Count-table columns must exactly match sample IDs; an
optional `count_column` maps each sample to the original featureCounts column
name. The first count-table column contains unique gene IDs. CSV and TSV, including
gzip, are supported. Standard featureCounts annotation columns and comment lines
are accepted. Unmapped sample columns are rejected. Counts must be nonnegative
integers; TPM, FPKM, VST and library-normalized values are not raw gene counts.

## Model specification

Declare every model variable as categorical (`design.factors`, with explicit
level order) or continuous (`design.numeric`). Formula operators support additive
terms, interactions, no-intercept coding and ordinary fixed-effect expansions;
arbitrary R function calls are not accepted. Multiple condition levels are
supported. This is a fixed-effect model, not a random-effects model.

Contrasts are named numeric weights on exact `model.matrix` coefficient names.
The `--check` output lists these names, as does `native/design.tsv`. With
`~ batch + condition * time`, the condition coefficient is its effect at time
zero; the interaction coefficient is the change in time slope. To compare KO
versus WT at time t, combine `conditionKO: 1` with `conditionKO:time: t` in the
weights object (see the executable example). Continuous covariates are not
silently centered/scaled. Decide meaningful units and zero points beforehand.

DGE and DEJU use the same supplied contrasts. FDR is calculated separately for
each contrast and analysis family, not globally across every module/comparison.
Rank-deficient designs and designs with no residual degrees of freedom fail.
Sufficient independent biological replication remains the investigator's
responsibility; technical lanes must not be treated as biological replicates.

## BAM-based DGE and DEJU

Remove `gene_counts` and/or `dju_counts` when the corresponding module should
count BAMs. Supply `gtf`, `bam` and a `bam` column in the sample sheet:

```json
"bam": {"paired": true, "strand": 2, "min_mapq": 10, "remove_duplicates": false},
"threads": 8,
"dju": {"engine": "modern", "min_anchor": 8, "min_intron": 20,
        "max_intron": 1000000, "min_count": 10, "min_total_count": 15},
"dge": {"min_count": 10, "min_samples": 2, "fit_type": "parametric"}
```

These are configuration members, not a complete JSON document. Strand is
0=unstranded, 1=forward, 2=reverse. `paired` is required. One library orientation
applies to all samples; mixed library types require separate handling.

Native BAM filters retain mapped, primary, non-supplementary, QC-passing
alignments with `NH:i:1` and adequate MAPQ. Duplicate removal uses existing flags;
it does not discover duplicates. `native/alignment_qc.tsv` reports total/passing
alignments, missing NH, duplicate flags and other categories. These categories
can overlap and must not be summed as mutually exclusive exclusions.

DGE makes a **separate gene-level featureCounts pass** using all eligible spliced
and non-spliced reads and gene meta-features. Paired-end DGE counts fragments,
requires both ends mapped, excludes chimeric fragments, and does not impose a
genomic fragment-length cutoff. Never sum the DEJU exon+junction matrix to obtain
DGE counts. DEJU counts aligned reads and each accepted N gap, with conservative
contiguous splice-anchor and intron-length filters. Soft clips/insertions/deletions
interrupt that contiguous anchor. The count units are recorded here deliberately.

DESeq2 keeps raw integer input, estimates size factors and dispersion, and uses
Wald contrasts. Low-count filtering and DESeq2 independent filtering are distinct.
Automatic outlier count replacement is disabled; Cook's-distance evidence is
exported and default Cook's testing exclusion remains enabled. Unavailable
p-values remain NA. A failed dispersion fit remains an error, rather than silently
switching statistical methods. Log2FC values are unshrunk and labeled accordingly.

## Event classification and PSI

Add `events` and configure:

```json
"events": {
  "read_length": 150,
  "variable_read_length": false,
  "min_anchor": 8,
  "delta_filter": 0.1,
  "comparisons": [{"name": "KO_vs_WT", "column": "condition",
    "treatment": "KO", "reference": "WT", "model": "unpaired",
    "dju_contrast": "KO_vs_WT"}]
}
```

Use `model: "paired"` and `subject: "subject_column"` for biological matched
pairs; the wrapper sorts and validates matching subject sets. Read pairs and
matched biological subjects are different. rMATS models **do not inherit** the
DESeq2/DEJU design. They do not adjust an arbitrary batch/continuous-interaction
design here. Its outputs are explicitly two-group event analyses; do not report
them as adjusted effects from the general design.

The wrapper uses treatment as rMATS b1, so delta PSI is treatment minus reference.
It preserves both JC and JCEC outputs, recomputes BH across all five event types
separately for each count type/comparison, validates reported PSI against inclusion
and exclusion counts/effective lengths, and exports replicate PSI. Pick a primary
count type before interpretation. `passes_delta_filter` is a descriptive filter,
not a formal test against that effect-size bound. Intronic retention evidence is
provided by rMATS; missing junction reads alone are not treated as retention.

Where requested, DEJU junctions are linked to event structures by gene, strand
and genomic endpoints. A junction can belong to multiple events; unmatched
junctions remain unclassified. An associated significant DEJU junction does not
make the event significant, or vice versa. rMATS/MXE inclusion labels follow the
caller's event definition, not a claim that one exon is universally dominant.

rMATS uses its own alignment policies on original BAMs. Its read-length setting
influences effective-length PSI; `variable_read_length` does not make that setting
irrelevant. Check its raw logs/read outcomes and match settings to the experiment.

## Transcript isoforms and novel loci

`isoforms` runs StringTie per sample, merges assemblies across samples, compares
the merged annotation with GffCompare, and re-quantifies every sample against the
same merged transcript reference. It exports per-transcript TPM and the original
quantified GTF/gene abundance files. Missing transcript records stay NA, not
invented zero counts. This is assembly/quantification, not differential isoform
usage testing.

GffCompare class `u` is exported as an **intergenic transcript / novel gene
candidate**. Other novel transcript structures retain their class codes.
Short-read assembly does not prove a full-length structure, an independent gene,
translation or biological function. Validate candidate loci and structures
before making those claims. TPM is never sent directly to DESeq2.

## Isoform-switch analysis

Add `isoform_switch` with `"isoform_switch": {"config": "isoform_switch.json"}`.
It invokes the standalone Salmon/DEXSeq module using the same sample sheet and
reference GTF. Its output path is set to the main run's `isoform_switch/` folder.
The explicit switch comparisons/covariates are separate from native DGE/DJU
contrasts; they are never silently copied between models. See
[the module guide](ISOFORM_SWITCH.md) for inputs, outputs and limitations.

## Fusion candidates

Add `fusions`, `fasta` and:

```json
"fusions": {
  "star_chimeric_bam": true,
  "blacklist": "arriba_matching_assembly_blacklist.tsv.gz",
  "known_fusions": "known_fusions.tsv.gz"
}
```

`star_chimeric_bam: true` declares the required upstream preparation; it does
not transform an ordinary BAM into a fusion-ready alignment. Follow Arriba's
STAR alignment instructions. Original BAMs, including chimeric/supplementary
alignment evidence, go directly to Arriba. Reference assembly/resources must
match. Candidates, caller confidence and discarded calls are retained. These
are not confirmed genomic rearrangements or established novel genes.

## Locus visualizations

```json
"loci": [{"name": "example_gene", "chr": "chr1", "start": 100000,
          "end": 110000, "strand": "+"}]
```

Each region produces a PDF containing sample-level raw coverage, junction arcs
and merged gene exons. Coordinates are 1-based, inclusive; intervals are limited
to 1 Mb. Native alignment filters apply. Optional locus strand restricts reads
when library strandedness is known. Arc labels are local raw junction support;
they are not the filtered differential-testing count matrix and not PSI. The
module indexes BAM aliases in its output folder without changing source BAMs.

## External executable configuration

Install rMATS-turbo, StringTie, GffCompare and Arriba on an appropriate analysis
system (typically Linux/conda/container). No automatic system installation or
reference downloads occur. Override executables with argument arrays:

```json
"tools": {
  "Rscript": ["Rscript"],
  "rmats": ["python3", "/opt/rmats/rmats.py"],
  "stringtie": ["stringtie"],
  "gffcompare": ["gffcompare"],
  "arriba": ["arriba"]
}
```

All commands run as argument arrays, without a shell. Missing required tools fail
preflight. Stage commands/stdout/stderr are saved in logs. Native R session
versions, input checksums, BAM size/mtime manifest and resolved configuration
are saved. Full BAM contents are not hashed. An output directory's COMPLETE.txt
is written only when every requested module succeeds; failures retain
INCOMPLETE.txt, stage status and partial outputs for diagnosis.

## Main outputs

* `report.html`, `run_manifest.json`: status, result links, analysis caveats.
* `dge/`: raw-input filter decisions, normalized/VST counts, PCA, sample
  correlations, dispersion/MA/volcano plots, Cook's distances, model RDS and a
  table per contrast (`baseMean`, log2FC, SE, statistic, p-value, adjusted p-value).
* `native/`: sample alignment QC, independently counted gene matrix, gene
  assignment statistics, model matrix/contrast weights and R session versions.
* `junction_counts/`, `dju/<contrast>/`: feature counts, assignment audit,
  junction/exon/gene usage tests, diagnostic PDF, model and transcript positions.
* `events/<comparison>/`: original caller results, event classes, PSI by sample,
  delta PSI, joint event-type FDR, sample order and optional DEJU-event links.
* `isoforms/`: assemblies, merged reference, quantified transcripts, TPM matrix,
  structural classifications and novel-locus candidates.
* `isoform_switch/`: complete tested isoform results, IF/dIF, gene summaries, QC,
  and requested sequence/consequence/GO/plot outputs.
* `fusions/`: per-sample calls, discarded calls and combined candidate table.
* `loci/`: locus PDFs.

## Isoform-switch analysis

The original Salmon/DEXSeq/IsoformSwitchAnalyzeR scripts have been reviewed in
[ISOFORM_SWITCH_REVIEW.md](ISOFORM_SWITCH_REVIEW.md). A standalone module callable
from this orchestrator is recommended. It is deliberately not enabled during
this inspection-first step. Salmon quantifications and their matching transcript
annotation are needed; a gene-count matrix cannot recover isoform fractions.

## Method sources

[DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html),
[rMATS-turbo](https://github.com/Xinglab/rmats-turbo),
[StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual),
[GffCompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml),
[Arriba](https://github.com/suhrig/arriba/wiki).

## Reference resources

See [reference downloads](../references/README.md) for the verified GENCODE v48
manifest, matching-assembly guidance, checksum-verifying fetcher and third-party
licensing/citation guidance. Full genomes are downloaded explicitly, outside Git.
