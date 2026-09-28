# Transcriptome_Analysis

**An RNA-seq toolkit for studying gene expression, exon and junction usage,
alternative splicing, and transcript isoforms.**

Transcriptome_Analysis brings complementary analyses into one configurable
workflow. It can start from aligned reads, existing count tables, or Salmon
transcript quantifications, depending on the analysis selected. Results include
statistical tables, quality-control plots, and a report linking the outputs.

The workflow is modular: run only the analyses needed for an experiment.
Gene-expression changes, junction-usage changes, event PSI and isoform fractions
answer different questions and are reported separately.

## Contents

1. [Choose an analysis](#1-choose-an-analysis)
2. [Download and install](#2-download-and-install)
3. [Run a complete example](#3-run-a-complete-example)
4. [Analyze your own BAM files](#4-analyze-your-own-bam-files)
5. [Analyze an existing gene-count table](#5-analyze-an-existing-gene-count-table)
6. [Run standalone DEJU](#6-run-standalone-deju)
7. [Analyze isoform switching with Salmon](#7-analyze-isoform-switching-with-salmon)
8. [Add other analyses](#8-add-other-analyses)
9. [Find and interpret results](#9-find-and-interpret-results)
10. [Troubleshooting and further documentation](#10-troubleshooting-and-further-documentation)

## 1. Choose an analysis

| Question | Analysis | Input | Main result |
|---|---|---|---|
| Which genes change expression? | `dge`: DESeq2 | Gene counts, or BAM files + GTF | Gene log2 fold changes and adjusted p-values |
| Which exon/junction features change relative to their gene? | `dju`: DEJU | Feature counts + feature annotation, or BAM files + GTF | Relative usage effects and FDR |
| Which splicing events change? | `events`: rMATS | BAM files + GTF + read-length information | Exon skipping, alternative donor/acceptor, mutually exclusive exons, intron retention; PSI and ΔPSI |
| Which transcript structures are supported? | `isoforms`: StringTie/GffCompare | BAM files + GTF | Reconstructed transcripts, abundance and novel-locus candidates |
| Which transcript isoforms change their share of gene expression? | `isoform_switch`: IsoformSwitchAnalyzeR/DEXSeq | Salmon quantifications + matching GTF | Isoform fractions, dIF and switch statistics |
| What do reads look like at a locus? | `loci` | BAM files + GTF + region coordinates | Coverage and junction-arc plots |
| Are fusion transcripts supported? | `fusions`: Arriba | Suitable STAR chimeric alignments + matched reference resources | Fusion candidates |

**Start here:** for gene expression and junction usage together, follow section 4.
For an existing gene-count table, follow section 5. For Salmon isoform switching,
follow section 7. Gene-count tables alone cannot recover event PSI or isoform usage.

## 2. Download and install

### Open a terminal and download the repository

The commands below are for Bash or Zsh on Linux/macOS. Run them in a terminal,
not inside the R console. Copy only the text inside each code block; do not add
a `$` prompt. Run commands from the repository directory unless stated otherwise.

```bash
mkdir -p "$HOME/projects"
cd "$HOME/projects"
git clone https://github.com/shahin-sharif/Transcriptome_Analysis.git
cd Transcriptome_Analysis
```

If the repository is private, your GitHub account needs access and Git must be
authenticated. If you already downloaded it, enter its existing directory rather
than cloning it again. Check your location and files:

```bash
pwd
ls
```

You should see `TranscriptomePipeline.py`, `DEJUPipeline.R`,
`IsoformSwitchAnalysis.R`, and the `examples` directory. Keep the repository's
files together: the entry scripts load helper files from the same repository.

### Check Python and R

```bash
python3 --version
Rscript --version
```

The integrated workflow needs Python 3.9 or newer. Python uses only its standard
library. R analyses require R and compatible Bioconductor packages. If either
command is missing, install that language before continuing; on a managed server,
use the software environment supplied by the administrator.

Install the core R dependencies:

```bash
Rscript install_dependencies.R
```

This installs packages, not R itself, and requires internet access. Use a current,
coherent R/Bioconductor environment. Modern DEJU requires edgeR ≥ 4.6 and compatible
limma. The installer does not automatically upgrade every existing package.

For isoform switching, also run:

```bash
Rscript install_dependencies.R --isoform-switch
```

For optional human GO enrichment, add `--go`. Other organisms need their matching
annotation package. The installer does **not** install rMATS, StringTie,
GffCompare, Arriba or Salmon; those are needed only for the corresponding analyses.
See [module requirements](docs/INTEGRATED_WORKFLOW.md#external-executable-configuration).

## 3. Run a complete example

This example uses the small synthetic gene-count table supplied with the
repository. No BAM files, reference genome or downloads of experimental data
are needed after installing the R dependencies.

First check the configuration, input files and dependencies:

```bash
python3 TranscriptomePipeline.py --config examples/integrated/config.json --check
```

Then run the analysis:

```bash
python3 TranscriptomePipeline.py --config examples/integrated/config.json
```

The example fits DESeq2 with condition, batch, time and a condition-by-time
interaction. It generates three comparisons. Its synthetic signals demonstrate
software behavior; they are not biological discoveries.

Find the results:

```bash
ls results/example_dge
ls results/example_dge/dge
cat results/example_dge/COMPLETE.txt
```

Open `results/example_dge/report.html` in a web browser. On macOS:

```bash
open results/example_dge/report.html
```

On a remote server, download the **whole output folder** to your computer before
opening the report so its links to result tables and PDFs continue to work.

**To repeat a run:** change `out` in the configuration to a new folder, for example
`../../results/example_dge_run2`. Existing output folders are not overwritten.
Paths in a JSON configuration are relative to that configuration's directory.

## 4. Analyze your own BAM files

This walkthrough runs **DESeq2 gene expression and DEJU junction usage together**
for three WT and three KO biological replicates. Replace these names and paths
with your experiment's actual samples. The example assumes paired-end,
reverse-stranded libraries; change those settings if your preparation differs.

### Step A — prepare the input files

You need:

* One aligned BAM file per biological sample.
* A GTF matching the genome assembly and chromosome names used for alignment.
* Known paired/single-end and library-strand settings.

Native BAM counting requires the aligner's `NH:i:1` tag for retained reads.
Reads without NH are excluded. Technical lanes are not independent biological
replicates. DESeq2 counts genes separately from DEJU features; it does not sum
junction counts to estimate gene expression.

### Step B — create the sample sheet

Create a folder for configuration files and open a new file with nano:

```bash
mkdir -p analysis
nano analysis/samples.tsv
```

Paste the following **tab-separated** table. Replace every `/absolute/path/...`
with a real path on your computer or server. Use actual tabs between columns,
not spaces. Do not include the code-block markers.

```tsv
sample_id	condition	bam
WT1	WT	/absolute/path/bams/WT1.bam
WT2	WT	/absolute/path/bams/WT2.bam
WT3	WT	/absolute/path/bams/WT3.bam
KO1	KO	/absolute/path/bams/KO1.bam
KO2	KO	/absolute/path/bams/KO2.bam
KO3	KO	/absolute/path/bams/KO3.bam
```

In nano, press **Ctrl+O**, then **Enter** to save, and **Ctrl+X** to exit.
`sample_id` identifies a sample; `condition` identifies its experimental group;
`bam` points to its alignment file. Relative BAM paths are resolved against the
sample-sheet directory. Absolute paths are easier to follow when starting out.

### Step C — create the analysis configuration

A JSON configuration is a text file containing the input paths, selected modules
and analysis settings. Create it:

```bash
nano analysis/bam_analysis.json
```

Paste this complete configuration, replacing the GTF path:

```json
{
  "samples": "samples.tsv",
  "gtf": "/absolute/path/reference/annotation.gtf",
  "out": "../results/WT_KO",
  "modules": ["dge", "dju"],
  "bam": {
    "paired": true,
    "strand": 2,
    "min_mapq": 10,
    "remove_duplicates": false
  },
  "threads": 4,
  "design": {
    "formula": "~ condition",
    "factors": {"condition": ["WT", "KO"]},
    "numeric": []
  },
  "contrasts": [
    {"name": "KO_vs_WT", "weights": {"conditionKO": 1}}
  ],
  "fdr": 0.05
}
```

Save and exit nano as above. JSON requires double quotes around text and no
trailing comma after the final item. Important settings:

| Setting | Meaning |
|---|---|
| `samples` | The sample sheet; here it is beside the JSON file |
| `gtf` | The matching annotation file; use an uncompressed GTF for BAM modules |
| `out` | A new results directory; this example creates `results/WT_KO` |
| `modules` | Analyses to run; choose `["dge"]` or `["dju"]` to run only one |
| `paired` | `true` for paired-end sequencing, `false` for single-end |
| `strand` | `0`: unstranded; `1`: forward; `2`: reverse. For paired reads, orientation is relative to read 1 |
| `threads` | CPU threads available for applicable counting/external steps |
| `formula` | Statistical model; `~ condition` compares groups without extra covariates |
| `conditionKO` | KO compared with the first factor level, WT |
| `fdr` | False-discovery-rate threshold |

Positive effects in this comparison mean **KO relative to WT**. DESeq2's gene
log2 fold change and DEJU's relative feature-usage effect are different measures.
For batch adjustment, matched subjects, more conditions or interactions, follow
[the model guide](docs/INTEGRATED_WORKFLOW.md#model-specification); adding a metadata
column alone does not add it to the statistical model.

### Step D — check, run and inspect

```bash
python3 TranscriptomePipeline.py --config analysis/bam_analysis.json --check
python3 TranscriptomePipeline.py --config analysis/bam_analysis.json
```

The first command checks the setup; it does not run the analysis. Run the second
only after correcting any reported errors. Check the completed output:

```bash
cat results/WT_KO/COMPLETE.txt
ls results/WT_KO/dge
ls results/WT_KO/dju/KO_vs_WT
```

The overview is `results/WT_KO/report.html`. A failed run retains an
`INCOMPLETE.txt` marker and logs; the presence of some result files does not mean
all requested analyses finished.

## 5. Analyze an existing gene-count table

Use this route when raw gene counts are already available. No BAM or GTF is
needed for this DGE-only analysis. Counts must be nonnegative integers—not TPM,
FPKM, normalized counts or transformed values.

The table's first column contains unique gene IDs; subsequent columns identify
samples. This two-row illustration shows the format, not a sufficient dataset
for differential-expression analysis:

```tsv
gene_id	WT1	WT2	WT3	KO1	KO2	KO3
GENE_A	120	135	118	260	245	271
GENE_B	80	91	85	78	87	82
```

Save your **complete** count table as `analysis/gene_counts.tsv`. Create
`analysis/count_samples.tsv` with the same sample IDs and their conditions:

```tsv
sample_id	condition
WT1	WT
WT2	WT
WT3	WT
KO1	KO
KO2	KO
KO3	KO
```

Use `nano analysis/count_analysis.json` to create this configuration:

```json
{
  "samples": "count_samples.tsv",
  "gene_counts": "gene_counts.tsv",
  "out": "../results/gene_expression",
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

Run:

```bash
python3 TranscriptomePipeline.py --config analysis/count_analysis.json --check
python3 TranscriptomePipeline.py --config analysis/count_analysis.json
```

Standard featureCounts annotation columns are supported. If its count-column
names are BAM paths rather than sample IDs, add a `count_column` column to the
sample sheet containing each exact count-column name. See
[count-table requirements](docs/INTEGRATED_WORKFLOW.md#start-with-the-count-only-example).

## 6. Run standalone DEJU

`DEJUPipeline.R` is a separate command-line entry point for a two-condition exon/
junction-usage analysis. It does not run DESeq2. Keep `DEJUHelpers.R` beside it.

Try the included feature-count example:

```bash
Rscript DEJUPipeline.R \
  --samples examples/counts/samples.tsv \
  --counts examples/counts/counts.tsv \
  --features examples/counts/features.tsv \
  --reference WT \
  --treatment KO \
  --covariates batch \
  --out results/example_deju
```

A backslash `\` continues the same shell command on the next line. It must be the
last character on the line, with no spaces after it. Alternatively, put the
entire command on one line without backslashes.

To use the real BAM sample sheet from section 4:

```bash
Rscript DEJUPipeline.R \
  --samples analysis/samples.tsv \
  --gtf /absolute/path/reference/annotation.gtf \
  --paired true \
  --strand 2 \
  --reference WT \
  --treatment KO \
  --threads 4 \
  --out results/standalone_deju
```

Replace the GTF path and sequencing settings. Add `--covariates batch` only if
that column exists and belongs in the model. Standalone covariates are categorical;
use the integrated entry point for continuous covariates or interaction contrasts.

```bash
Rscript DEJUPipeline.R --help
```

For exact counting rules, thresholds, count-table formats, GO options and all
output columns, see the [standalone DEJU reference](docs/DEJU.md).

## 7. Analyze isoform switching with Salmon

`IsoformSwitchAnalysis.R` starts from existing Salmon `quant.sf` files. It tests
whether transcripts change their fraction of a gene's expression. It does not
run Salmon or infer transcript counts from a gene-count matrix.

Install the additional packages, then try the synthetic example:

```bash
Rscript install_dependencies.R --isoform-switch
Rscript IsoformSwitchAnalysis.R examples/isoform_switch/config.json --check
Rscript IsoformSwitchAnalysis.R examples/isoform_switch/config.json
```

Outputs go to `results/example_isoform_switch`. This example disables optional
consequence prediction and switch-gene plots. The module requires the current
count-based IsoformSwitchAnalyzeR filter API; an incompatible older installation
is reported during preflight.

For real data, create `analysis/isoform_samples.tsv`:

```tsv
sample_id	condition	quant
WT1	WT	/absolute/path/salmon/WT1/quant.sf
WT2	WT	/absolute/path/salmon/WT2/quant.sf
WT3	WT	/absolute/path/salmon/WT3/quant.sf
KO1	KO	/absolute/path/salmon/KO1/quant.sf
KO2	KO	/absolute/path/salmon/KO2/quant.sf
KO3	KO	/absolute/path/salmon/KO3/quant.sf
```

Create `analysis/isoform_analysis.json`:

```json
{
  "samples": "isoform_samples.tsv",
  "gtf": "/absolute/path/reference/annotation.gtf",
  "out": "../results/isoform_switch",
  "comparisons": [
    {"name": "KO_vs_WT", "reference": "WT", "treatment": "KO"}
  ],
  "covariates": {},
  "alpha": 0.05,
  "delta_if": 0.1,
  "strip_pipe": false,
  "consequences": false,
  "predict_novel_orfs": false,
  "plots": 10
}
```

Run:

```bash
Rscript IsoformSwitchAnalysis.R analysis/isoform_analysis.json --check
Rscript IsoformSwitchAnalysis.R analysis/isoform_analysis.json
```

All quantifications must use a consistent transcript reference matching the GTF.
`alpha` controls the significance threshold; `delta_if` sets the minimum absolute
change in isoform fraction. Positive dIF means greater usage in KO than WT.
`strip_pipe` can be enabled for pipe-delimited transcript identifiers; version
suffixes are preserved and ambiguous ID collisions are rejected.

For annotated coding/NMD consequences, add a matching **transcript FASTA** and
configure the optional stages as described in the [isoform-switch guide](docs/ISOFORM_SWITCH.md).
That guide also covers batch/subject covariates, GO enrichment, sequence exports,
and calling the module from `TranscriptomePipeline.py`.

## 8. Add other analyses

The integrated workflow selects modules through its JSON file. Adding a module
also requires its inputs, settings and software; simply adding its name is not
enough. These guides provide the additional configuration and interpretation:

| Analysis | Setup guide |
|---|---|
| Event classification, PSI and ΔPSI | [rMATS event analysis](docs/INTEGRATED_WORKFLOW.md#event-classification-and-psi) |
| Transcript assembly and novel-locus candidates | [StringTie/GffCompare analysis](docs/INTEGRATED_WORKFLOW.md#transcript-isoforms-and-novel-loci) |
| Fusion candidates | [Arriba inputs and settings](docs/INTEGRATED_WORKFLOW.md#fusion-candidates) |
| Coverage and junction-arc plots | [Locus coordinates and plotting](docs/INTEGRATED_WORKFLOW.md#locus-visualizations) |
| Multiple conditions, interactions and continuous covariates | [Designs and contrasts](docs/INTEGRATED_WORKFLOW.md#model-specification) |

After saving an expanded configuration, use the same check/run commands from
section 4 with that file's path. General DGE/DEJU designs are not automatically
applied to the separate rMATS or isoform-switch models. Novel loci and fusions
are candidates for further investigation, not confirmed discoveries by themselves.

To obtain reference files, list the supplied GENCODE v48 resources:

```bash
python3 references/DownloadReferences.py --list
```

Read the [reference guide](references/README.md) before downloading. It explains
CHR/PRI/ALL annotations, genome versus transcript FASTA, checksums and matching
the reference already used for alignment or quantification.

## 9. Find and interpret results

For an integrated run, begin with `report.html`, then inspect QC and the module
folders. The following paths are relative to the configured output directory:

| Output | What to inspect |
|---|---|
| `report.html` | Run overview, module status and links to tables/PDFs |
| `COMPLETE.txt` / `INCOMPLETE.txt` | Whether every requested stage finished |
| `logs/` | Commands and execution messages; start here after an error |
| `dge/` | Gene statistics, normalized counts, PCA and diagnostic plots |
| `dju/COMPARISON_NAME/` | Junction/exon usage statistics and diagnostics |
| `native/` | Design and, for BAM analyses, alignment/assignment QC and gene counts |
| `events/` | Splicing-event results, replicate PSI and ΔPSI |
| `isoforms/` | Transcript assemblies, abundance and structural classifications |
| `isoform_switch/` | Isoform-usage statistics and requested optional outputs |
| `loci/` | Coverage and junction-arc PDFs |
| `fusions/` | Fusion candidates and caller evidence |

Standalone DEJU places files such as `junctions_results.tsv.gz` directly in its
output directory. Standalone isoform switching similarly writes files such as
`all_tested_isoforms.tsv` directly in its own output directory. Their layouts
are documented in their respective guides.

Read QC before interpreting significance. FDR adjustments belong to the stated
analysis family/comparison; testing several modules does not create one combined
FDR guarantee. Predicted ORF/NMD consequences and reconstructed structures require
appropriate biological supporting evidence.

## 10. Troubleshooting and further documentation

| Message or symptom | What to check |
|---|---|
| `python3`, `Rscript` or another command is missing | The required program is installed and available in the current terminal environment |
| `can't open file` or missing script | Run `pwd` and `ls`; enter the repository directory and use the exact filename/capitalization |
| Missing input file | Replace placeholder paths; remember that JSON paths are relative to the JSON file and BAM/quant paths to the sample sheet |
| Invalid JSON | Use double quotes, balanced brackets and no trailing commas; run `python3 -m json.tool analysis/bam_analysis.json` to check syntax |
| Output directory exists | Choose a new `out` path; do not reuse a partial run's folder |
| Count columns do not match samples | Match exact sample IDs or provide `count_column` mappings |
| Confounded/rank-deficient design | Check condition/batch/subject assignments and replication; effects cannot be separated when the design confounds them |
| Very few BAM alignments retained | Inspect NH tags, mapping quality, alignment flags and the alignment-QC table |
| Missing/incompatible R package | Use a coherent R/Bioconductor installation and rerun the appropriate dependency command |

**Validation status:** native DESeq2/DEJU analyses have synthetic count and BAM
test coverage. External caller execution and the full IsoformSwitchAnalyzeR
workflow still require end-to-end validation in a suitably provisioned environment.
See [validation coverage and limitations](docs/VALIDATION_V0.2.md) for the exact scope.

* [Detailed integrated configuration](docs/INTEGRATED_WORKFLOW.md)
* [Standalone DEJU reference](docs/DEJU.md)
* [Isoform-switch methods and outputs](docs/ISOFORM_SWITCH.md)
* [Reference downloads and third-party reuse](references/README.md)
* [Tests and synthetic examples](tests/README.md)
* [Version history](CHANGES.md)

Please cite the analysis methods used in a study. Method and software links are
provided in the detailed guides and [reference resources](references/README.md).
