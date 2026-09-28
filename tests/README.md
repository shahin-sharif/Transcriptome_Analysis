# Validation

Run from the repository root after installing dependencies:

```bash
Rscript tests/test_pipeline.R
Rscript tests/test_pipeline.R --bam
Rscript tests/test_pipeline.R --bam --go
```

The test runner creates deterministic synthetic data in a temporary directory.
`--bam` adds six-sample single-end and paired-end BAM generation/counting/fitting
and a full paired-end CLI run. No patient data or external download is needed.
`--go` additionally requires clusterProfiler and org.Hs.eg.db. Its arbitrary
synthetic gene mapping checks execution only, not biological enrichment.

Covered behavior:

* Modern and legacy statistical engines; expected signs and detection of planted
  usage changes; an expression-only control; exact junction BH family.
* Sample row/column order invariance; mismatched sample IDs, negative counts and
  confounded designs are rejected.
* Sample-column matching preserves the intended sample order.
* Minus-strand genomic junction coordinates and transcript ordering; annotated,
  novel, unassigned and ambiguous gene assignments.
* CIGARs containing multiple splice gaps, soft clipping and deletions.
* Primary/unique/MAPQ/QC/duplicate filtering before exon and junction counting.
* Both mates crossing the same junction count twice under the declared read policy.
* Forward/reverse orientation, chunk boundaries and pairedness validation.
* Full CLI outputs, plots, completion markers, transcript-position export,
  cleanup of filtered BAMs and refusal to overwrite existing outputs.
* Optional GO execution and explicit mapped-background exports.

The shipped count fixture can be regenerated with:

```r
source('tests/make_fixtures.R')
make_count_fixture('examples/counts')
```

Initial local validation environment: R 4.2.1 on Intel macOS, edgeR 4.10.5,
limma 3.68.5, Rsubread 2.10.5, Rsamtools 2.12.0, GenomicAlignments 1.32.1,
GenomicRanges 1.48.0. edgeR/limma were built in an isolated library; this mixed-age
local test environment is not an installation recommendation. Use a current,
coherent R/Bioconductor environment for new analyses. The runtime saves
`sessionInfo.txt` for each successful analysis.

A separate local check exercised optional GO with the installed organism database.
No full experimental RNA-seq dataset has been validated with this release.
Passing fixtures demonstrate tested software behavior, not absence of every bug
or biological correctness of an unspecified experiment.

## Integrated v0.2 tests

See [v0.2 validation](../docs/VALIDATION_V0.2.md). Run `python3 tests/test_integrated.py` for adapter/report contracts and `Rscript tests/test_native_v2.R --bam` for the integrated native statistics/counting/plotting checks. External caller tests are explicitly distinct from actual caller execution.

## Isoform and reference tests

```bash
Rscript tests/test_isoform_helpers.R
python3 tests/test_references.py
# Requires a current compatible IsoformSwitchAnalyzeR installation:
Rscript tests/test_isoform_helpers.R --full
```

The default isoform tests execute tximport and annotation checks, not the absent
IsoformSwitchAnalyzeR package. The optional full test asserts five planted gene
switches and their directions. Regenerate the deterministic example inputs with
`python3 tests/make_isoform_fixture.py`. They are synthetic Salmon-shaped data.
Reference tests use an in-memory stream, never download a human genome.
