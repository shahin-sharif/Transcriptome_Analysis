# Validation status — v0.2.0

This revision has not been validated against the user's real experimental count
matrix or original isoform-switch results. Implementation, synthetic test
coverage and full biological validation are distinct, as detailed below.

## Executed native analyses

* DESeq2 on 160 synthetic genes and 12 samples with a batch effect in the design,
  a continuous time variable, a two-condition comparison and an interaction.
  Baseline, later-time and interaction contrasts completed. The test asserted
  significance and direction for the five planted interaction genes.
* DEJU on independently simulated exon/junction counts with that interaction
  contrast; the planted junction's direction and significance passed.
* Three-condition model/contrast construction, continuous-variable handling,
  rank-deficiency/confounding rejection and exact coefficient mapping passed.
* Full native BAM workflow on six synthetic paired-end libraries: filtering/QC,
  independent gene-level fragment counts, DESeq2, exon/junction read counts,
  DEJU, transcript positions, model exports and locus PDF.
* Gene-count totals in the unambiguous synthetic fixture equal half of the
  retained alignment count, confirming fragments rather than both mates are
  counted in the gene matrix. This equality is fixture-specific, not expected
  for real BAMs with unassigned/ambiguous reads.
* Original standalone DEJU regression tests passed after adding anchor/intron
  filters and the optional general-design interface.
* Locus and DESeq2 PCA PDFs were rendered and inspected. They are raw-coverage
  and exploratory QC views, not independent statistical evidence.

One deliberately under-dispersed, nearly deterministic BAM fixture caused DESeq2
to reject dispersion fitting. The workflow correctly retained an incomplete
status and the error log. The statistical test fixture was then generated with
biological count variation. No silent dispersion-model fallback was introduced.

## External adapters: different validation level

Nine offline Python tests exercise five event types, effective-length PSI and
its direction, joint event-type BH, event/junction coordinate links, unmatched
junctions, paired-subject ordering, malformed-input rejection, transcript candidate
labels, missing quantification handling and fusion schemas. The rMATS adapter's
actual argument assembly/group ordering and parser, plus isoform-switch configuration
path resolution and shared-sample identity, are tested with a simulated
executable boundary.

**rMATS-turbo, StringTie, GffCompare and Arriba have not been run on biological
BAMs in this local environment; those executables are not installed here.**
The adapters are implemented, but caller execution, dependency compatibility and
biological performance require an appropriately provisioned analysis environment
and suitable inputs. Passing parser/adapter tests must not be described as
validation of those callers or their biological discoveries.

## Isoform-switch module

The standalone/integrated Salmon isoform-switch module is implemented. Real
`tximport` import and GTF/FASTA matching passed on six synthetic quantifications
(120 transcripts). Tests also cover ID collisions, unknown conditions, comparison
separation, effect/FDR filters and dIF sign. The module's named package arguments
were checked against upstream source. Full DEXSeq switch fitting, consequence
prediction, sequence export and switch plots through IsoformSwitchAnalyzeR have
**not** been executed here because that package is absent. Use
`Rscript tests/test_isoform_helpers.R --full` in a compatible current environment;
it asserts recovery and direction of five planted switching genes.

## Reference downloads

The release-48 file names and MD5 checksums were retrieved from GENCODE's
upstream release manifest. Offline tests verify successful download integrity,
reuse of a valid file, refusal to overwrite a mismatched file, and cleanup after
checksum failure. Full human references were not downloaded during these tests.

## Reproduce tests

```bash
python3 tests/test_integrated.py
Rscript tests/test_native_v2.R
Rscript tests/test_native_v2.R --bam
Rscript tests/test_pipeline.R
Rscript tests/test_isoform_helpers.R
python3 tests/test_references.py
```

The native BAM test generates small alignments locally; no real/sample-identifying
data is shipped. Local R testing used R 4.2.1, DESeq2 1.36.0, edgeR 4.10.5,
limma 3.68.5, Rsubread 2.10.5, Rsamtools 2.12.0, GenomicAlignments 1.32.1 and
GenomicRanges 1.48.0. edgeR/limma used an isolated library. Use a coherent current
R/Bioconductor installation for production, not this mixed-age local setup.

## Before interpreting real discoveries

Check sample labels and independent biological replicates; aligner/NH tagging;
read length, library orientation and assembly consistency; batch/subject/time
model definitions; count units; outlier/assignment reports; and representative
locus evidence. Predicted structures, novel loci, NMD effects and fusions require
appropriate supporting evidence and are not established solely by this pipeline.
