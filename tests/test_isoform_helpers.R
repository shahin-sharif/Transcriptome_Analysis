#!/usr/bin/env Rscript
args <- commandArgs(FALSE);script <- sub('^--file=','',grep('^--file=',args,value=TRUE)[1])
root <- dirname(dirname(normalizePath(script)))
source(file.path(root,'IsoformSwitchHelpers.R'))
source(file.path(root,'IsoformSwitchAnalysis.R'))
fails <- function(expr) inherits(tryCatch({force(expr);NULL},error=identity),'error')
cfg <- is_config(file.path(root,'examples/isoform_switch/config.json'))
# A release-style import function deliberately lacks the development-only option.
import_args <- is_import_args(cfg)
stopifnot(!'autoCastDesignCol' %in% names(import_args))
release_import <- function() NULL
formals(release_import) <- as.pairlist(setNames(rep(list(NULL),length(import_args)),names(import_args)))
stopifnot(is_check_api('importRdata',names(import_args),release_import))
stopifnot(fails(is_check_api('importRdata',c(names(import_args),'autoCastDesignCol'),release_import)))
d <- is_design(cfg)
stopifnot(nrow(d$meta)==6,d$comparisons$condition_1=='control',d$comparisons$condition_2=='treated')
stopifnot(fails(is_id(c('a|one','a|two'),TRUE)))
txi <- is_quant(d$meta,FALSE)
a <- is_annotation(cfg,rownames(txi$counts))
stopifnot(nrow(a)==120,all(colnames(txi$counts)==d$meta$sample_id),all(abs(colSums(txi$abundance)-1e6)<.01))
bad <- cfg;bad$comparisons[[1]]$reference <- 'absent';stopifnot(fails(is_design(bad)))
confounded <- d$meta;confounded$batch <- confounded$condition
sheet <- tempfile(fileext='.tsv');is_write(confounded,sheet)
bad <- cfg;bad$samples <- sheet;bad$covariates <- list(batch='factor')
stopifnot(fails(is_design(bad)));unlink(sheet)
stopifnot(fails(is_annotation(cfg,c(rownames(txi$counts),'wrong_transcript'))))
tab <- data.frame(isoform_id=c('a','b','a','b'),gene_id='g',condition_1=c('c','c','c','c'),condition_2=c('t','t','u','u'),
  IF1=.5,IF2=c(.8,.2,.55,.45),dIF=c(.3,-.3,.05,-.05),isoform_switch_q_value=c(.01,.01,.01,.01),gene_switch_q_value=c(.02,.02,.03,.03))
comp <- data.frame(name=c('t_vs_c','u_vs_c'),condition_1='c',condition_2=c('t','u'))
r <- is_results(tab,comp,.05,.1)
stopifnot(nrow(r$all)==4,nrow(r$significant)==2,nrow(r$genes)==2,all(r$significant$comparison=='t_vs_c'))
wrong <- tab;wrong$dIF <- -wrong$dIF;stopifnot(fails(is_results(wrong,comp,.05,.1)))
# Real tximport/GTF/FASTA validation above. Optional full package test is explicit.
if('--full' %in% commandArgs(TRUE)) {
  cfg$out <- tempfile('isoform-switch-test-')
  r <- is_run(cfg)
  planted <- sprintf('SYN_G%03d',1:5)
  stopifnot(all(planted %in% r$significant$gene_id))
  t <- r$all[r$all$gene_id %in% planted,,drop=FALSE]
  stopifnot(all(t$dIF[grepl('_A',t$isoform_id)]<0),all(t$dIF[grepl('_B',t$isoform_id)]>0))
  cat('Full IsoformSwitchAnalyzeR/DEXSeq fixture completed:',cfg$out,'\n')
}
cat('Isoform helper tests passed (full package analysis only with --full).\n')
