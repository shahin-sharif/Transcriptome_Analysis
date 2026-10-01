#!/usr/bin/env Rscript
arg <- grep('^--file=',commandArgs(FALSE),value=TRUE)
root <- dirname(dirname(normalizePath(sub('^--file=','',arg[1]))))
source(file.path(root,'IsoformSwitchHelpers.R'))
fails <- function(expr) inherits(tryCatch({force(expr);NULL},error=identity),'error')
ids <- sprintf('SYN_TX%02d',1:14)
tab <- data.frame(isoform_id=rep(ids,2),gene_id='SYN_GENE',
  condition_1='WT',condition_2=rep(c('KO1','KO2'),each=14),
  IF1=rep(c(.35,.35,rep(.025,12)),2),IF2=c(.55,.15,rep(.025,12),.35,.35,rep(.025,12)),
  isoform_switch_q_value=c(.001,.001,rep(.9,12),rep(.9,14)),gene_switch_q_value=.001)
tab$dIF <- tab$IF2-tab$IF1
ex <- GenomicRanges::GRanges('chr1',IRanges::IRanges(rep(c(101,301),14),width=50),strand='-')
ex$isoform_id <- rep(ids,each=2)
sw <- list(isoformFeatures=tab,exons=ex,orfAnalysis=NULL)
cfg <- list(comparisons=list(list(name='KO1_vs_WT',reference='WT',treatment='KO1'),
                             list(name='KO2_vs_WT',reference='WT',treatment='KO2')),
            alpha=.05,delta_if=.1,plots=10)
snapshot <- serialize(sw,NULL)
tmp <- tempfile('isoform-plot-recovery-');dir.create(tmp)
input <- file.path(tmp,'saved');dir.create(input)
saveRDS(sw,file.path(input,'all_tested_switches.rds'))
jsonlite::write_json(cfg,file.path(input,'resolved_config.json'),auto_unbox=TRUE)
output <- file.path(tmp,'plots')
index <- is_plot_switches(sw,cfg,output)
stopifnot(nrow(index)==1,index$comparison=='KO1_vs_WT',index$pages==2,
          file.info(file.path(output,index$file))$size>1000,
          identical(snapshot,serialize(sw,NULL)),file.exists(file.path(output,'PLOTS_COMPLETE.txt')),
          !file.exists(file.path(output,'INCOMPLETE.txt')),fails(is_plot_switches(sw,cfg,output)))
plotted <- read.delim(file.path(output,'plotted_isoform_statistics.tsv'))
stopifnot(nrow(plotted)==14,all(plotted$comparison=='KO1_vs_WT'),
          all(abs(plotted$dIF-(plotted$IF2-plotted$IF1))<1e-8))
# End-to-end CLI reads saved statistics only; fake input has no BAMs, quant or GTF paths.
cli_out <- file.path(tmp,'cli')
status <- system2(file.path(R.home('bin'),'Rscript'),
  c('--vanilla',shQuote(file.path(root,'PlotIsoformSwitches.R')),shQuote(input),shQuote(cli_out)),
  stdout=file.path(tmp,'cli.log'),stderr=file.path(tmp,'cli.log'))
if(status!=0) stop(paste(readLines(file.path(tmp,'cli.log')),collapse='\n'))
stopifnot(file.exists(file.path(cli_out,'source_checksums.tsv')),
          file.exists(file.path(cli_out,'PLOTS_COMPLETE.txt')))
zero <- cfg;zero$plots<-0
stopifnot(nrow(is_plot_core_switches(sw,zero,file.path(tmp,'zero')))==0)
cat('PASS: missing ORF fallback, per-comparison selection, multipage plots, unchanged model, overwrite protection and saved-result CLI\n')
cat('Plot fixture:',file.path(output,index$file),'\n')

if(nzchar(Sys.getenv('PLOT_TEST_PDF'))) file.copy(file.path(output,index$file),Sys.getenv('PLOT_TEST_PDF'),overwrite=TRUE)
