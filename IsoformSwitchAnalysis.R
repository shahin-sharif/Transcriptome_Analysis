#!/usr/bin/env Rscript
# Standalone and integrated Salmon isoform-switch analysis.
is_run <- function(cfg,check=FALSE) {
  d <- is_design(cfg)
  pkgs <- is_dependencies(cfg)
  if(check) {message('Isoform-switch configuration, design and dependencies passed. Quantification/annotation content is checked when running.');return(invisible(NULL))}
  is_assert(!file.exists(cfg$out),'Output exists; choose a new isoform-switch output directory')
  dir.create(cfg$out,recursive=TRUE)
  out <- normalizePath(cfg$out,mustWork=TRUE)
  writeLines('Analysis incomplete. Inspect log and saved core results.',file.path(out,'INCOMPLETE.txt'))
  on.exit(writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt')),add=TRUE)
  jsonlite::write_json(cfg,file.path(out,'resolved_config.json'),auto_unbox=TRUE,pretty=TRUE)
  is_write(data.frame(package=pkgs,version=vapply(pkgs,function(p) as.character(utils::packageVersion(p)),character(1))),file.path(out,'package_versions.tsv'))
  is_write(d$meta,file.path(out,'samples.tsv'));is_write(d$design,file.path(out,'design.tsv'));is_write(d$comparisons,file.path(out,'comparisons.tsv'))
  inputs <- unique(c(cfg$samples,cfg$gtf,cfg$transcript_fasta,d$meta$quant,if(!is.null(cfg$go)) cfg$go$mapping))
  is_write(data.frame(path=inputs,md5=unname(tools::md5sum(inputs))),file.path(out,'input_checksums.tsv'))
  txi <- is_quant(d$meta,cfg$strip_pipe)
  annot <- is_annotation(cfg,rownames(txi$counts))
  is_write(annot,file.path(out,'transcript_to_gene.tsv'))
  matrix_table <- function(x) data.frame(isoform_id=rownames(x),x,check.names=FALSE)
  for(k in c('counts','abundance','length')) is_write(matrix_table(txi[[k]]),file.path(out,paste0('salmon_',k,'.tsv')))
  is_write(data.frame(sample_id=d$meta$sample_id,scaled_count_total=colSums(txi$counts),
    tpm_total=colSums(txi$abundance),expressed_transcripts=colSums(txi$abundance>0)),file.path(out,'sample_qc.tsv'))
  logtpm <- log2(txi$abundance+1)
  variable <- apply(logtpm,1,stats::sd)>0
  qc <- file.path(out,'sample_qc.pdf');grDevices::pdf(qc,width=9,height=7)
  tryCatch({
    boxplot(logtpm,names=d$meta$sample_id,las=2,ylab='log2(TPM + 1)',main='Transcript abundance (unadjusted)')
    if(sum(variable)>=2) {
      pc <- prcomp(t(logtpm[variable,,drop=FALSE]),center=TRUE,scale.=FALSE)
      pv <- 100*pc$sdev^2/sum(pc$sdev^2)
      plot(pc$x[,1:2],pch=19,col=as.integer(d$design$condition),xlab=sprintf('PC1 (%.1f%%)',pv[1]),ylab=sprintf('PC2 (%.1f%%)',pv[2]),main='Transcript TPM PCA')
      text(pc$x[,1:2],labels=d$meta$sample_id,pos=3,cex=.7)
      legend('topright',legend=levels(d$design$condition),col=seq_along(levels(d$design$condition)),pch=19)
      is_write(data.frame(sample_id=d$meta$sample_id,pc$x),file.path(out,'pca.tsv'))
    }
  },finally=grDevices::dev.off())
  sw <- is_api('importRdata',is_import_args(cfg,d$design,
    d$comparisons[,c('condition_1','condition_2')],matrix_table(txi$counts),matrix_table(txi$abundance)))
  imported <- unique(sw$isoformFeatures$isoform_id)
  f <- cfg$filter
  sw <- is_api('preFilter',list(switchAnalyzeRlist=sw,isoCount=f$iso_count,min.Count.prop=f$count_proportion,
    IFcutoff=f$if_cutoff,min.IF.prop=f$if_proportion,removeSingleIsoformGenes=TRUE,
    reduceToSwitchingGenes=FALSE,reduceFurtherToGenesWithConsequencePotential=FALSE,
    keepIsoformInAllConditions=TRUE,alpha=cfg$alpha,dIFcutoff=cfg$delta_if))
  retained <- unique(sw$isoformFeatures$isoform_id)
  is_write(data.frame(isoform_id=rownames(txi$counts),imported=rownames(txi$counts) %in% imported,
    retained_for_testing=rownames(txi$counts) %in% retained),file.path(out,'filter_audit.tsv'))
  is_assert(length(retained)>1,'No testable multi-isoform genes after filtering')
  sw <- is_api('isoformSwitchTestDEXSeq',list(switchAnalyzeRlist=sw,alpha=cfg$alpha,dIFcutoff=cfg$delta_if,
    reduceToSwitchingGenes=FALSE,reduceFurtherToGenesWithConsequencePotential=FALSE,keepIsoformInAllConditions=TRUE))
  saveRDS(sw,file.path(out,'all_tested_switches.rds'))
  if(is.data.frame(sw$isoformSwitchAnalysis)) is_write(sw$isoformSwitchAnalysis,file.path(out,'isoform_switch_statistics.tsv'))
  results <- is_results(sw$isoformFeatures,d$comparisons,cfg$alpha,cfg$delta_if)
  is_write(results$all,file.path(out,'all_tested_isoforms.tsv'))
  is_write(results$significant,file.path(out,'significant_isoforms.tsv'))
  is_write(results$genes,file.path(out,'gene_switch_summary.tsv'))
  grDevices::pdf(file.path(out,'switch_volcano.pdf'),width=8,height=7)
  tryCatch({for(cmp in d$comparisons$name) {
    t <- results$all[results$all$comparison==cmp,,drop=FALSE]
    y <- -log10(pmax(t$isoform_switch_q_value,.Machine$double.xmin));ok <- is.finite(t$dIF) & is.finite(y)
    if(!any(ok)) {plot.new();title(paste(cmp,'— no testable isoforms'));next}
    plot(t$dIF[ok],y[ok],pch=16,cex=.6,col=ifelse(t$significant[ok],'firebrick','grey50'),
      xlab='dIF (treatment minus reference)',ylab='-log10 isoform q-value',main=cmp)
    abline(v=c(-cfg$delta_if,cfg$delta_if),h=-log10(cfg$alpha),lty=2)
  }},finally=grDevices::dev.off())
  if(!is.null(cfg$go)) is_go(results$all,cfg,out)
  # Core results are durable before expensive optional annotation and plots.
  has_switches <- nrow(results$significant)>0
  if(cfg$consequences && has_switches) {
    genes <- unique(results$significant$gene_id)
    swc <- is_api('subsetSwitchAnalyzeRlist',list(switchAnalyzeRlist=sw,subset=sw$isoformFeatures$gene_id %in% genes))
    swc <- is_api('analyzeAlternativeSplicing',list(switchAnalyzeRlist=swc,onlySwitchingGenes=FALSE,alpha=cfg$alpha,dIFcutoff=cfg$delta_if))
    if(cfg$predict_novel_orfs && any(swc$orfAnalysis$orf_origin=='not_annotated_yet',na.rm=TRUE)) {
      swc <- is_api('analyzeNovelIsoformORF',list(switchAnalyzeRlist=swc,analysisAllIsoformsWithoutORF=FALSE))
    }
    swc <- is_api('extractSequence',list(switchAnalyzeRlist=swc,onlySwitchingGenes=TRUE,alpha=cfg$alpha,dIFcutoff=cfg$delta_if,
      extractNTseq=TRUE,extractAAseq=TRUE,addToSwitchAnalyzeRlist=TRUE,writeToFile=TRUE,pathToOutput=out,outputPrefix='switching_isoforms'))
    swc <- is_api('analyzeSwitchConsequences',list(switchAnalyzeRlist=swc,consequencesToAnalyze=c('NMD_status','ORF_seq_similarity','intron_retention'),
      alpha=cfg$alpha,dIFcutoff=cfg$delta_if,removeNonConseqSwitches=FALSE))
    saveRDS(swc,file.path(out,'switches_with_predicted_consequences.rds'))
    for(k in c('orfAnalysis','AlternativeSplicingAnalysis','switchConsequence')) {
      x <- swc[[k]]
      if(is.data.frame(x)) is_write(x,file.path(out,paste0(k,'.tsv')))
    }
    sw <- swc
  }
  if(cfg$plots>0 && has_switches) {
    is_plot_switches(sw,cfg,file.path(out,'switch_plots'))
  }
  is_write(data.frame(stage=c('core_dexseq','consequences','switch_plots','go'),status=c('completed',
    if(!cfg$consequences) 'not_requested' else if(!has_switches) 'no_significant_switches' else 'completed',
    if(cfg$plots==0) 'not_requested' else if(!has_switches) 'no_significant_switches' else 'completed',
    if(is.null(cfg$go)) 'not_requested' else 'completed')),file.path(out,'stage_status.tsv'))
  writeLines('Completed requested isoform-switch stages. Review QC and predicted-consequence limitations.',file.path(out,'COMPLETE.txt'))
  unlink(file.path(out,'INCOMPLETE.txt'))
  invisible(results)
}

is_main <- function() {
  args <- commandArgs(TRUE)
  if(!length(args) || length(args)>2 || (length(args)==2 && args[2]!='--check'))
    stop('Usage: Rscript IsoformSwitchAnalysis.R config.json [--check]',call.=FALSE)
  is_run(is_config(args[1]),check='--check' %in% args)
}
if(sys.nframe()==0) {
  script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
  source(file.path(dirname(normalizePath(script)),'IsoformSwitchHelpers.R'))
  is_main()
}
