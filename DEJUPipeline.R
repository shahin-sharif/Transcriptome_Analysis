#!/usr/bin/env Rscript
# Differential exon/junction usage: joint exon + junction negative-binomial model.
DEJU_VERSION <- '0.2.0'
`%||DEJU%` <- function(x,y) if (is.null(x)) y else x

# Locate the helper beside this file both for Rscript and source().
.deju_file <- tryCatch(sys.frame(1)$ofile, error=function(e) NULL)
if (is.null(.deju_file)) {
  .deju_arg <- grep('^--file=', commandArgs(FALSE), value=TRUE)
  .deju_file <- if (length(.deju_arg)) sub('^--file=', '', .deju_arg[1]) else 'DEJUPipeline.R'
}
source(file.path(dirname(normalizePath(.deju_file, mustWork=TRUE)), 'DEJUHelpers.R'))

# Aligned reads are the counting unit in both exon and junction channels.
# Pre-filtering also applies to junctions, which featureCounts may collect independently
# of exon assignment. Unique mapping still requires the aligner's NH tag.
deju_count_bams <- function(samples, annotation, outdir, paired, strand, threads=1L,
                            min_mapq=10L, remove_duplicates=FALSE, pre_filtered=FALSE,
                            min_anchor=8L, min_intron=20L, max_intron=1000000L) {
  deju_require(c('Rsubread','Rsamtools'))
  scratch <- file.path(outdir, 'counting_tmp'); dir.create(scratch)
  on.exit(unlink(scratch, recursive=TRUE), add=TRUE)
  filtered <- file.path(scratch, paste0(samples$sample_id,'.bam'))
  flags <- Rsamtools::scanBamFlag(isUnmappedQuery=FALSE, isSecondaryAlignment=FALSE,
                                 isSupplementaryAlignment=FALSE, isNotPassingQualityControls=FALSE,
                                 isDuplicate=if (remove_duplicates) FALSE else NA)
  if (pre_filtered) filtered <- samples$bam
  for (i in if (pre_filtered) integer() else seq_len(nrow(samples))) {
    # No genomic insert-size filter: it would reject ordinary spliced RNA pairs.
    Rsamtools::filterBam(samples$bam[i], filtered[i], index=character(),
                        param=Rsamtools::ScanBamParam(flag=flags, mapqFilter=min_mapq, tagFilter=list(NH=1L)))
  }
  fc <- Rsubread::featureCounts(files=filtered, annot.ext=annotation$saf,
    isGTFAnnotationFile=FALSE, useMetaFeatures=FALSE, strandSpecific=strand,
    isPairedEnd=paired, countReadPairs=FALSE, requireBothEndsMapped=FALSE,
    checkFragLength=FALSE, countMultiMappingReads=FALSE, primaryOnly=FALSE,
    allowMultiOverlap=FALSE, nonSplitOnly=TRUE, juncCounts=FALSE, nthreads=threads)
  ex <- fc$annotation
  ecounts <- fc$counts[,deju_sample_columns(colnames(fc$counts),filtered,samples$sample_id),drop=FALSE]
  ef <- data.frame(feature_id=paste0('E:',ex$GeneID,':',ex$Chr,':',ex$Start,':',ex$End,':',ex$Strand),
                   gene_id=as.character(ex$GeneID), feature_type='exon', chr=as.character(ex$Chr),
                   left=ex$Start, right=ex$End, strand=ex$Strand, stringsAsFactors=FALSE)
  junction_data <- deju_junction_counts(filtered,samples$sample_id,paired,strand,min_anchor=min_anchor,min_intron=min_intron,max_intron=max_intron)
  jc <- junction_data$counts
  ja <- deju_assign_junctions(junction_data$annotation,annotation)
  deju_write(cbind(ja,setNames(as.data.frame(jc),samples$sample_id)), file.path(outdir,'junction_assignment.tsv.gz'))
  keep <- !is.na(ja$gene_id)
  deju_assert(any(keep), 'No unambiguous intragenic junctions remain')
  ja <- ja[keep,,drop=FALSE]; jc <- jc[keep,,drop=FALSE]
  jf <- data.frame(feature_id=paste0('J:',ja$gene_id,':',ja$chr,':',ja$left,':',ja$right,':',ja$strand),
                   gene_id=ja$gene_id, feature_type='junction', chr=ja$chr,
                   left=ja$left, right=ja$right, strand=ja$strand, stringsAsFactors=FALSE)
  deju_assert(!anyDuplicated(jf$feature_id), 'Duplicate assigned junction coordinates')
  f <- rbind(ef,jf); counts <- rbind(ecounts,jc)
  rownames(counts) <- f$feature_id; colnames(counts) <- samples$sample_id
  stat <- fc$stat; names(stat)[-1] <- samples$sample_id
  deju_write(stat,file.path(outdir,'featurecounts_assignment.tsv'))
  list(counts=counts,features=f)
}

deju_fit <- function(counts, features, samples, reference, treatment, covariates=character(),
                     engine='modern', min_count=10, min_total_count=15, fdr=0.05, min_abs_log2=0, design_override=NULL) {
  deju_require(c('edgeR','limma','statmod'))
  deju_assert(engine %in% c('modern','legacy'), 'engine must be modern or legacy')
  deju_assert(utils::packageVersion('edgeR') >= '4.0.0', 'edgeR >= 4.0.0 is required')
  if (engine=='modern') deju_assert(!is.null(getS3method('diffSplice','DGEGLM',optional=TRUE,envir=asNamespace('limma'))),
    'Modern engine requires edgeR with diffSplice.DGEGLM; install current edgeR/limma or explicitly choose --engine legacy')
  counts <- deju_validate_counts(counts,features,samples)
  design <- design_override %||DEJU% deju_design(samples,reference,treatment,covariates)
  deju_assert(all(colSums(counts)>0), 'A sample has a zero count library')
  y <- edgeR::DGEList(counts=counts,genes=features)
  expressed <- edgeR::filterByExpr(y,design=design$matrix,min.count=min_count,min.total.count=min_total_count)
  ng <- table(features$gene_id[expressed])
  gj <- unique(features$gene_id[expressed & features$feature_type=='junction'])
  eligible <- intersect(names(ng)[ng>=2L],gj)
  keep <- expressed & features$gene_id %in% eligible
  deju_assert(sum(keep)>=10 && length(eligible)>=3, 'Too few expressed features/genes for dispersion modelling (need >=10 features and >=3 genes)')
  filtering <- features[,c('feature_id','gene_id','feature_type')]
  filtering$expressed <- expressed; filtering$retained <- keep
  filtering$reason <- ifelse(!expressed,'low_expression',ifelse(!keep,'insufficient_gene_features_or_no_junction','retained'))
  y <- y[keep,,keep.lib.sizes=FALSE]
  deju_assert(all(colSums(y$counts)>0), 'A sample has no counts after filtering')
  normalize <- if ('normLibSizes' %in% getNamespaceExports('edgeR')) edgeR::normLibSizes else edgeR::calcNormFactors
  y <- normalize(y,method='TMM')
  y <- edgeR::estimateDisp(y,design$matrix,robust=TRUE)
  fit <- edgeR::glmQLFit(y,design$matrix,robust=TRUE,legacy=(engine=='legacy'))
  if (engine=='modern') {
    sp <- limma::diffSplice(fit,coef=design$coef %||DEJU% ncol(design$matrix),contrast=design$contrast,geneid='gene_id',exonid='feature_id',robust=TRUE,verbose=FALSE)
    pv <- as.numeric(sp$p.value)
    gs <- as.numeric(sp$gene.simes.p.value); gf <- as.numeric(sp$gene.F.p.value)
  } else {
    sp <- edgeR::diffSpliceDGE(fit,coef=design$coef %||DEJU% ncol(design$matrix),contrast=design$contrast,geneid='gene_id',exonid='feature_id',robust=TRUE,verbose=FALSE)
    pv <- as.numeric(sp$exon.p.value)
    gs <- as.numeric(sp$gene.Simes.p.value); gf <- as.numeric(sp$gene.p.value)
  }
  # edgeR's GLM and diffSplice coefficients are natural-log values, despite the
  # topSpliceDGE table calling them logFC. Convert explicitly, once, to log2.
  result <- sp$genes
  result$usage_log2FC <- as.numeric(sp$coefficients)/log(2)
  result$P.Value <- pv; result$all_feature_FDR <- p.adjust(pv,'BH')
  j <- result[result$feature_type=='junction',,drop=FALSE]
  j$junction_FDR <- p.adjust(j$P.Value,'BH')
  j$significant <- !is.na(j$junction_FDR) & j$junction_FDR<=fdr
  j$passes_effect_filter <- j$significant & abs(j$usage_log2FC)>=min_abs_log2
  e <- result[result$feature_type=='exon',,drop=FALSE]
  e$exon_FDR <- p.adjust(e$P.Value,'BH')
  genes <- sp$gene.genes
  names(genes)[names(genes)=='NExons'] <- 'n_features'
  genes$Simes_P.Value <- gs; genes$Simes_FDR <- p.adjust(gs,'BH')
  genes$F_P.Value <- gf; genes$F_FDR <- p.adjust(gf,'BH')
  # Membership in this list is not an additional gene-level hypothesis test.
  genes$has_significant_junction <- genes$gene_id %in% j$gene_id[j$significant]
  list(y=y,fit=fit,splice=sp,design=design,filtering=filtering,features=result,junctions=j,exons=e,genes=genes)
}

deju_positions <- function(junctions,annotation) {
  if (!all(c('chr','left','right','strand') %in% names(junctions))) return(NULL)
  keys <- c('gene_id','chr','left','right','strand')
  merge(junctions,annotation$junctions,by=keys,all.x=TRUE,sort=FALSE)
}

deju_go <- function(result,mapping_path,orgdb,outdir,fdr) {
  deju_require(c('clusterProfiler',orgdb))
  map <- read.delim(mapping_path,colClasses='character',check.names=FALSE)
  deju_assert(all(c('gene_id','ENTREZID') %in% names(map)) && !anyNA(map[,c('gene_id','ENTREZID')]),
              'GO mapping requires nonmissing gene_id and ENTREZID columns')
  map <- unique(map[,c('gene_id','ENTREZID')])
  universe <- unique(map$ENTREZID[map$gene_id %in% result$genes$gene_id])
  deju_assert(length(universe)>0, 'No tested genes map to ENTREZID')
  selections <- list(gene_Simes=result$genes$gene_id[result$genes$Simes_FDR<=fdr & !is.na(result$genes$Simes_FDR)],
                      significant_junction=result$junctions$gene_id[result$junctions$significant])
  db <- getExportedValue(orgdb,orgdb)
  for (name in names(selections)) {
    selected <- unique(map$ENTREZID[map$gene_id %in% selections[[name]]])
    deju_write(data.frame(ENTREZID=universe,selected=universe %in% selected),file.path(outdir,paste0('GO_',name,'_universe.tsv')))
    go <- if (length(selected)) clusterProfiler::enrichGO(gene=selected,universe=universe,OrgDb=db,keyType='ENTREZID',
              ont='BP',pAdjustMethod='BH',pvalueCutoff=1,qvalueCutoff=1,readable=FALSE) else NULL
    tab <- if (is.null(go)) data.frame(ID=character(),Description=character(),p.adjust=numeric()) else as.data.frame(go)
    deju_write(tab,file.path(outdir,paste0('GO_',name,'.tsv')))
  }
}

deju_plots <- function(r,outdir) {
  grDevices::pdf(file.path(outdir,'diagnostics.pdf'),width=8,height=6)
  on.exit(grDevices::dev.off())
  limma::plotMDS(r$y,labels=colnames(r$y),main='Joint exon/junction count MDS')
  edgeR::plotBCV(r$y)
  j <- r$junctions
  plot(j$usage_log2FC,-log10(pmax(j$junction_FDR,.Machine$double.xmin)),pch=16,cex=.5,
       col=ifelse(j$significant,'firebrick','grey50'),xlab='Relative usage change (log2; treatment - reference)',
       ylab='-log10 junction BH FDR',main='Differential junction usage')
}

deju_help <- function() cat('DEJUPipeline 0.2.0\n',
 'BAM analysis:\n  Rscript DEJUPipeline.R --samples samples.tsv --gtf annotation.gtf --paired true --strand 2 --reference WT --treatment KO --out results\n',
 'Count-table analysis:\n  Rscript DEJUPipeline.R --samples samples.tsv --counts counts.tsv --features features.tsv --reference WT --treatment KO --out results\n',
 'Options (each needs a value):\n',
 '  --engine modern|legacy (modern)   --covariates batch,subject\n',
 '  --threads 1 --min-mapq 10 --remove-duplicates false\n',
 '  --min-count 10 --min-total-count 15 --fdr 0.05 --min-abs-log2 0\n',
 '  --gene-map gene_to_entrez.tsv --go-orgdb org.Hs.eg.db (optional GO pair)\n',
 'strand: 0 unstranded, 1 forward, 2 reverse; paired: true or false.\n',
 'Output directory must not exist. See README for assumptions and file schemas.\n',sep='')

deju_cli <- function(args=commandArgs(TRUE)) {
  if (!length(args) || identical(args,'--help')) { deju_help(); return(invisible(NULL)) }
  allowed <- c('samples','gtf','paired','strand','reference','treatment','out','counts','features','engine','covariates',
               'threads','min-mapq','remove-duplicates','min-count','min-total-count','fdr','min-abs-log2','gene-map','go-orgdb')
  deju_assert(length(args)%%2==0,'Options require --name value pairs; use --help')
  keys <- sub('^--','',args[seq(1,length(args),2)])
  deju_assert(all(startsWith(args[seq(1,length(args),2)],'--')) && all(keys %in% allowed) && !anyDuplicated(keys),'Unknown or repeated option')
  o <- as.list(args[seq(2,length(args),2)]); names(o)<-keys
  defaults <- list(engine='modern',threads='1','min-mapq'='10','remove-duplicates'='false','min-count'='10',
                   'min-total-count'='15',fdr='0.05','min-abs-log2'='0')
  for (k in names(defaults)) if (is.null(o[[k]])) o[[k]]<-defaults[[k]]
  deju_assert(all(c('samples','reference','treatment','out') %in% names(o)),'Missing samples, reference, treatment or out')
  numeric_option <- function(k,low,high=Inf,integer=FALSE) {
    v <- suppressWarnings(as.numeric(o[[k]]))
    deju_assert(length(v)==1 && is.finite(v) && v>=low && v<=high && (!integer || v==round(v)),paste('Invalid',k))
    v
  }
  boolean <- function(k) { deju_assert(o[[k]] %in% c('true','false'),paste(k,'must be true or false')); o[[k]]=='true' }
  threads<-numeric_option('threads',1,64,TRUE); mapq<-numeric_option('min-mapq',0,255,TRUE)
  mincount<-numeric_option('min-count',0); mintotal<-numeric_option('min-total-count',0)
  fdr<-numeric_option('fdr',0,1); effect<-numeric_option('min-abs-log2',0)
  rmdup<-boolean('remove-duplicates')
  deju_assert(o$engine %in% c('modern','legacy'),'Invalid engine')
  count_mode <- !is.null(o$counts)
  deju_assert(count_mode==!is.null(o$features),'counts and features must be supplied together')
  deju_assert(is.null(o[['gene-map']])==is.null(o[['go-orgdb']]),'gene-map and go-orgdb must be supplied together')
  if (count_mode) deju_assert(!any(c('paired','strand','min-mapq','remove-duplicates','threads') %in% keys),
    'BAM counting options cannot change supplied count tables; remove them in counts mode')
  if (!count_mode) {
    deju_assert(all(c('gtf','paired','strand') %in% names(o)),'BAM mode requires gtf, paired and strand explicitly')
    paired<-boolean('paired'); strand<-numeric_option('strand',0,2,TRUE)
  }
  covars <- if (is.null(o$covariates)) character() else strsplit(o$covariates,',',fixed=TRUE)[[1]]
  samples <- deju_samples(o$samples,o$reference,o$treatment,covars,require_bams=!count_mode)
  deju_design(samples,o$reference,o$treatment,covars)
  deju_assert(!file.exists(o$out),'Output path already exists; choose a new output directory')
  deju_assert(dir.create(o$out,recursive=TRUE),'Cannot create output directory')
  writeLines('Analysis started; COMPLETE.txt is written only after success.',file.path(o$out,'INCOMPLETE.txt'))
  annotation <- if (!is.null(o$gtf)) prepare_deju_annotation(o$gtf) else NULL
  if (!is.null(annotation)) {
    deju_write(annotation$saf,file.path(o$out,'flattened_exons.saf.tsv.gz'))
    deju_write(annotation$junctions,file.path(o$out,'annotated_junctions.tsv.gz'))
  }
  if (count_mode) {
    d <- read.delim(o$counts,check.names=FALSE,row.names=1)
    f <- read.delim(o$features,check.names=FALSE,colClasses='character')
    input <- list(counts=as.matrix(d),features=f)
  } else input <- deju_count_bams(samples,annotation,o$out,paired,strand,threads,mapq,rmdup)
  r <- deju_fit(input$counts,input$features,samples,o$reference,o$treatment,covars,o$engine,mincount,mintotal,fdr,effect)
  deju_write(cbind(feature_id=rownames(input$counts),as.data.frame(input$counts)),file.path(o$out,'counts.tsv.gz'))
  deju_write(input$features,file.path(o$out,'features.tsv.gz'))
  for (name in c('filtering','features','junctions','exons','genes')) deju_write(r[[name]],file.path(o$out,paste0(name,'_results.tsv.gz')))
  deju_write(data.frame(sample_id=rownames(r$y$samples),r$y$samples),file.path(o$out,'normalization.tsv'))
  deju_write(data.frame(sample_id=rownames(r$design$matrix),r$design$matrix,check.names=FALSE),file.path(o$out,'design.tsv'))
  if (!is.null(annotation)) {
    p <- deju_positions(r$junctions,annotation)
    if (!is.null(p)) deju_write(p,file.path(o$out,'junction_transcript_positions.tsv.gz'))
  }
  deju_plots(r,o$out)
  if (!is.null(o[['gene-map']])) deju_go(r,o[['gene-map']],o[['go-orgdb']],o$out,fdr)
  saveRDS(list(version=DEJU_VERSION,settings=o,samples=samples,result=r),file.path(o$out,'analysis.rds'))
  capture.output(sessionInfo(),file=file.path(o$out,'sessionInfo.txt'))
  dput(o,file=file.path(o$out,'settings.R'))
  paths <- unlist(o[intersect(c('samples','gtf','counts','features','gene-map'),names(o))],use.names=FALSE)
  deju_write(data.frame(path=paths,md5=unname(tools::md5sum(paths))),file.path(o$out,'input_checksums.tsv'))
  if (!count_mode) deju_write(data.frame(sample_id=samples$sample_id,path=samples$bam,size_bytes=file.info(samples$bam)$size,
                                        modified=as.character(file.info(samples$bam)$mtime)),file.path(o$out,'bam_manifest.tsv'))
  msg <- sprintf('Tested %d genes and %d junctions; %d junctions at FDR <= %g. Contrast: %s - %s; engine: %s.',
    nrow(r$genes),nrow(r$junctions),sum(r$junctions$significant),fdr,o$treatment,o$reference,o$engine)
  writeLines(msg,file.path(o$out,'COMPLETE.txt')); unlink(file.path(o$out,'INCOMPLETE.txt'))
  message(msg); invisible(r)
}
if (sys.nframe()==0L) deju_cli()
