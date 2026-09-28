# Native statistical and counting modules for TranscriptomePipeline.py.
`%||%` <- function(x,y) if (is.null(x)) y else x

ta_samples <- function(cfg) {
  s <- read.delim(cfg$samples,check.names=FALSE,colClasses='character',na.strings=c('','NA'),quote='',comment.char='')
  deju_assert('sample_id' %in% names(s) && nrow(s)>=1 && !anyNA(s$sample_id) &&
    !anyDuplicated(s$sample_id) && all(grepl('^[A-Za-z0-9][A-Za-z0-9_.-]*$',s$sample_id)), 'Invalid/duplicated sample IDs or empty sample sheet')
  rownames(s)<-s$sample_id
  if ('bam' %in% names(s)) {
    relative<-!grepl('^(/|~)',s$bam)
    s$bam[relative]<-file.path(dirname(cfg$samples),s$bam[relative])
    s$bam<-path.expand(s$bam)
  }
  s
}

ta_design <- function(s,cfg) {
  spec<-cfg$design; deju_assert(!is.null(spec$formula),'A design formula is required')
  numeric_vars<-unlist(spec$numeric %||% list(),use.names=FALSE)
  factors<-spec$factors %||% list()
  deju_assert(!any(numeric_vars %in% names(factors)),'A variable cannot be both factor and numeric')
  # Only ordinary fixed-effect formula operators; no executable function calls.
  deju_assert(grepl('^[A-Za-z0-9_ .~+*:()/^-]+$',spec$formula) && !grepl('[A-Za-z_][A-Za-z0-9_.]*[[:space:]]*\\(',spec$formula),
              'Design supports fixed-effect terms/operators, not R function calls')
  form<-as.formula(spec$formula,env=baseenv()); vars<-all.vars(form)
  deju_assert(length(vars)>0 && all(vars %in% names(s)),'Design variables are missing from the sample sheet')
  deju_assert(setequal(vars,c(numeric_vars,names(factors))),'Declare every design variable in design.factors or design.numeric, and no unused variables')
  for (v in vars) {
    deju_assert(!anyNA(s[[v]]),paste('Missing values for',v))
    if (v %in% numeric_vars) {
      s[[v]]<-suppressWarnings(as.numeric(s[[v]]))
      deju_assert(all(is.finite(s[[v]])) && stats::sd(s[[v]])>0,paste('Invalid/constant numeric covariate',v))
    } else {
      levels<-unlist(factors[[v]],use.names=FALSE)
      deju_assert(length(levels)>=2 && !anyDuplicated(levels) && setequal(unique(s[[v]]),levels),paste('Factor levels do not match observed values:',v))
      s[[v]]<-factor(s[[v]],levels=levels)
    }
  }
  x<-model.matrix(form,s)
  deju_assert(qr(x)$rank==ncol(x),'Design is confounded or rank deficient')
  deju_assert(nrow(x)>ncol(x),'No residual degrees of freedom; biological replication is required')
  contrasts<-cfg$contrasts
  deju_assert(length(contrasts)>0,'Supply at least one named contrast')
  cn<-vapply(contrasts,function(z) z$name %||% '',character(1))
  deju_assert(!anyDuplicated(cn) && all(grepl('^[A-Za-z0-9][A-Za-z0-9_.-]*$',cn)),'Invalid/duplicate contrast names')
  vectors<-lapply(contrasts,function(z) {
    w<-unlist(z$weights)
    deju_assert(length(w)>0 && !is.null(names(w)) && !anyDuplicated(names(w)) && all(names(w) %in% colnames(x)) && is.numeric(w) && all(is.finite(w)),
      paste('Invalid contrast weights:',z$name,'; available coefficients:',paste(colnames(x),collapse=', ')))
    v<-setNames(rep(0,ncol(x)),colnames(x)); v[names(w)]<-w
    deju_assert(any(v!=0),paste('Zero contrast:',z$name));v
  });names(vectors)<-cn
  list(samples=s,matrix=x,contrasts=vectors,formula=form)
}

ta_read_counts <- function(path,samples) {
  sep<-if (grepl('[.]csv([.]gz)?$',path,ignore.case=TRUE)) ',' else '\t'
  d<-read.table(path,header=TRUE,sep=sep,check.names=FALSE,comment.char='#',quote='"',stringsAsFactors=FALSE)
  deju_assert(ncol(d)>1 && nrow(d)>0,'Empty count table')
  ids<-as.character(d[[1]])
  deju_assert(!anyNA(ids) && all(nzchar(ids)) && !anyDuplicated(ids),'Missing/duplicate count feature IDs')
  wanted<-if ('count_column' %in% names(samples)) samples$count_column else samples$sample_id
  deju_assert(!anyNA(wanted) && !anyDuplicated(wanted) && all(wanted %in% names(d)),
    'Count columns do not match sample_id; provide exact count_column values in the sample sheet for featureCounts/BAM column names')
  metadata<-c(names(d)[1],'Chr','Start','End','Strand','Length')
  deju_assert(all(names(d) %in% c(metadata,wanted)),'Unmapped columns in count table; include every sample or explicitly subset the input table')
  c<-as.matrix(d[,wanted,drop=FALSE]); storage.mode(c)<-'double'
  deju_assert(all(is.finite(c)) && all(c>=0) && all(c==round(c)) && all(c<=.Machine$integer.max),'Counts must be finite, nonnegative integers within R integer range; no TPM/FPKM or transformed values')
  rownames(c)<-ids;colnames(c)<-samples$sample_id
  deju_assert(all(colSums(c)>0),'Zero-count sample');c
}

ta_filter_bams <- function(s,cfg,out) {
  deju_require(c('Rsamtools','Rsubread'))
  deju_assert(all(file.exists(s$bam)) && !anyDuplicated(normalizePath(s$bam)),'Missing or duplicated BAMs')
  b<-cfg$bam; deju_assert(is.logical(b$paired) && length(b$paired)==1 && b$strand %in% 0:2,'Set bam.paired true/false and bam.strand 0/1/2')
  minq<-b$min_mapq %||% 10; rmdup<-b$remove_duplicates %||% FALSE
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  dest<-file.path(out,paste0(s$sample_id,'.bam')); qc<-list()
  for (i in seq_len(nrow(s))) {
    bf<-Rsamtools::BamFile(s$bam[i],yieldSize=250000L); open(bf)
    metrics<-c(total=0,mapped=0,secondary=0,supplementary=0,qc_fail=0,duplicate_marked=0,missing_NH=0,unique_NH=0,below_MAPQ=0,passing=0)
    tryCatch(repeat {
      z<-Rsamtools::scanBam(bf,param=Rsamtools::ScanBamParam(what=c('flag','mapq'),tag='NH'))[[1]]
      n<-length(z$flag);if (!n) break
      flag<-z$flag; mapped<-bitwAnd(flag,4L)==0
      paired<-bitwAnd(flag,1L)!=0
      deju_assert(all(paired[mapped]==b$paired),paste('Pairedness mismatch in',s$sample_id[i]))
      nh<-z$tag$NH;if (is.null(nh)) nh<-rep(NA_integer_,n)
      secondary<-bitwAnd(flag,256L)!=0;supp<-bitwAnd(flag,2048L)!=0
      badqc<-bitwAnd(flag,512L)!=0;dup<-bitwAnd(flag,1024L)!=0
      unique<-!is.na(nh) & nh==1
      passing<-mapped & !secondary & !supp & !badqc & unique & z$mapq>=minq & (!rmdup | !dup)
      metrics<-metrics+c(n,sum(mapped),sum(secondary),sum(supp),sum(badqc),sum(dup),sum(mapped & is.na(nh)),sum(mapped & unique),sum(mapped & z$mapq<minq),sum(passing))
    },finally=close(bf))
    deju_assert(metrics['passing']>0,paste('No alignments pass filters for',s$sample_id[i],'; inspect NH tags and MAPQ'))
    flags<-Rsamtools::scanBamFlag(isUnmappedQuery=FALSE,isSecondaryAlignment=FALSE,isSupplementaryAlignment=FALSE,
      isNotPassingQualityControls=FALSE,isDuplicate=if(rmdup) FALSE else NA)
    Rsamtools::filterBam(s$bam[i],dest[i],index=character(),param=Rsamtools::ScanBamParam(flag=flags,mapqFilter=minq,tagFilter=list(NH=1L)))
    qc[[i]]<-data.frame(sample_id=s$sample_id[i],as.list(metrics),pass_fraction=metrics['passing']/metrics['total'],row.names=NULL)
  }
  qc<-do.call(rbind,qc);deju_write(qc,file.path(dirname(out),'alignment_qc.tsv'))
  if(any(qc$pass_fraction<.1)) warning('Less than 10% of alignments pass in one or more samples; review alignment_qc.tsv')
  filtered<-s;filtered$bam<-dest;list(samples=filtered,qc=qc)
}

ta_gene_counts <- function(s,annotation,cfg,out) {
  b<-cfg$bam
  fc<-Rsubread::featureCounts(files=s$bam,annot.ext=annotation$saf,useMetaFeatures=TRUE,
    isGTFAnnotationFile=FALSE,strandSpecific=b$strand,isPairedEnd=b$paired,countReadPairs=b$paired,
    requireBothEndsMapped=b$paired,checkFragLength=FALSE,countChimericFragments=FALSE,
    countMultiMappingReads=FALSE,allowMultiOverlap=FALSE,nonSplitOnly=FALSE,juncCounts=FALSE,
    nthreads=cfg$threads %||% 1)
  c<-fc$counts[,deju_sample_columns(colnames(fc$counts),s$bam,s$sample_id),drop=FALSE]
  rownames(c)<-fc$annotation$GeneID;colnames(c)<-s$sample_id
  deju_write(cbind(gene_id=rownames(c),as.data.frame(c)),file.path(out,'gene_counts.tsv.gz'))
  stat<-fc$stat;names(stat)[-1]<-s$sample_id;deju_write(stat,file.path(out,'gene_assignment.tsv'))
  c
}

ta_dge <- function(c,design,cfg,out) {
  deju_require('DESeq2');dir.create(out,recursive=TRUE,showWarnings=FALSE)
  deju_assert(!any(grepl('^[EJ]:',rownames(c))),'DEJU feature counts are not gene counts. Use gene-level featureCounts output for DESeq2.')
  c<-c[,rownames(design$matrix),drop=FALSE]
  opt<-cfg$dge %||% list(); min_count<-opt$min_count %||% 10; min_samples<-opt$min_samples %||% 2
  keep<-rowSums(c>=min_count)>=min_samples
  deju_write(data.frame(gene_id=rownames(c),retained=keep,total_count=rowSums(c)),file.path(out,'gene_filtering.tsv'))
  deju_assert(sum(keep)>=10,'Fewer than 10 genes remain for DESeq2; review count units/filter settings')
  deju_assert(all(colSums(c[keep,,drop=FALSE])>0),'A sample is empty after gene filtering')
  storage.mode(c)<-'integer'
  dds<-DESeq2::DESeqDataSetFromMatrix(c[keep,,drop=FALSE],colData=design$samples,design=design$matrix)
  # Retain outlier evidence; no automatic count replacement.
  dds<-DESeq2::DESeq(dds,betaPrior=FALSE,minReplicatesForReplace=Inf,fitType=opt$fit_type %||% 'parametric',quiet=TRUE)
  normalized<-DESeq2::counts(dds,normalized=TRUE)
  deju_write(cbind(gene_id=rownames(normalized),as.data.frame(normalized)),file.path(out,'normalized_counts.tsv.gz'))
  sf<-DESeq2::sizeFactors(dds)
  deju_write(data.frame(sample_id=colnames(c),raw_gene_counts=colSums(c),detected_genes=colSums(c>0),size_factor=sf),file.path(out,'sample_qc.tsv'))
  deju_write(data.frame(matrix_coefficient=colnames(design$matrix),DESeq2_name=DESeq2::resultsNames(dds)),file.path(out,'coefficients.tsv'))
  vst<-DESeq2::varianceStabilizingTransformation(dds,blind=FALSE)
  v<-SummarizedExperiment::assay(vst)
  deju_write(cbind(gene_id=rownames(v),as.data.frame(v)),file.path(out,'vst_counts.tsv.gz'))
  deju_write(data.frame(sample_id=colnames(v),cor(v),check.names=FALSE),file.path(out,'sample_correlations.tsv'))
  covar<-SummarizedExperiment::assays(dds)[['cooks']]
  deju_write(cbind(gene_id=rownames(covar),as.data.frame(covar)),file.path(out,'cooks_distances.tsv.gz'))
  pdf(file.path(out,'QC.pdf'),width=9,height=7);on.exit(dev.off(),add=TRUE)
  barplot(colSums(c),names.arg=colnames(c),las=2,main='Raw gene-count library totals',ylab='Assigned counts')
  vars<-apply(v,1,var);top<-order(vars,decreasing=TRUE)[seq_len(min(500,nrow(v)))]
  pca<-prcomp(t(v[top,,drop=FALSE])); pct<-100*pca$sdev^2/sum(pca$sdev^2)
  plot(pca$x[,1:2],pch=16,xlab=sprintf('PC1 (%.1f%%)',pct[1]),ylab=sprintf('PC2 (%.1f%%)',pct[2]),main='VST sample PCA (unadjusted visualization)')
  text(pca$x[,1:2],labels=rownames(pca$x),pos=3,cex=.7)
  deju_write(data.frame(sample_id=rownames(pca$x),pca$x),file.path(out,'PCA_scores.tsv'))
  heatmap(as.matrix(dist(t(v))),symm=TRUE,main='VST sample distances')
  DESeq2::plotDispEsts(dds)
  summaries<-list();fdr<-cfg$fdr %||% .05
  for (name in names(design$contrasts)) {
    res<-DESeq2::results(dds,contrast=unname(design$contrasts[[name]]),alpha=fdr,independentFiltering=TRUE,cooksCutoff=TRUE)
    tab<-data.frame(gene_id=rownames(res),as.data.frame(res),row.names=NULL)
    tab$significant<-!is.na(tab$padj)&tab$padj<=fdr
    tab$status<-ifelse(is.na(tab$pvalue),'pvalue_unavailable_check_Cooks',ifelse(is.na(tab$padj),'independent_filter','tested'))
    deju_write(tab,file.path(out,paste0(name,'.tsv.gz')))
    DESeq2::plotMA(res,main=name)
    ok<-is.finite(tab$log2FoldChange)&!is.na(tab$padj)
    if(any(ok)) plot(tab$log2FoldChange[ok],-log10(pmax(tab$padj[ok],.Machine$double.xmin)),pch=16,cex=.5,
      col=ifelse(tab$significant[ok],'firebrick','grey50'),xlab='Gene expression log2 fold change (unshrunk)',ylab='-log10 BH adjusted p',main=name)
    summaries[[name]]<-data.frame(contrast=name,genes_retained=nrow(tab),genes_tested=sum(!is.na(tab$pvalue)),significant=sum(tab$significant),
      up=sum(tab$significant & tab$log2FoldChange>0,na.rm=TRUE),down=sum(tab$significant & tab$log2FoldChange<0,na.rm=TRUE))
  }
  deju_write(do.call(rbind,summaries),file.path(out,'summary.tsv'))
  saveRDS(dds,file.path(out,'DESeq2_model.rds'))
  invisible(dds)
}

ta_loci <- function(samples,cfg,annotation,out) {
  deju_require(c('Rsamtools','GenomicAlignments','GenomicRanges'));dir.create(out,recursive=TRUE,showWarnings=FALSE)
  loci<-cfg$loci; deju_assert(length(loci)>0,'loci module requires regions with name, chr, start, end')
  # Index aliases under the output directory, leaving source BAMs untouched.
  aliases<-file.path(out,paste0(samples$sample_id,'.bam'))
  for(i in seq_len(nrow(samples))) {
    deju_assert(file.symlink(normalizePath(samples$bam[i]),aliases[i]),'Cannot create BAM alias for locus indexing')
    Rsamtools::indexBam(aliases[i])
  }
  b<-cfg$bam;flags<-Rsamtools::scanBamFlag(isUnmappedQuery=FALSE,isSecondaryAlignment=FALSE,isSupplementaryAlignment=FALSE,
    isNotPassingQualityControls=FALSE,isDuplicate=if(isTRUE(b$remove_duplicates)) FALSE else NA)
  for(locus in loci) {
    deju_assert(locus$start>=1 && locus$end>locus$start && locus$end-locus$start<=1000000,'Locus must be 1-based and <= 1 Mb')
    region<-GenomicRanges::GRanges(locus$chr,IRanges::IRanges(locus$start,locus$end))
    pdf(file.path(out,paste0(locus$name,'.pdf')),width=12,height=max(5,2*nrow(samples)+2))
    tryCatch({
      par(mfrow=c(nrow(samples)+1,1),mar=c(2,4,2,1))
      for(i in seq_len(nrow(samples))) {
        param<-Rsamtools::ScanBamParam(which=region,flag=flags,mapqFilter=b$min_mapq %||% 10,tagFilter=list(NH=1L),what='flag')
        ga<-GenomicAlignments::readGAlignments(aliases[i],param=param)
        flag<-S4Vectors::mcols(ga)$flag
        if(!is.null(locus$strand) && b$strand!=0 && length(ga)) {
          neg<-bitwAnd(flag,16L)!=0
          if(b$paired) neg<-xor(neg,bitwAnd(flag,128L)!=0)
          if(b$strand==2) neg<-!neg
          ga<-ga[if(locus$strand=='+') !neg else neg]
        }
        blocks<-GenomicAlignments::grglist(ga)
        gr<-unlist(blocks,use.names=FALSE)
        coverage<-GenomicRanges::coverage(gr)
        cv<-numeric(locus$end-locus$start+1L)
        if(locus$chr %in% names(coverage)) {
          r<-coverage[[locus$chr]];hi<-min(length(r),locus$end)
          if(hi>=locus$start) cv[seq_len(hi-locus$start+1L)]<-as.numeric(r[locus$start:hi])
        }
        plot(seq(locus$start,locus$end),cv,type='l',xlim=c(locus$start,locus$end),ylim=c(-max(1,max(cv))*.7,max(1,max(cv))),
          xlab='',ylab='Read depth',yaxt='n',main=paste(samples$sample_id[i],'- raw aligned-read coverage'))
        ticks<-pretty(c(0,max(1,max(cv))));axis(2,at=ticks[ticks>=0 & ticks<=max(1,max(cv))])
        intr<-GenomicAlignments::cigarRangesAlongReferenceSpace(GenomicAlignments::cigar(ga),pos=GenomicRanges::start(ga),ops='N')
        ir<-unlist(intr,use.names=FALSE)
        if(length(ir)) {
          tab<-table(paste(IRanges::start(ir)-1L,IRanges::end(ir)+1L,sep=':'))
          for(k in names(tab)) {
            ends<-as.integer(strsplit(k,':',fixed=TRUE)[[1]])
            if(ends[1]<locus$start || ends[2]>locus$end) next
            x<-seq(ends[1],ends[2],length.out=50);h<-max(1,max(cv))*.45
            y<- -h*sin(seq(0,pi,length.out=50))
            lines(x,y,lwd=1+log1p(tab[[k]])/2,col='steelblue')
            text(mean(ends),-h,tab[[k]],cex=.6,pos=1)
          }
        }
      }
      ex<-annotation$saf;ex<-ex[ex$Chr==locus$chr & ex$End>=locus$start & ex$Start<=locus$end,,drop=FALSE]
      genes<-unique(ex$GeneID)
      plot(NA,xlim=c(locus$start,locus$end),ylim=c(0,max(1,length(genes))+1),xlab=paste(locus$chr,'(1-based)'),ylab='',yaxt='n',main='Merged annotated exons by gene')
      for(g in seq_along(genes)) {
        e<-ex[ex$GeneID==genes[g],];rect(e$Start,g-.15,e$End,g+.15,col='grey30');text(locus$start,g,genes[g],pos=4,cex=.6)
      }
    },finally=dev.off())
  }
  unlink(c(aliases,paste0(aliases,'.bai')))
}

ta_main <- function(path,check=FALSE) {
  cfg<-jsonlite::fromJSON(path,simplifyVector=FALSE)
  modules<-unlist(cfg$modules,use.names=FALSE);s<-ta_samples(cfg)
  statistical<-any(c('dge','dju') %in% modules)
  design<-if(statistical) ta_design(s,cfg) else NULL
  if(check) {
    if(('dge' %in% modules && is.null(cfg$gene_counts)) || ('dju' %in% modules && is.null(cfg$dju_counts)) || 'loci' %in% modules)
      deju_require(c('Rsubread','Rsamtools','GenomicRanges','GenomicAlignments','IRanges'))
    if(statistical) cat('Design coefficients:',paste(colnames(design$matrix),collapse=', '),'\n')
    for(m in intersect(modules,c('dge','dju'))) deju_require(if(m=='dge') 'DESeq2' else c('edgeR','limma','statmod'))
    if('gene_counts' %in% names(cfg)) ta_read_counts(cfg$gene_counts,s)
    if('dju_counts' %in% names(cfg)) ta_read_counts(cfg$dju_counts,s)
    cat('Native preflight passed\n');return(invisible(NULL))
  }
  out<-cfg$out;dir.create(file.path(out,'native'),recursive=TRUE,showWarnings=FALSE)
  if(statistical) {
    deju_write(data.frame(sample_id=rownames(design$matrix),design$matrix,check.names=FALSE),file.path(out,'native','design.tsv'))
    deju_write(data.frame(coefficient=colnames(design$matrix),do.call(cbind,design$contrasts),check.names=FALSE),file.path(out,'native','contrasts.tsv'))
  }
  needs_bam<-('dge' %in% modules && is.null(cfg$gene_counts)) || ('dju' %in% modules && is.null(cfg$dju_counts))
  annotation<-if(needs_bam || 'loci' %in% modules || ('dju' %in% modules && !is.null(cfg$gtf))) prepare_deju_annotation(cfg$gtf) else NULL
  if(needs_bam) {
    scratch<-file.path(out,'native','filtered_bams');on.exit(unlink(scratch,recursive=TRUE),add=TRUE)
    filtered<-ta_filter_bams(s,cfg,scratch)
  }
  if('dge' %in% modules) {
    c<-if(!is.null(cfg$gene_counts)) ta_read_counts(cfg$gene_counts,s) else ta_gene_counts(filtered$samples,annotation,cfg,file.path(out,'native'))
    ta_dge(c,design,cfg,file.path(out,'dge'))
  }
  if('dju' %in% modules) {
    opt<-cfg$dju %||% list()
    if(!is.null(cfg$dju_counts)) {
      input<-list(counts=ta_read_counts(cfg$dju_counts,s),features=read.delim(cfg$dju_features,check.names=FALSE,colClasses='character'))
    } else {
      path<-file.path(out,'junction_counts');dir.create(path)
      input<-deju_count_bams(filtered$samples,annotation,path,cfg$bam$paired,cfg$bam$strand,threads=cfg$threads %||% 1,pre_filtered=TRUE,
        min_anchor=opt$min_anchor %||% 8,min_intron=opt$min_intron %||% 20,max_intron=opt$max_intron %||% 1000000)
      deju_write(cbind(feature_id=rownames(input$counts),as.data.frame(input$counts)),file.path(path,'counts.tsv.gz'))
      deju_write(input$features,file.path(path,'features.tsv.gz'))
    }
    for(name in names(design$contrasts)) {
      path<-file.path(out,'dju',name);dir.create(path,recursive=TRUE)
      r<-deju_fit(input$counts,input$features,s,NULL,NULL,engine=opt$engine %||% 'modern',min_count=opt$min_count %||% 10,
        min_total_count=opt$min_total_count %||% 15,fdr=cfg$fdr %||% .05,
        design_override=list(matrix=design$matrix,contrast=unname(design$contrasts[[name]])))
      for(k in c('junctions','exons','genes','filtering')) deju_write(r[[k]],file.path(path,paste0(k,'.tsv.gz')))
      deju_plots(r,path);saveRDS(r,file.path(path,'model.rds'))
      if(!is.null(annotation)) {p<-deju_positions(r$junctions,annotation);if(!is.null(p)) deju_write(p,file.path(path,'transcript_positions.tsv.gz'))}
    }
  }
  if('loci' %in% modules) ta_loci(s,cfg,annotation,file.path(out,'loci'))
  capture.output(sessionInfo(),file=file.path(out,'native','sessionInfo.txt'))
  writeLines('Native modules completed',file.path(out,'native','COMPLETE.txt'))
}
