#!/usr/bin/env Rscript
arg<-grep('^--file=',commandArgs(FALSE),value=TRUE)
root<-dirname(dirname(normalizePath(sub('^--file=','',arg[1]))))
source(file.path(root,'DEJUPipeline.R'));source(file.path(root,'modules','AnalysisCore.R'))
source(file.path(root,'tests','make_fixtures.R'))
assert_error<-function(expr) stopifnot(inherits(tryCatch({force(expr);NULL},error=identity),'error'))
tmp<-tempfile('transcriptome-v2-');dir.create(tmp)
cfg<-jsonlite::fromJSON(file.path(root,'examples/integrated/config.json'),simplifyVector=FALSE)
cfg$samples<-file.path(root,'examples/integrated/samples.tsv')
s<-ta_samples(cfg);d<-ta_design(s,cfg)
stopifnot(identical(colnames(d$matrix),c('(Intercept)','batchB','conditionKO','time','conditionKO:time')))
c<-ta_read_counts(file.path(root,'examples/integrated/gene_counts.tsv'),s)
r<-ta_dge(c,d,cfg,file.path(tmp,'dge'))
res<-read.delim(gzfile(file.path(tmp,'dge/time_by_condition_interaction.tsv.gz')))
stopifnot(all(res$padj[match(sprintf('Gene%03d',11:15),res$gene_id)]<.05),
          all(res$log2FoldChange[match(sprintf('Gene%03d',11:15),res$gene_id)]>1))
# Column order invariant on real DESeq2 fitted object construction.
c2<-c[,rev(seq_len(ncol(c))),drop=FALSE]
stopifnot(identical(c2[,rownames(d$matrix),drop=FALSE],c))
bad<-cfg;bad$design$formula<-'~ batch + condition + time';assert_error(ta_design(s,bad))
bad<-cfg;bad$design$formula<-'~ condition + condition2';bad$design$factors<-list(condition=c('WT','KO'),condition2=c('WT','KO'));bad$design$numeric<-list()
s2<-s;s2$condition2<-s2$condition;assert_error(ta_design(s2,bad))
# Shared interaction contrast in DEJU on independently simulated feature counts.
set.seed(229)
f<-expand.grid(part=c('E1','E2','J1','J2'),gene_id=sprintf('G%03d',1:100),stringsAsFactors=FALSE)
f$feature_id<-paste(f$gene_id,f$part,sep=':');f$feature_type<-ifelse(startsWith(f$part,'J'),'junction','exon')
mu<-matrix(rep(runif(nrow(f),200,500),nrow(s)),nrow(f),nrow(s))
mu[f$gene_id=='G001' & f$part=='J1',s$condition=='KO' & s$time=='1']<-4000
fc<-matrix(rnbinom(length(mu),mu=mu,size=90),nrow(f),nrow(s),dimnames=list(f$feature_id,s$sample_id))
ju<-deju_fit(fc,f,s,NULL,NULL,design_override=list(matrix=d$matrix,contrast=unname(d$contrasts$time_by_condition_interaction)))
j<-ju$junctions[ju$junctions$feature_id=='G001:J1',]
stopifnot(j$junction_FDR<.05,j$usage_log2FC>1)
# Three-condition design, explicit continuous covariate and named contrast.
multi<-data.frame(sample_id=paste0('M',1:9),group=rep(c('A','B','C'),each=3),age=rep(c(20,30,40),3))
mc<-list(design=list(formula='~ group + age',factors=list(group=c('A','B','C')),numeric=list('age')),
         contrasts=list(list(name='C_vs_B',weights=list(groupC=1,groupB=-1))))
md<-ta_design(multi,mc);stopifnot(identical(unname(md$contrasts$C_vs_B),c(0,-1,1,0)))
cat('PASS: DESeq2 condition/interaction, shared DEJU contrast, three conditions, continuous covariates and confounding rejection\n')
if('--bam' %in% commandArgs(TRUE)) {
  path<-file.path(tmp,'bams');make_bam_fixture(path,TRUE,TRUE)
  bc<-list(samples=file.path(path,'samples.tsv'),gtf=file.path(path,'annotation.gtf'),out=file.path(tmp,'bam_output'),
    modules=list('dge','dju','loci'),bam=list(paired=TRUE,strand=1,min_mapq=10,remove_duplicates=TRUE),
    design=list(formula='~ condition',factors=list(condition=list('WT','KO')),numeric=list()),
    contrasts=list(list(name='KO_vs_WT',weights=list(conditionKO=1))),dge=list(fit_type='mean'),
    loci=list(list(name='G001',chr='chrTest',start=990,end=1400,strand='+')))
  file<-file.path(tmp,'bam.json');jsonlite::write_json(bc,file,auto_unbox=TRUE)
  ta_main(file)
  count<-as.matrix(read.delim(gzfile(file.path(bc$out,'native/gene_counts.tsv.gz')),row.names=1))
  qc<-read.delim(file.path(bc$out,'native/alignment_qc.tsv'))
  stopifnot(all(colSums(count)==qc$passing/2),file.exists(file.path(bc$out,'loci/G001.pdf')),
            file.exists(file.path(bc$out,'native/COMPLETE.txt')),!dir.exists(file.path(bc$out,'native/filtered_bams')))
  cat('PASS: independent fragment gene counts, read-level junctions, shared BAM QC, full native workflow and locus PDF\n')
}
unlink(tmp,recursive=TRUE)
