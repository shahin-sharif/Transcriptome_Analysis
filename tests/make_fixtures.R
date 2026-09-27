# Deterministic synthetic inputs only; no biological/private data.
make_count_fixture <- function(path) {
  dir.create(path,recursive=TRUE,showWarnings=FALSE)
  set.seed(20260927)
  samples <- data.frame(sample_id=paste0('S',1:8),condition=rep(c('WT','KO'),each=4),batch=rep(c('A','B'),4))
  f <- expand.grid(part=c('E1','E2','J1','J2'),gene_id=sprintf('G%03d',1:120),stringsAsFactors=FALSE)
  f$feature_id<-paste(f$gene_id,f$part,sep=':'); f$feature_type<-ifelse(startsWith(f$part,'E'),'exon','junction')
  f<-f[,c('feature_id','gene_id','feature_type','part')]
  baseline<-rep(runif(120,150,500),each=4)*rep(c(1,1.2,.7,.9),120)
  mu<-outer(baseline,c(.8,1,1.1,.9,1,.85,1.2,1.05))
  # Relative usage changes and whole-gene expression changes are distinct controls.
  mu[f$gene_id=='G001' & f$part=='J1',5:8]<-mu[f$gene_id=='G001' & f$part=='J1',5:8]*12
  mu[f$gene_id=='G002' & f$part=='J1',5:8]<-mu[f$gene_id=='G002' & f$part=='J1',5:8]/12
  mu[f$gene_id=='G003',5:8]<-mu[f$gene_id=='G003',5:8]*5
  counts<-matrix(rnbinom(length(mu),mu=mu,size=80),nrow(mu),dimnames=list(f$feature_id,samples$sample_id))
  write.table(samples,file.path(path,'samples.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  write.table(f,file.path(path,'features.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  write.table(data.frame(feature_id=rownames(counts),counts),file.path(path,'counts.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
}

make_bam_fixture <- function(path,paired=FALSE) {
  dir.create(path,recursive=TRUE,showWarnings=FALSE)
  gtf <- character()
  for (g in 1:24) {
    start<-g*1000
    for (e in 0:2) gtf<-c(gtf,paste('chrTest','fixture','exon',start+e*150,start+e*150+49,'.','+','.',
      sprintf('gene_id "G%03d"; transcript_id "T%03d";',g,g),sep='\t'))
  }
  writeLines(gtf,file.path(path,'annotation.gtf'))
  samples<-data.frame(sample_id=paste0('S',1:6),condition=rep(c('WT','KO'),each=3),bam=paste0('S',1:6,'.bam'))
  set.seed(131)
  for (s in 1:6) {
    lines<-c('@HD\tVN:1.6\tSO:unsorted','@SQ\tSN:chrTest\tLN:100000')
    k<-0L
    add <- function(pos,cigar,nh=1,flag=0L,mapq=60L) {
      k<<-k+1L
      if (paired) {
        # Both mates support the SAME junction; verifies read vs fragment counting.
        lines<<-c(lines,paste(paste0('r',k),99+flag,'chrTest',pos,mapq,cigar,'=',pos,20,strrep('A',20),strrep('I',20),paste0('NH:i:',nh),sep='\t'),
                       paste(paste0('r',k),147+flag,'chrTest',pos,mapq,cigar,'=',pos,-20,strrep('T',20),strrep('I',20),paste0('NH:i:',nh),sep='\t'))
      } else lines<<-c(lines,paste(paste0('r',k),flag,'chrTest',pos,mapq,cigar,'*',0,0,strrep('A',20),strrep('I',20),paste0('NH:i:',nh),sep='\t'))
    }
    for (g in 1:24) {
      start<-g*1000
      for (e in 0:2) for (z in seq_len(30L+(g+s+e)%%8)) add(start+e*150+10,'20M')
      for (e in 0:1) {
        n<-30L+(g+s+e)%%8
        if (g==1 && e==0 && s>3) n<-n*8L
        for (z in seq_len(n)) add(start+e*150+40,'10M100N10M')
      }
    }
    # Uniqueness, primary/supplementary, quality and duplicate policy controls.
    for (z in 1:5) add(1040,'10M100N10M',nh=2)
    for (flag in c(256L,2048L,512L,1024L)) add(1040,'10M100N10M',flag=flag)
    add(1040,'10M100N10M',mapq=0L)
    sam<-file.path(path,paste0('S',s,'.sam')); writeLines(lines,sam)
    Rsamtools::asBam(sam,destination=file.path(path,paste0('S',s)),overwrite=TRUE)
    unlink(sam)
  }
  write.table(samples,file.path(path,'samples.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
}
