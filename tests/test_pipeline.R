#!/usr/bin/env Rscript
arg <- grep('^--file=',commandArgs(FALSE),value=TRUE)
root <- dirname(dirname(normalizePath(sub('^--file=','',arg[1]))))
source(file.path(root,'DEJUPipeline.R'))
source(file.path(root,'tests','make_fixtures.R'))
assert_error <- function(expr) stopifnot(inherits(tryCatch({force(expr);NULL},error=identity),'error'))
tmp<-tempfile('deju-tests-'); dir.create(tmp)
make_count_fixture(tmp)
s<-deju_samples(file.path(tmp,'samples.tsv'),'WT','KO',require_bams=FALSE)
f<-read.delim(file.path(tmp,'features.tsv'))
c<-as.matrix(read.delim(file.path(tmp,'counts.tsv'),row.names=1,check.names=FALSE))
r<-deju_fit(c,f,s,'WT','KO',covariates='batch')
j<-r$junctions
stopifnot(j$usage_log2FC[j$feature_id=='G001:J1']>1,j$junction_FDR[j$feature_id=='G001:J1']<.05,
          j$usage_log2FC[j$feature_id=='G002:J1']< -1,j$junction_FDR[j$feature_id=='G002:J1']<.05,
          all(j$junction_FDR[j$gene_id=='G003']>.05),
          isTRUE(all.equal(j$junction_FDR,p.adjust(j$P.Value,'BH'))),nrow(r$genes)==120)
permuted<-deju_fit(c[rev(seq_len(nrow(c))),rev(seq_len(ncol(c)))],f,s,'WT','KO',covariates='batch')
stopifnot(isTRUE(all.equal(r$junctions,permuted$junctions)))
legacy<-deju_fit(c,f,s,'WT','KO',engine='legacy')
stopifnot(legacy$junctions$usage_log2FC[legacy$junctions$feature_id=='G001:J1']>1)
bad<-c; bad[1,1]<- -1; assert_error(deju_validate_counts(bad,f,s))
bad<-c; colnames(bad)[1]<-'unknown'; assert_error(deju_validate_counts(bad,f,s))
bad<-s; bad$batch<-bad$condition; assert_error(deju_design(bad,'WT','KO','batch'))
stopifnot(identical(deju_sample_columns(c('meta','/b/S2.bam','/b/S1.bam'),c('/a/S1.bam','/a/S2.bam'),c('S1','S2')),c(3L,2L)))
# Minus-strand coordinate identity and all transcript positions.
gtf<-file.path(tmp,'tiny.gtf')
writeLines(c('chr1\tx\texon\t101\t150\t.\t-\t.\tgene_id "M"; transcript_id "MT";',
             'chr1\tx\texon\t201\t250\t.\t-\t.\tgene_id "M"; transcript_id "MT";',
             'chr1\tx\texon\t301\t350\t.\t-\t.\tgene_id "M"; transcript_id "MT";'),gtf)
a<-prepare_deju_annotation(gtf)
stopifnot(identical(a$junctions$left,c(150L,250L)),identical(a$junctions$donor,c(201L,301L)),
          identical(a$junctions$junction_rank,c(2L,1L)))
x<-data.frame(chr='chr1',chr2='chr1',left=c(150,160,20),right=c(201,240,80),primary=NA_character_,secondary=NA_character_)
ja<-deju_assign_junctions(x,a)
stopifnot(identical(ja$assignment,c('annotated','novel_intragenic','unassigned')))
# Stranded intergenic entries must not delete list slots or shift later genes.
stranded <- x[c(3,1,2,3),]
stranded$read_strand <- '-'
assigned <- deju_assign_junctions(stranded,a)
stopifnot(identical(assigned$assignment,c('unassigned','annotated','novel_intragenic','unassigned')),
          identical(assigned$gene_id,c(NA_character_,'M','M',NA_character_)))
stranded$read_strand <- '+'
stopifnot(all(deju_assign_junctions(stranded,a)$assignment=='unassigned'))
a$genes<-rbind(a$genes,transform(a$genes,gene_id='OTHER'))
stopifnot(deju_assign_junctions(x,a)$assignment[2]=='ambiguous_gene')
cli_out<-file.path(tmp,'cli_counts')
script<-file.path(root,'DEJUPipeline.R')
run_cli<-function(args) {
  status<-system2(file.path(R.home('bin'),'Rscript'),c('--vanilla',shQuote(script),vapply(args,shQuote,character(1))),
                  stdout=file.path(tmp,'cli.log'),stderr=file.path(tmp,'cli.log'))
  if (status!=0) stop(paste(readLines(file.path(tmp,'cli.log')),collapse='\n'))
}
run_cli(c('--samples',file.path(tmp,'samples.tsv'),'--counts',file.path(tmp,'counts.tsv'),
          '--features',file.path(tmp,'features.tsv'),'--reference','WT','--treatment','KO','--out',cli_out))
stopifnot(file.exists(file.path(cli_out,'COMPLETE.txt')),file.exists(file.path(cli_out,'analysis.rds')),
          file.info(file.path(cli_out,'diagnostics.pdf'))$size>1000,!file.exists(file.path(cli_out,'INCOMPLETE.txt')))
assert_error(deju_cli(c('--samples',file.path(tmp,'samples.tsv'),'--counts',file.path(tmp,'counts.tsv'),
          '--features',file.path(tmp,'features.tsv'),'--reference','WT','--treatment','KO','--out',cli_out)))
cat('PASS: full count-table CLI, output artifacts and overwrite protection\n')
# Multiple splice gaps, soft clipping and deletions must advance the reference correctly.
sam<-file.path(tmp,'cigar.sam')
writeLines(c('@HD\tVN:1.6\tSO:unsorted','@SQ\tSN:chr1\tLN:1000',
  paste('complex',16,'chr1',100,60,'5S10M100N5M2D5M100N10M','*',0,0,strrep('A',35),strrep('I',35),'NH:i:1',sep='\t')),sam)
Rsamtools::asBam(sam,destination=file.path(tmp,'cigar'),overwrite=TRUE)
cigar<-deju_junction_counts(file.path(tmp,'cigar.bam'),'one',FALSE,2,yield_size=1L,min_anchor=0L)
stopifnot(identical(cigar$annotation$left,c(109L,221L)),identical(cigar$annotation$right,c(210L,322L)),
          all(cigar$counts==1),all(cigar$annotation$read_strand=='+'))
cat('PASS: multi-junction CIGAR, clipping, deletion and reverse-strand coordinates\n')

cat('PASS: count models, effect direction, expression-only control, BH families, sample order, validation and annotation\n')
if ('--bam' %in% commandArgs(TRUE)) {
  for (paired in c(FALSE,TRUE)) {
    path<-file.path(tmp,if (paired) 'paired' else 'single'); make_bam_fixture(path,paired)
    samples<-deju_samples(file.path(path,'samples.tsv'),'WT','KO')
    ann<-prepare_deju_annotation(file.path(path,'annotation.gtf'))
    out<-file.path(path,'counted');dir.create(out)
    d<-deju_count_bams(samples,ann,out,paired,0,remove_duplicates=TRUE)
    expected<-(30+(1+1)%%8)*if (paired) 2 else 1
    observed<-d$counts['J:G001:chrTest:1049:1150:+','S1']
    stopifnot(observed==expected)
    forward<-deju_junction_counts(samples$bam,samples$sample_id,paired,1,yield_size=37L)
    reverse<-deju_junction_counts(samples$bam,samples$sample_id,paired,2,yield_size=37L)
    stopifnot(all(forward$annotation$read_strand=='+'),all(reverse$annotation$read_strand=='-'),
              identical(forward$counts,reverse$counts))
    stopifnot(all(is.na(deju_assign_junctions(reverse$annotation,ann)$gene_id)))
    assert_error(deju_junction_counts(samples$bam,samples$sample_id,!paired,0))
    r<-deju_fit(d$counts,d$features,samples,'WT','KO')
    stopifnot(nrow(r$junctions)==48,r$junctions$usage_log2FC[r$junctions$feature_id=='J:G001:chrTest:1049:1150:+']>1)
    if (paired) {
      bam_out<-file.path(path,'cli_bams')
      run_cli(c('--samples',file.path(path,'samples.tsv'),'--gtf',file.path(path,'annotation.gtf'),
        '--paired','true','--strand','1','--remove-duplicates','true','--reference','WT','--treatment','KO','--out',bam_out))
      stopifnot(file.exists(file.path(bam_out,'COMPLETE.txt')),!dir.exists(file.path(bam_out,'counting_tmp')))
      positions<-read.delim(gzfile(file.path(bam_out,'junction_transcript_positions.tsv.gz')))
      stopifnot(nrow(positions)==48,all(!is.na(positions$transcript_id)))
      cat('PASS: full paired BAM CLI and transcript positions\n')
    }
    cat('PASS: native BAM counting, filters and fit; paired =',paired,'\n')
  }
}
if ('--go' %in% commandArgs(TRUE)) {
  path<-file.path(tmp,'go');dir.create(path)
  # Arbitrary synthetic mapping: execution test, not an enrichment result to interpret.
  deju_write(data.frame(gene_id=sprintf('G%03d',1:120),ENTREZID=as.character(1:120)),file.path(path,'mapping.tsv'))
  deju_go(permuted,file.path(path,'mapping.tsv'),'org.Hs.eg.db',path,.05)
  stopifnot(file.exists(file.path(path,'GO_gene_Simes.tsv')),file.exists(file.path(path,'GO_significant_junction_universe.tsv')))
  cat('PASS: optional GO and universe exports (synthetic mapping only)\n')
}
unlink(tmp,recursive=TRUE)
