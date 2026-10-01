# Helpers for Salmon -> IsoformSwitchAnalyzeR/DEXSeq. Sourcing does not run analysis.
is_assert <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call.=FALSE)
is_default <- function(x,y) if (is.null(x)) y else x
is_write <- function(x,p) utils::write.table(x,p,sep='\t',quote=FALSE,row.names=FALSE,na='NA')
is_id <- function(x, bar=FALSE) {
  if (bar) x <- sub('\\|.*$', '', x)
  is_assert(!anyNA(x) && all(nzchar(x)) && !any(grepl('[[:space:]]',x)) && !anyDuplicated(x),
            'Missing/duplicate transcript IDs, including collisions after configured pipe removal')
  x
}
is_path <- function(p, base) {
  is_assert(is.character(p) && length(p)==1 && nzchar(p),'Invalid path')
  p <- path.expand(p)
  if (!grepl('^/',p)) p <- file.path(base,p)
  normalizePath(p,mustWork=TRUE)
}
is_config <- function(path) {
  is_assert(requireNamespace('jsonlite',quietly=TRUE),'Install jsonlite')
  cfg <- jsonlite::fromJSON(path,simplifyVector=FALSE)
  allowed <- c('samples','out','gtf','transcript_fasta','comparisons','covariates','alpha','delta_if',
               'filter','strip_pipe','consequences','predict_novel_orfs','plots','go')
  is_assert(!length(setdiff(names(cfg),allowed)),'Unknown isoform-switch configuration option')
  base <- dirname(normalizePath(path,mustWork=TRUE))
  for (k in c('samples','gtf','transcript_fasta')) if (!is.null(cfg[[k]])) cfg[[k]] <- is_path(cfg[[k]],base)
  is_assert(!is.null(cfg$samples) && !is.null(cfg$gtf) && !is.null(cfg$out),'samples, gtf and out are required')
  cfg$out <- path.expand(cfg$out)
  if (!grepl('^/',cfg$out)) cfg$out <- file.path(base,cfg$out)
  cfg$alpha <- is_default(cfg$alpha,.05);cfg$delta_if <- is_default(cfg$delta_if,.1)
  is_assert(is.numeric(cfg$alpha) && cfg$alpha>0 && cfg$alpha<1,'alpha must be in (0,1)')
  is_assert(is.numeric(cfg$delta_if) && cfg$delta_if>=0 && cfg$delta_if<=1,'delta_if must be in [0,1]')
  for (k in c('strip_pipe','consequences','predict_novel_orfs')) {
    cfg[[k]] <- is_default(cfg[[k]],FALSE)
    is_assert(is.logical(cfg[[k]]) && length(cfg[[k]])==1 && !is.na(cfg[[k]]),paste('Invalid boolean:',k))
  }
  cfg$plots <- is_default(cfg$plots,10L)
  is_assert(is.numeric(cfg$plots) && cfg$plots>=0 && cfg$plots==as.integer(cfg$plots),'plots must be a nonnegative integer')
  if (cfg$consequences || cfg$predict_novel_orfs) is_assert(!is.null(cfg$transcript_fasta),
    'Consequence/novel-ORF analysis needs the matched transcript FASTA, not the genome FASTA')
  is_assert(!cfg$predict_novel_orfs || cfg$consequences,'predict_novel_orfs requires consequences=true')
  filt <- is_default(cfg$filter,list())
  is_assert(!length(setdiff(names(filt),c('iso_count','count_proportion','if_cutoff','if_proportion'))),'Unknown filter option')
  cfg$filter <- modifyList(list(iso_count=10,count_proportion=.7,if_cutoff=.01,if_proportion=.5),filt)
  is_assert(is.numeric(cfg$filter$iso_count) && cfg$filter$iso_count>=0,'Invalid iso_count')
  for(k in c('count_proportion','if_proportion')) is_assert(is.numeric(cfg$filter[[k]]) && cfg$filter[[k]]>0 && cfg$filter[[k]]<=1,paste('Invalid',k))
  is_assert(is.numeric(cfg$filter$if_cutoff) && cfg$filter$if_cutoff>=0 && cfg$filter$if_cutoff<=1,'Invalid if_cutoff')
  if(!is.null(cfg$go)) {
    is_assert(all(c('orgdb','mapping') %in% names(cfg$go)) && !length(setdiff(names(cfg$go),c('orgdb','mapping','ontology'))),
              'go needs orgdb and mapping (optional ontology)')
    cfg$go$mapping <- is_path(cfg$go$mapping,base)
    cfg$go$ontology <- is_default(cfg$go$ontology,'BP')
    is_assert(cfg$go$ontology %in% c('BP','MF','CC'),'Invalid GO ontology')
  }
  cfg
}
is_design <- function(cfg) {
  meta <- utils::read.delim(cfg$samples,check.names=FALSE,stringsAsFactors=FALSE,colClasses='character')
  is_assert(all(c('sample_id','condition','quant') %in% names(meta)),'Samples need sample_id, condition, quant columns')
  is_assert(!anyDuplicated(names(meta)) && nrow(meta)>=4,'Duplicate columns or insufficient samples')
  is_assert(all(grepl('^[A-Za-z0-9][A-Za-z0-9_.-]*$',meta$sample_id)) && !anyDuplicated(meta$sample_id),'Invalid sample IDs')
  meta$quant <- vapply(meta$quant,is_path,character(1),base=dirname(cfg$samples))
  is_assert(!anyDuplicated(meta$quant),'A quantification cannot be reused as a biological replicate')
  is_assert(!anyNA(meta$condition) && all(nzchar(meta$condition)),'Missing condition')
  is_assert(length(cfg$comparisons)>0,'Explicit comparisons are required')
  comparisons <- do.call(rbind,lapply(cfg$comparisons,function(x) {
    is_assert(setequal(names(x),c('name','reference','treatment')),'Each comparison needs name, reference, treatment')
    is_assert(grepl('^[A-Za-z0-9][A-Za-z0-9_.-]*$',x$name) && x$reference!=x$treatment,'Invalid comparison')
    is_assert(all(c(x$reference,x$treatment) %in% meta$condition),'Unknown comparison condition')
    data.frame(name=x$name,condition_1=x$reference,condition_2=x$treatment,stringsAsFactors=FALSE)
  }))
  is_assert(!anyDuplicated(comparisons$name) && !anyDuplicated(comparisons[,2:3]),'Duplicate comparisons')
  needed <- unique(c(comparisons$condition_1,comparisons$condition_2))
  is_assert(setequal(needed,unique(meta$condition)),'Every sample condition must appear in a requested comparison')
  is_assert(all(table(meta$condition)>=2),'At least two independent biological replicates per condition are required')
  design <- data.frame(sampleID=meta$sample_id,condition=factor(meta$condition,levels=needed))
  cv <- is_default(cfg$covariates,list())
  for(k in names(cv)) {
    is_assert(k %in% names(meta) && make.names(k)==k && !k %in% c('sampleID','condition','sample','exon','quant','sample_id'),paste('Invalid covariate',k))
    is_assert(cv[[k]] %in% c('factor','numeric'),paste('Covariate type must be factor or numeric:',k))
    v <- meta[[k]]
    is_assert(!anyNA(v) && all(nzchar(v)),paste('Missing covariate:',k))
    design[[k]] <- if(cv[[k]]=='factor') factor(v) else suppressWarnings(as.numeric(v))
    is_assert(!anyNA(design[[k]]) && (cv[[k]]!='numeric' || all(is.finite(design[[k]]))),paste('Invalid covariate values:',k))
  }
  for(i in seq_len(nrow(comparisons))) {
    d <- droplevels(design[design$condition %in% unlist(comparisons[i,2:3]),,drop=FALSE])
    for(k in names(cv)) {
      is_assert(length(unique(d[[k]]))>1,paste('Covariate constant in comparison:',k,comparisons$name[i]))
      if(cv[[k]]=='numeric') is_assert(length(unique(d[[k]]))*2>nrow(d),
        paste('IsoformSwitchAnalyzeR would recast this numeric covariate as categorical:',k,'; explicitly choose factor or use the general DEJU/DGE model'))
    }
    mm <- model.matrix(reformulate(c('condition',names(cv))),d)
    is_assert(qr(mm)$rank==ncol(mm) && nrow(mm)>ncol(mm),paste('Confounded/saturated comparison design:',comparisons$name[i]))
  }
  list(meta=meta,design=design,comparisons=comparisons)
}
# Keep import arguments shared between preflight and execution. The release
# importRdata API does not have the development-only autoCastDesignCol argument.
is_import_args <- function(cfg,design=NULL,comparisons=NULL,counts=NULL,abundance=NULL) {
  list(isoformCountMatrix=counts,isoformRepExpression=abundance,
    designMatrix=design,isoformExonAnnoation=cfg$gtf,isoformNtFasta=cfg$transcript_fasta,
    comparisonsToMake=comparisons,detectUnwantedEffects=FALSE,
    # Reference Salmon quantifications must retain the supplied GTF gene IDs.
    # StringTie repair can replace distinct IDs with a shared gene symbol.
    fixStringTieAnnotationProblem=FALSE,
    addAnnotatedORFs=cfg$consequences,removeNonConvensionalChr=FALSE,removeTECgenes=FALSE,
    ignoreAfterBar=cfg$strip_pipe,ignoreAfterSpace=TRUE,ignoreAfterPeriod=FALSE,
    ignoreSurplusIsoforms=FALSE,estimateDifferentialGeneRange=FALSE)
}
is_import_gene_audit <- function(sw,annotation) {
  original <- unique(annotation[,c('isoform_id','gene_id')])
  actual <- unique(sw$isoformFeatures[,c('isoform_id','gene_id')])
  is_assert(!anyDuplicated(original$isoform_id),'Ambiguous original transcript-to-gene mapping')
  is_assert(!anyDuplicated(actual$isoform_id),'Imported transcript assigned to multiple genes')
  is_assert(all(actual$isoform_id %in% original$isoform_id),'Imported transcript absent from original mapping')
  names(original)[2] <- 'original_gene_id'
  names(actual)[2] <- 'imported_gene_id'
  audit <- merge(original,actual,by='isoform_id',all.x=TRUE,sort=FALSE)
  audit$status <- ifelse(is.na(audit$imported_gene_id),'not_imported',
    ifelse(audit$original_gene_id==audit$imported_gene_id,'preserved','changed'))
  audit
}
is_check_api <- function(name,arg_names,fn=getExportedValue('IsoformSwitchAnalyzeR',name)) {
  missing <- setdiff(arg_names,names(formals(fn)))
  is_assert(!length(missing),paste('Unsupported IsoformSwitchAnalyzeR API for',name,
    '; wrapper arguments absent from the installed function:',paste(missing,collapse=', '),
    '. Check wrapper/package compatibility; upgrading alone may not fix this.'))
  invisible(TRUE)
}
is_api <- function(name,args) {
  fn <- getExportedValue('IsoformSwitchAnalyzeR',name)
  is_check_api(name,names(args),fn)
  do.call(fn,args)
}
is_prepare_runtime <- function() {
  is_assert(requireNamespace('dplyr',quietly=TRUE),'Install dplyr for isoform filtering')
  required <- c('rename_with','inner_join','pull','rowwise','across','do','ungroup')
  is_assert(all(required %in% getNamespaceExports('dplyr')),
            'Installed dplyr lacks functions required by IsoformSwitchAnalyzeR filtering')
  # The 2.12 release calls these verbs without importing them in NAMESPACE.
  # Attach the dependency normally; do not alter installed package namespaces.
  if(!'package:dplyr' %in% search()) suppressPackageStartupMessages(library('dplyr',character.only=TRUE))
  invisible(TRUE)
}
is_dependencies <- function(cfg) {
  pkgs <- c('tximport','IsoformSwitchAnalyzeR','DEXSeq','rtracklayer','Biostrings','jsonlite','dplyr')
  if(!is.null(cfg$go)) pkgs <- c(pkgs,'clusterProfiler','AnnotationDbi',cfg$go$orgdb)
  missing <- pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
  is_assert(!length(missing),paste('Install required packages:',paste(missing,collapse=', ')))
  is_prepare_runtime()
  is_check_api('importRdata',names(is_import_args(cfg)))
  # Explicitly require the count-based filter API; never reinterpret old TPM cutoffs.
  is_assert(all(c('isoCount','min.Count.prop','min.IF.prop') %in% names(formals(IsoformSwitchAnalyzeR::preFilter))),
    'This module requires the current count-based IsoformSwitchAnalyzeR preFilter API (documented in 2.12); older releases need a coherent R/Bioconductor upgrade')
  invisible(pkgs)
}
is_quant <- function(meta,strip_pipe) {
  reference <- NULL
  for(p in meta$quant) {
    x <- utils::read.delim(p,check.names=FALSE,stringsAsFactors=FALSE)
    is_assert(all(c('Name','Length','EffectiveLength','TPM','NumReads') %in% names(x)),'Invalid Salmon quant.sf schema')
    ids <- is_id(as.character(x$Name),strip_pipe)
    for(k in c('Length','EffectiveLength','TPM','NumReads')) is_assert(is.numeric(x[[k]]) && all(is.finite(x[[k]])) && all(x[[k]]>=0),paste('Invalid quantification values:',p,k))
    is_assert(all(x$Length>0) && all(x$EffectiveLength>0) && sum(x$NumReads)>0,'Zero length/library quantification')
    if (is.null(reference)) reference <- ids else is_assert(identical(ids,reference),'Salmon transcript IDs/order differ between samples; use one reference index')
  }
  txi <- tximport::tximport(setNames(meta$quant,meta$sample_id),type='salmon',txOut=TRUE,
    countsFromAbundance='scaledTPM',dropInfReps=TRUE)
  is_assert(identical(colnames(txi$counts),meta$sample_id),'Imported samples are out of order')
  for(k in c('counts','abundance','length')) rownames(txi[[k]]) <- is_id(rownames(txi[[k]]),strip_pipe)
  txi
}
is_annotation <- function(cfg,ids) {
  anno <- rtracklayer::import(cfg$gtf)
  is_assert(all(c('type','transcript_id','gene_id') %in% names(S4Vectors::mcols(anno))),'GTF lacks transcript/gene IDs')
  ex <- anno[anno$type=='exon']
  pairs <- unique(data.frame(isoform_id=as.character(ex$transcript_id),gene_id=as.character(ex$gene_id)))
  is_assert(nrow(pairs)>0 && !anyNA(pairs) && all(nzchar(pairs$isoform_id)) && all(nzchar(pairs$gene_id)),'Missing exon transcript/gene identifiers')
  is_assert(!anyDuplicated(pairs$isoform_id),'One transcript ID maps to multiple genes')
  is_assert(all(ids %in% pairs$isoform_id),'Salmon transcript IDs absent from GTF: check annotation release, versions and pipe policy')
  if(!is.null(cfg$transcript_fasta)) {
    seqs <- Biostrings::readDNAStringSet(cfg$transcript_fasta)
    names(seqs) <- is_id(sub('[[:space:]].*$', '',names(seqs)),cfg$strip_pipe)
    is_assert(all(ids %in% names(seqs)),'Transcript FASTA is missing quantified transcripts')
    # Exact exon-spliced length catches a common wrong FASTA/GTF pairing.
    widths <- tapply(IRanges::width(ex),as.character(ex$transcript_id),sum)
    is_assert(all(BiocGenerics::width(seqs[match(ids,names(seqs))])==widths[ids]),'Transcript FASTA lengths differ from GTF exon lengths')
  }
  pairs[match(ids,pairs$isoform_id),,drop=FALSE]
}
is_results <- function(tab,comparisons,alpha,delta) {
  need <- c('isoform_id','gene_id','condition_1','condition_2','IF1','IF2','dIF','isoform_switch_q_value','gene_switch_q_value')
  is_assert(all(need %in% names(tab)),'Unexpected IsoformSwitchAnalyzeR result schema')
  is_assert(all(is.na(tab$dIF) | abs(tab$dIF-(tab$IF2-tab$IF1))<1e-7),'Unexpected dIF direction: expected IF2 minus IF1')
  key <- paste(tab$condition_1,tab$condition_2,sep='\r')
  idx <- match(key,paste(comparisons$condition_1,comparisons$condition_2,sep='\r'))
  is_assert(!anyNA(idx),'Unexpected comparison in result')
  tab$comparison <- comparisons$name[idx]
  is_assert(!anyDuplicated(tab[,c('comparison','isoform_id')]),'Duplicate isoform/comparison result')
  tab$passes_fdr <- !is.na(tab$isoform_switch_q_value) & tab$isoform_switch_q_value<alpha
  tab$passes_effect <- !is.na(tab$dIF) & abs(tab$dIF)>delta
  tab$significant <- tab$passes_fdr & tab$passes_effect
  summaries <- lapply(split(seq_len(nrow(tab)),paste(tab$comparison,tab$gene_id,sep='\r')),function(i) {
    t <- tab[i,,drop=FALSE];q <- unique(t$gene_switch_q_value[!is.na(t$gene_switch_q_value)])
    is_assert(length(q)<=1,'Inconsistent gene q-values within comparison')
    data.frame(comparison=t$comparison[1],gene_id=t$gene_id[1],tested_isoform_rows=sum(!is.na(t$isoform_switch_q_value)),
      significant_isoform_rows=sum(t$significant),gene_switch_q_value=if(length(q)) q else NA_real_,
      significant_up=sum(t$significant & t$dIF>0),significant_down=sum(t$significant & t$dIF<0))
  })
  list(all=tab,significant=tab[tab$significant,,drop=FALSE],genes=do.call(rbind,summaries))
}
is_go <- function(tab,cfg,out) {
  opt <- cfg$go
  mapping <- utils::read.delim(opt$mapping,colClasses='character',check.names=FALSE)
  is_assert(all(c('gene_id','entrez_id') %in% names(mapping)) && !anyNA(mapping[,c('gene_id','entrez_id')]),'GO mapping needs gene_id and entrez_id')
  mapping <- unique(mapping[,c('gene_id','entrez_id')])
  is_assert(all(nzchar(mapping$gene_id)) && all(nzchar(mapping$entrez_id)),'Blank GO mapping IDs')
  ambiguous <- unique(mapping$gene_id[duplicated(mapping$gene_id) | duplicated(mapping$gene_id,fromLast=TRUE)])
  valid <- mapping[!mapping$gene_id %in% ambiguous,,drop=FALSE]
  tested <- unique(tab$gene_id[!is.na(tab$isoform_switch_q_value)])
  audit <- data.frame(gene_id=tested,status=ifelse(tested %in% ambiguous,'ambiguous_excluded',ifelse(tested %in% valid$gene_id,'mapped','unmapped')))
  is_write(audit,file.path(out,'go_mapping_audit.tsv'))
  db <- getExportedValue(opt$orgdb,opt$orgdb)
  for(cmp in unique(tab$comparison)) {
    t <- tab[tab$comparison==cmp,,drop=FALSE]
    universe <- unique(valid$entrez_id[valid$gene_id %in% t$gene_id[!is.na(t$isoform_switch_q_value)]])
    selected <- unique(valid$entrez_id[valid$gene_id %in% t$gene_id[t$significant]])
    is_write(data.frame(entrez_id=universe),file.path(out,paste0(cmp,'.go_universe.tsv')))
    if(!length(selected) || !length(universe)) {
      is_write(data.frame(status='No mapped significant/tested genes'),file.path(out,paste0(cmp,'.go_status.tsv')))
      next
    }
    result <- clusterProfiler::enrichGO(gene=selected,universe=universe,OrgDb=db,keyType='ENTREZID',ont=opt$ontology,
      pAdjustMethod='BH',pvalueCutoff=1,qvalueCutoff=1,readable=FALSE)
    is_write(as.data.frame(result),file.path(out,paste0(cmp,'.go.tsv')))
  }
}

# Plot retained exon structures and group isoform fractions without ORF assumptions.
# This path does not fit a model, annotate coding potential, or modify the object.
is_plot_core_switches <- function(sw,cfg,out) {
  comparisons <- do.call(rbind,lapply(cfg$comparisons,function(x)
    data.frame(name=x$name,condition_1=x$reference,condition_2=x$treatment)))
  results <- is_results(sw$isoformFeatures,comparisons,cfg$alpha,cfg$delta_if)
  is_assert(!dir.exists(out),'Plot output already exists; choose a new directory')
  dir.create(out,recursive=TRUE)
  writeLines('Plot generation incomplete.',file.path(out,'INCOMPLETE.txt'))
  ex <- as.data.frame(sw$exons)
  is_assert(all(c('isoform_id','seqnames','start','end','strand') %in% names(ex)),
            'Saved object lacks exon coordinates required for structure plots')
  index <- data.frame(comparison=character(),gene_id=character(),file=character(),pages=integer())
  plotted <- list()
  for (cmp in comparisons$name) {
    significant <- results$significant[results$significant$comparison==cmp,,drop=FALSE]
    if (!nrow(significant) || cfg$plots==0) next
    # Rank by minimum isoform q-value for display only, not gene-level inference.
    scores <- tapply(significant$isoform_switch_q_value,significant$gene_id,min)
    genes <- names(scores)[order(scores,names(scores))]
    genes <- head(genes,cfg$plots)
    for (gene in genes) {
      tab <- results$all[results$all$comparison==cmp & results$all$gene_id==gene,,drop=FALSE]
      tab <- tab[order(-pmax(tab$IF1,tab$IF2),tab$isoform_id),,drop=FALSE]
      coords <- ex[ex$isoform_id %in% tab$isoform_id,,drop=FALSE]
      is_assert(all(tab$isoform_id %in% coords$isoform_id),'Missing exon structures for plotted isoforms')
      chr <- unique(as.character(coords$seqnames))
      is_assert(length(chr)==1,paste('Cannot plot a gene on multiple sequences:',gene))
      is_assert(all(is.finite(tab$IF1)&is.finite(tab$IF2)&tab$IF1>=0&tab$IF1<=1&tab$IF2>=0&tab$IF2<=1),
                'Invalid isoform fractions for plotting')
      filename <- sprintf('%03d_%s_%s.pdf',nrow(index)+1L,
                          gsub('[^A-Za-z0-9_.-]','_',cmp),gsub('[^A-Za-z0-9_.-]','_',gene))
      pages <- split(seq_len(nrow(tab)),ceiling(seq_len(nrow(tab))/12))
      grDevices::pdf(file.path(out,filename),width=13,height=8,onefile=TRUE)
      tryCatch({
        for (ii in pages) {
          t <- tab[ii,,drop=FALSE]; y <- rev(seq_len(nrow(t)))
          graphics::layout(matrix(c(1,2),nrow=1),widths=c(2.1,1))
          graphics::par(mar=c(6,12,5,1))
          limits <- range(c(coords$start,coords$end))
          graphics::plot(NA,xlim=limits,ylim=c(.4,nrow(t)+.6),yaxt='n',
            xlab=paste(chr,'genomic position (increasing left to right)'),ylab='',
            main=paste(gene,cmp,sep=' | '))
          graphics::axis(2,at=y,labels=t$isoform_id,las=2,cex.axis=.65)
          for (k in seq_len(nrow(t))) {
            e <- coords[coords$isoform_id==t$isoform_id[k],,drop=FALSE]
            graphics::segments(min(e$start),y[k],max(e$end),y[k],col='grey60')
            graphics::rect(e$start,y[k]-.18,e$end,y[k]+.18,col='#376795',border=NA)
          }
          graphics::mtext(paste('Exons only; coding status not displayed. Strand:',
            paste(unique(as.character(coords$strand)),collapse=',')),side=3,line=.2,cex=.7)
          graphics::par(mar=c(6,1,5,1))
          graphics::plot(NA,xlim=c(0,1),ylim=c(.4,nrow(t)+.6),yaxt='n',
            xlab='Group isoform fraction',ylab='',main='Relative isoform usage')
          graphics::segments(t$IF1,y,t$IF2,y,col='grey55')
          graphics::points(t$IF1,y,pch=16,col='#2166AC')
          graphics::points(t$IF2,y,pch=17,col='#B35806')
          graphics::legend('top',inset=c(0,-.16),xpd=NA,bty='n',horiz=TRUE,
            legend=c(t$condition_1[1],t$condition_2[1]),pch=c(16,17),col=c('#2166AC','#B35806'),cex=.8)
          graphics::mtext('Group estimates; no confidence intervals shown',side=1,line=4.3,cex=.65)
        }
      },finally=grDevices::dev.off())
      index <- rbind(index,data.frame(comparison=cmp,gene_id=gene,file=filename,pages=length(pages)))
      plotted[[length(plotted)+1L]] <- tab
    }
  }
  is_write(index,file.path(out,'plot_index.tsv'))
  is_write(if(length(plotted)) do.call(rbind,plotted) else results$all[FALSE,,drop=FALSE],
           file.path(out,'plotted_isoform_statistics.tsv'))
  writeLines(c('Plots show retained/tested isoforms, not all reference transcripts.',
    'Uniform-height exon boxes indicate exon structure only, not noncoding status.',
    'Coordinates are genomic and not intron-rescaled; distant exons may appear small.',
    'Selection: up to configured plots per comparison, ranked by minimum significant isoform q-value.',
    'Ranking is descriptive, not a calibrated gene-level FDR. No statistics were refitted.',
    'Review original annotation exclusions and QC before biological interpretation.'),file.path(out,'README.txt'))
  writeLines('Requested core structure/fraction plots completed; original workflow status is unchanged.',
             file.path(out,'PLOTS_COMPLETE.txt'))
  unlink(file.path(out,'INCOMPLETE.txt'))
  invisible(index)
}

is_plot_switches <- function(sw,cfg,out) {
  if (!is.data.frame(sw$orfAnalysis) || !nrow(sw$orfAnalysis)) {
    message('No ORF annotation: plotting exon structures and isoform fractions without coding inference.')
    return(is_plot_core_switches(sw,cfg,out))
  }
  dir.create(out,recursive=TRUE)
  is_api('switchPlotTopSwitches',list(switchAnalyzeRlist=sw,alpha=cfg$alpha,dIFcutoff=cfg$delta_if,n=cfg$plots,
    filterForConsequences=FALSE,pathToOutput=out,splitComparison=TRUE,splitFunctionalConsequences=FALSE))
}
