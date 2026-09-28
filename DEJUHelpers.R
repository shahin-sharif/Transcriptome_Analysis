# DEJU annotation, sample validation and count-table helpers.
# Sourcing this file defines functions only; it does not run an analysis.
deju_assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}
deju_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  deju_assert(!length(missing), paste('Missing packages:', paste(missing, collapse = ', ')))
}
deju_write <- function(x, path) {
  con <- if (grepl('[.]gz$', path)) gzfile(path, 'wt') else file(path, 'wt')
  on.exit(close(con))
  utils::write.table(x, con, sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA')
}
deju_tokens <- function(x) {
  if (is.na(x) || !nzchar(x) || x %in% c('NA', '.', '-')) return(character())
  unique(trimws(strsplit(x, '[,;|]')[[1]]))
}
deju_samples <- function(path, reference, treatment, covariates = character(), require_bams = TRUE) {
  deju_assert(file.exists(path), 'Sample sheet does not exist')
  s <- utils::read.delim(path, colClasses = 'character', check.names = FALSE,
                        na.strings = '', comment.char = '', quote = '')
  need <- c('sample_id', 'condition', covariates, if (require_bams) 'bam')
  deju_assert(all(need %in% names(s)), paste('Sample sheet requires:', paste(need, collapse = ', ')))
  deju_assert(nrow(s) > 0 && !anyNA(s[, need, drop = FALSE]), 'Missing sample metadata')
  deju_assert(!anyDuplicated(s$sample_id) && all(grepl('^[A-Za-z0-9][A-Za-z0-9_.-]*$', s$sample_id)),
              'Sample IDs must be unique and use letters, digits, dot, dash or underscore')
  deju_assert(reference != treatment && setequal(unique(s$condition), c(reference, treatment)),
              'Sample sheet must contain exactly the selected reference and treatment conditions')
  deju_assert(all(table(s$condition) >= 2), 'At least two biological replicates per condition are required')
  if (require_bams) {
    abs <- grepl('^(/|~)', s$bam)
    s$bam[!abs] <- file.path(dirname(normalizePath(path)), s$bam[!abs])
    s$bam <- normalizePath(path.expand(s$bam), mustWork = TRUE)
    deju_assert(!anyDuplicated(s$bam), 'The same BAM cannot represent multiple biological samples')
    deju_assert(all(grepl('[.]bam$', s$bam, ignore.case = TRUE)), 'Provide BAM files, not SAM or CRAM')
  }
  s
}
deju_design <- function(samples, reference, treatment, covariates = character()) {
  deju_assert(!anyDuplicated(covariates) && !any(covariates %in% c('condition', 'sample_id', 'bam')),
              'Covariates must be distinct additional sample-sheet columns')
  deju_assert(all(make.names(covariates) == covariates), 'Covariate names must be valid R names')
  s <- samples
  s$condition <- factor(s$condition, levels = c(reference, treatment))
  for (v in covariates) {
    deju_assert(v %in% names(s) && !anyNA(s[[v]]) && all(nzchar(as.character(s[[v]]))), paste('Missing covariate:', v))
    # Categorical batch/subject effects; no implicit numeric encoding.
    s[[v]] <- factor(s[[v]])
    deju_assert(nlevels(s[[v]]) > 1, paste('Constant covariate:', v))
  }
  design <- stats::model.matrix(stats::reformulate(c(covariates, 'condition')), s)
  rownames(design) <- s$sample_id
  deju_assert(qr(design)$rank == ncol(design), 'Design is not full rank: condition may be confounded with batch/subject')
  deju_assert(nrow(design) > ncol(design), 'Design has no residual degrees of freedom')
  # Condition is the final additive term and has exactly two levels.
  coefficient <- which(attr(design, 'assign') == length(covariates) + 1L)
  deju_assert(length(coefficient) == 1L, 'Could not identify treatment-minus-reference coefficient')
  list(matrix = design, coef = coefficient)
}

# Read only exon records in chunks; keep full versioned IDs for internal matching.
deju_read_exons <- function(gtf) {
  deju_assert(file.exists(gtf), 'GTF file does not exist')
  con <- if (grepl('[.]gz$', gtf)) gzfile(gtf, 'rt') else file(gtf, 'rt')
  on.exit(close(con))
  chunks <- list()
  attr_value <- function(x, key) {
    pat <- paste0('(?:^|;)[[:space:]]*', key, '[[:space:]]+"([^"]+)"')
    hit <- regexec(pat, x, perl = TRUE)
    vapply(regmatches(x, hit), function(z) if (length(z) >= 2) z[2] else NA_character_, character(1))
  }
  repeat {
    z <- readLines(con, n = 50000L, warn = FALSE)
    if (!length(z)) break
    z <- z[!startsWith(z, '#') & grepl('^[^\t]+\t[^\t]+\texon\t', z)]
    if (!length(z)) next
    fields <- strsplit(z, '\t', fixed = TRUE)
    deju_assert(all(lengths(fields) == 9L), 'Malformed exon record in GTF')
    m <- do.call(rbind, fields)
    d <- data.frame(gene_id = attr_value(m[,9], 'gene_id'), transcript_id = attr_value(m[,9], 'transcript_id'),
                    chr = m[,1], start = suppressWarnings(as.integer(m[,4])), end = suppressWarnings(as.integer(m[,5])),
                    strand = m[,7], stringsAsFactors = FALSE)
    deju_assert(!anyNA(d) && all(d$start > 0 & d$end >= d$start) && all(d$strand %in% c('+','-')),
                'Exons require valid coordinates, strand, gene_id and transcript_id')
    chunks[[length(chunks)+1L]] <- d
  }
  deju_assert(length(chunks) > 0L, 'No GTF exon records found')
  unique(do.call(rbind, chunks))
}
prepare_deju_annotation <- function(gtf) {
  deju_require(c('GenomicRanges', 'IRanges'))
  ex <- deju_read_exons(gtf)
  # One locus per gene ID; never silently combine chromosomes/strands.
  locus <- unique(ex[,c('gene_id','chr','strand')])
  deju_assert(!anyDuplicated(locus$gene_id), 'A gene ID occurs on multiple chromosomes/strands; disambiguate the annotation')
  tx_locus <- unique(ex[,c('transcript_id','gene_id','chr','strand')])
  deju_assert(!anyDuplicated(tx_locus$transcript_id), 'A transcript ID has conflicting gene/locus annotations')
  gr <- GenomicRanges::GRanges(ex$chr, IRanges::IRanges(ex$start, ex$end), strand = ex$strand)
  by_gene <- GenomicRanges::GRangesList(split(gr, ex$gene_id))
  merged <- GenomicRanges::reduce(by_gene)
  saf <- as.data.frame(unlist(merged, use.names = FALSE))
  saf <- data.frame(GeneID = rep(names(merged), lengths(merged)), Chr = as.character(saf$seqnames),
                    Start = saf$start, End = saf$end, Strand = as.character(saf$strand), stringsAsFactors = FALSE)
  saf <- saf[order(saf$Chr, saf$Start, saf$End, saf$GeneID), ]; rownames(saf) <- NULL
  genes <- do.call(rbind,lapply(split(ex,ex$gene_id),function(d)
    data.frame(gene_id=d$gene_id[1],chr=d$chr[1],start=min(d$start),end=max(d$end),strand=d$strand[1])))
  rownames(genes) <- NULL
  empty <- data.frame(gene_id=character(), transcript_id=character(), chr=character(), left=integer(), right=integer(),
                       strand=character(), donor=integer(), acceptor=integer(), junction_rank=integer(),
                       n_junctions=integer(), scaled_position=numeric())
  rows <- lapply(split(ex, ex$transcript_id), function(d) {
    d <- d[order(d$start, d$end), ]
    if (nrow(d) < 2L) return(NULL)
    deju_assert(all(head(d$end,-1L) < tail(d$start,-1L)), paste('Overlapping exons within transcript', d$transcript_id[1]))
    left <- head(d$end,-1L); right <- tail(d$start,-1L)
    keep <- right > left + 1L; left <- left[keep]; right <- right[keep]
    n <- length(left); if (!n) return(NULL)
    plus <- d$strand[1] == '+'; rank <- if (plus) seq_len(n) else rev(seq_len(n))
    data.frame(gene_id=d$gene_id[1], transcript_id=d$transcript_id[1], chr=d$chr[1], left=left, right=right,
               strand=d$strand[1], donor=if (plus) left else right, acceptor=if (plus) right else left,
               junction_rank=rank, n_junctions=n, scaled_position=(rank-0.5)/n, stringsAsFactors=FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  reference <- if (length(rows)) do.call(rbind, rows) else empty
  rownames(reference) <- NULL
  list(saf=saf, genes=genes, junctions=reference, gtf_md5=unname(tools::md5sum(gtf)))
}

deju_sample_columns <- function(columns, targets, sample_ids) {
  deju_assert(length(targets) == length(sample_ids), 'Target/sample length mismatch')
  idx <- match(targets, columns)
  if (anyNA(idx)) {
    bn <- basename(columns); wanted <- basename(targets)
    deju_assert(!anyDuplicated(wanted) && all(vapply(wanted, function(z) sum(bn==z)==1L, logical(1))),
                'Junction sample columns cannot be matched uniquely to BAM targets')
    idx <- match(wanted, bn)
  }
  deju_assert(!anyNA(idx) && !anyDuplicated(idx), 'Ambiguous junction sample columns')
  idx
}
deju_standardize_junctions <- function(d) {
  pick <- function(aliases, required=TRUE) {
    hits <- match(tolower(aliases), tolower(names(d))); hits <- hits[!is.na(hits)]
    if (!length(hits)) {
      deju_assert(!required, paste('Missing junction annotation; expected one of', paste(aliases,collapse=', ')))
      return(rep(NA_character_, nrow(d)))
    }
    as.character(d[[hits[1]]])
  }
  chr1 <- pick(c('Site1_chr','Chr_SP1','Chr1'))
  chr2 <- pick(c('Site2_chr','Chr_SP2','Chr2'))
  p1 <- suppressWarnings(as.integer(pick(c('Site1_location','Location_SP1','SP1'))))
  p2 <- suppressWarnings(as.integer(pick(c('Site2_location','Location_SP2','SP2'))))
  data.frame(chr=chr1, chr2=chr2, left=pmin(p1,p2), right=pmax(p1,p2),
              primary=pick(c('PrimaryGene','PrimaryGeneID','Gene_SP1'),FALSE),
              secondary=pick(c('SecondaryGenes','SecondaryGene','SecondaryGeneID','Gene_SP2'),FALSE), stringsAsFactors=FALSE)
}
deju_assign_junctions <- function(j, annotation) {
  deju_require(c('GenomicRanges','IRanges'))
  j$gene_id <- NA_character_; j$strand <- NA_character_; j$assignment <- 'invalid_coordinates'
  j$candidate_genes <- ''; j$annotated <- FALSE
  valid <- !is.na(j$chr) & !is.na(j$chr2) & nzchar(j$chr) & j$chr==j$chr2 &
    !is.na(j$left) & !is.na(j$right) & j$left>0L & j$right>j$left+1L
  key <- function(chr,left,right) paste(chr,left,right,sep=':')
  ref <- annotation$junctions
  if (!'read_strand' %in% names(j)) j$read_strand <- '*'
  reference <- split(ref$gene_id, key(ref$chr,ref$left,ref$right))
  genes <- annotation$genes
  region <- GenomicRanges::GRanges(genes$chr, IRanges::IRanges(genes$start,genes$end))
  ids <- which(valid); candidates <- vector('list', nrow(j))
  if (length(ids)) {
    q <- GenomicRanges::GRanges(j$chr[ids], IRanges::IRanges(j$left[ids],j$right[ids]))
    ov <- GenomicRanges::findOverlaps(q,region,type='within',ignore.strand=TRUE)
    hit <- split(genes$gene_id[S4Vectors::subjectHits(ov)], ids[S4Vectors::queryHits(ov)])
    for (k in names(hit)) candidates[[as.integer(k)]] <- unique(hit[[k]])
  }
  for (i in ids) {
    exact <- unique(reference[[key(j$chr[i],j$left[i],j$right[i])]])
    if (j$read_strand[i] != '*') {
      compatible <- genes$gene_id[genes$strand==j$read_strand[i]]
      exact <- intersect(exact,compatible)
      candidates[[i]] <- intersect(candidates[[i]],compatible)
    }
    if (length(exact)) {
      possible <- exact; status <- 'annotated'; j$annotated[i] <- TRUE
    } else {
      # Novel junction: unique containing gene with concordant available featureCounts hints.
      possible <- candidates[[i]]
      hints <- union(deju_tokens(j$primary[i]),deju_tokens(j$secondary[i]))
      status <- 'novel_intragenic'
      if (length(possible)==1L && length(hints) && !possible %in% hints) {
        j$assignment[i] <- 'conflicting_gene_hint'; j$candidate_genes[i] <- paste(possible,collapse=';'); next
      }
    }
    j$candidate_genes[i] <- paste(sort(possible),collapse=';')
    if (length(possible)!=1L) {
      j$assignment[i] <- if (length(possible)) 'ambiguous_gene' else 'unassigned'; next
    }
    j$gene_id[i] <- possible; j$strand[i] <- genes$strand[match(possible,genes$gene_id)]
    j$assignment[i] <- status
  }
  j$donor <- ifelse(j$strand=='+',j$left,j$right)
  j$acceptor <- ifelse(j$strand=='+',j$right,j$left)
  j
}

deju_validate_counts <- function(counts, features, samples) {
  counts <- as.matrix(counts)
  deju_assert(is.numeric(counts) && all(is.finite(counts)) && all(counts>=0) && all(counts==round(counts)),
              'Counts must be finite nonnegative integers, not CPM/TPM or transformed values')
  need <- c('feature_id','gene_id','feature_type')
  deju_assert(all(need %in% names(features)), 'Feature table requires feature_id, gene_id, feature_type')
  deju_assert(!anyNA(features[,need]) && all(nzchar(features$feature_id)) && all(nzchar(features$gene_id)) &&
                !anyDuplicated(features$feature_id), 'Feature IDs must be unique; gene/feature IDs cannot be missing')
  deju_assert(all(features$feature_type %in% c('exon','junction')), 'Unknown feature_type')
  deju_assert(!is.null(rownames(counts)) && !anyDuplicated(rownames(counts)) &&
                setequal(rownames(counts),features$feature_id), 'Count rows and feature IDs differ')
  deju_assert(!is.null(colnames(counts)) && !anyDuplicated(colnames(counts)) &&
                setequal(colnames(counts),samples$sample_id), 'Count columns and sample IDs differ')
  counts[features$feature_id,samples$sample_id,drop=FALSE]
}

# Count N CIGAR operations directly so strand and filtering apply equally to
# junctions and exon reads. One alignment contributes once per crossed intron.
deju_junction_counts <- function(files,sample_ids,paired,strand,yield_size=250000L,min_anchor=8L,min_intron=20L,max_intron=1000000L) {
  deju_require(c('Rsamtools','GenomicAlignments'))
  deju_assert(min_anchor>=0 && min_intron>=1 && max_intron>=min_intron,'Invalid junction filters')
  one <- function(path) {
    bf <- Rsamtools::BamFile(path,yieldSize=yield_size); open(bf); on.exit(close(bf))
    totals <- new.env(hash=TRUE,parent=emptyenv())
    repeat {
      b <- Rsamtools::scanBam(bf,param=Rsamtools::ScanBamParam(what=c('rname','pos','cigar','flag')))[[1]]
      if (!length(b$pos)) break
      is_pair <- bitwAnd(b$flag,1L)!=0L
      deju_assert(all(is_pair==paired),'--paired conflicts with BAM flags; mixed single/paired libraries are not supported')
      first <- bitwAnd(b$flag,64L)!=0L; second <- bitwAnd(b$flag,128L)!=0L
      if (paired) deju_assert(all(xor(first,second)),'Paired alignments require exactly one first/second-mate flag')
      take <- grepl('N',b$cigar,fixed=TRUE)
      if (!any(take)) next
      b <- lapply(b,function(x) x[take]); second <- second[take]
      introns <- GenomicAlignments::cigarRangesAlongReferenceSpace(b$cigar,pos=b$pos,ops='N')
      n <- lengths(introns); introns <- unlist(introns,use.names=FALSE)
      negative <- bitwAnd(b$flag,16L)!=0L
      if (paired) negative <- xor(negative,second)
      if (strand==2L) negative <- !negative
      orientation <- if (strand==0L) rep('*',length(negative)) else ifelse(negative,'-','+')
      key <- paste(rep(as.character(b$rname),n),IRanges::start(introns)-1L,IRanges::end(introns)+1L,rep(orientation,n),sep='\t')
      ops <- GenomicAlignments::explodeCigarOps(b$cigar)
      widths <- GenomicAlignments::explodeCigarOpLengths(b$cigar)
      anchors <- unlist(Map(function(op,w) {
        at<-which(op=='N')
        vapply(at,function(k) {
          left<-0;right<-0;j<-k-1L
          while(j>=1L && op[j] %in% c('M','=','X')) {left<-left+w[j];j<-j-1L}
          j<-k+1L
          while(j<=length(op) && op[j] %in% c('M','=','X')) {right<-right+w[j];j<-j+1L}
          min(left,right)>=min_anchor
        },logical(1))
      },ops,widths),use.names=FALSE)
      keep<-anchors & IRanges::width(introns)>=min_intron & IRanges::width(introns)<=max_intron
      tab <- table(key[keep])
      for (k in names(tab)) {
        prior <- totals[[k]]; if (is.null(prior)) prior<-0
        totals[[k]] <- prior+as.numeric(tab[[k]])
      }
    }
    keys <- ls(totals,all.names=TRUE)
    setNames(vapply(keys,function(k) totals[[k]],numeric(1)),keys)
  }
  per_sample <- lapply(files,one)
  keys <- sort(unique(unlist(lapply(per_sample,names),use.names=FALSE)))
  deju_assert(length(keys)>0,'No junctions remain after BAM filtering; verify NH tags, MAPQ and CIGAR strings')
  values <- matrix(0,length(keys),length(files),dimnames=list(NULL,sample_ids))
  for (i in seq_along(files)) values[match(names(per_sample[[i]]),keys),i] <- per_sample[[i]]
  fields <- do.call(rbind,strsplit(keys,'\t',fixed=TRUE))
  annotation <- data.frame(chr=fields[,1],chr2=fields[,1],left=as.integer(fields[,2]),right=as.integer(fields[,3]),
                            read_strand=fields[,4],primary=NA_character_,secondary=NA_character_)
  list(counts=values,annotation=annotation)
}
