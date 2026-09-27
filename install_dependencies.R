#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
if (length(setdiff(args,'--go'))) stop('Usage: Rscript install_dependencies.R [--go]')
if (!requireNamespace('BiocManager',quietly=TRUE)) install.packages('BiocManager',repos='https://cloud.r-project.org')
pkgs <- c('edgeR','limma','statmod','Rsubread','Rsamtools','GenomicRanges','GenomicAlignments','IRanges','S4Vectors')
if ('--go' %in% args) pkgs<-c(pkgs,'clusterProfiler','org.Hs.eg.db')
BiocManager::install(pkgs,ask=FALSE,update=FALSE)
if (utils::packageVersion('edgeR') < '4.6.0') warning('Modern mode needs edgeR >= 4.6. Use a newer R/Bioconductor environment, or explicitly choose legacy with edgeR >= 4.0.')
