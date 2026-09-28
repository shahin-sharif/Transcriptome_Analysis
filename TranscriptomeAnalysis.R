#!/usr/bin/env Rscript
args<-commandArgs(TRUE)
if(length(args)<1) stop('Usage: Rscript TranscriptomeAnalysis.R resolved_config.json [--check]',call.=FALSE)
file_arg<-grep('^--file=',commandArgs(FALSE),value=TRUE)
root<-dirname(normalizePath(sub('^--file=','',file_arg[1])))
source(file.path(root,'DEJUPipeline.R'))
source(file.path(root,'modules','AnalysisCore.R'))
deju_require('jsonlite')
ta_main(args[1],check='--check' %in% args)
