#!/usr/bin/env Rscript

# Set up R error handling to go to stderr
options(show.error.messages=FALSE,
        error=function() {
          cat(geterrmessage(), file=stderr())
          q(save="no", status=1, runLast=FALSE)
          }
        )

# Avoid crashing Galaxy with an UTF8 error on German LC settings
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

## Parse input arguments
library(optparse)
option_list <- list(
  make_option(c("-o", "--chromstarObject"), type = 'character',
              help="File with chromstaR object."),
  make_option(c("-a", "--annotationBed6"), type = 'character',
              help="Annotation file in BED-6 format."),
  make_option(c("-b", "--bpAroundAnnotation"), type = 'numeric', default = 10000,
              help="Base-pairs to consider around annotation. Default is %default."),
  make_option(c("-i", "--numIntervals"), type = 'numeric', default = 20,
              help="Intervals to consider inside annotation. Default is %default.")
)
opt <- parse_args(OptionParser(option_list=option_list))

## Sanity checks
message("Starting chromstaR.R with options")
message(paste0(paste0(names(opt), " = ", opt), "\n"), appendLF = FALSE)

#=====================
### Load libraries ###
#=====================
message("Loading libraries ...", appendLF=FALSE); ptm.start <- proc.time()
suppressPackageStartupMessages(library(chromstaR))
time <- proc.time() - ptm.start; message(" ",round(time[3],2),"s")

#===========
### Main ###
#===========

annotation <- readCustomBedFile(bedfile = opt$annotationBed6, col.classes = c('character', 'numeric', 'numeric', 'character', 'numeric', 'character'), chromosome.format = NULL)
model <- loadHmmsFromFiles(opt$chromstarObject)[[1]]

## Plot enrichment around annotation
savefolder <- "plotEnrichment"
if (!file.exists(savefolder)) { dir.create(savefolder) }
for (what in c('counts', 'peaks', 'combinations', 'transitions')) {
    ggplts <- plotEnrichment(hmm = model, annotation = annotation, bp.around.annotation = opt$bpAroundAnnotation, num.intervals = opt$numIntervals, what = what)
    for (i1 in 1:length(ggplts)) {
      ggsave(ggplts[[i1]], filename = file.path(savefolder, paste0(what, '_', names(ggplts)[i1], '.png')), width = 10, height = 7)
    }
}

total.time <- proc.time() - ptm.start; message("Total elapsed time: ",round(total.time[3],2),"s")

