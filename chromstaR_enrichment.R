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
              help="Annotation file in BED6 format."),
  make_option(c("-b", "--bpAroundAnnotation"), type = 'numeric', default = 10000,
              help="Base-pairs to consider around annotation. Default is %default."),
  make_option(c("-i", "--numIntervals"), type = 'numeric', default = 20,
              help="Intervals to consider inside annotation. Default is %default."),
  make_option(c("-s", "--statistic"), type = 'character', default = 'fold',
              help="Statistic to calculate. Either 'fold' or 'fraction'. Default is %default."),
  make_option(c("-l", "--numLoci"), type = 'numeric', default = Inf,
              help="Maximum number of rows for the count heatmap. Default is %default."),
  make_option(c("-S", "--sortBySample"), type = 'numeric', default = NULL,
              help="Number of the sample by which the count heatmap is to be sorted. Default is %default.")
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
ptm <- chromstaR:::startTimedMessage("Loading model ...")
model <- loadHmmsFromFiles(opt$chromstarObject)[[1]]
chromstaR:::stopTimedMessage(ptm)
savefolder <- "plotEnrichment"
if (!file.exists(savefolder)) { dir.create(savefolder) }

## Plot enrichment around annotation
for (what in c('counts', 'peaks', 'combinations')) {
    ggplts <- plotEnrichment(hmm = model, annotation = annotation, bp.around.annotation = opt$bpAroundAnnotation, num.intervals = opt$numIntervals, what = what, statistic = opt$statistic)
    for (i1 in 1:length(ggplts)) {
      ggsave(ggplts[[i1]], filename = file.path(savefolder, paste0(what, '_', names(ggplts)[i1], '.png')), width = 10, height = 7)
    }
}

## Count heatmap around start of annotation
ggplt <- plotEnrichCountHeatmap(hmm = model, annotation = annotation, bp.around.annotation = opt$bpAroundAnnotation, max.rows = opt$numLoci, sortByColumns = opt$sortBySample)
ggsave(ggplt, filename = file.path(savefolder, paste0('countHeatmap.png')), width = 10, height = 7)

total.time <- proc.time() - ptm.start; message("Total elapsed time: ",round(total.time[3],2),"s")

