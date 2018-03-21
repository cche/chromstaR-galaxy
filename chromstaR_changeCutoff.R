#!/usr/bin/env Rscript

# Set up R error handling to go to stderr
options(show.error.messages=FALSE,
        error=function() {
          cat(geterrmessage(), file=stderr())
          q(save="no", status=1, runLast=FALSE)
          }
        )

## Parse input arguments
library(optparse)
option_list <- list(
  make_option(c("-o", "--chromstarObject"), type = 'character',
              help="File with chromstaR object."),
  make_option(c("-w", "--changeWhat"), type = 'character', default = "changeMaxPostCutoff",
              help="Function to use. Either 'changePostCutoff' or 'changeMaxPostCutoff'. Default is %default."),
  make_option(c("-c", "--cutoff"), type = 'numeric', default = 0.9999,
              help="Cutoff to apply on the posteriors or maxPostInPeak. Default is %default.")
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

ptm <- chromstaR:::startTimedMessage("Loading model ...")
model <- loadHmmsFromFiles(opt$chromstarObject)[[1]]
chromstaR:::stopTimedMessage(ptm)

if (opt$changeWhat == 'changePostCutoff') {
    new.model <- changePostCutoff(model, post.cutoff = opt$cutoff)
} else if (opt$changeWhat == 'changeMaxPostCutoff') {
    new.model <- changeMaxPostCutoff(model, maxPost.cutoff = opt$cutoff)
}
ptm <- chromstaR:::startTimedMessage("Saving changed model to file ...")
save(new.model, file = "chromstaR-result.RData")
chromstaR:::stopTimedMessage(ptm)

## Export peak calls
folder <- 'peak-calls'
if (!file.exists(folder)) { dir.create(folder) }
exportPeaks(model = new.model, filename = file.path(folder,""), trackname = paste0("mode-", model$mode, ", ", opt$changeWhat, " = ", opt$cutoff))
files <- list.files(file.path(folder,""), pattern='_peaks_')
for (file in files) {
    filename.new <- sub("_peaks_", "peaks_", sub(".gz", "", file))
    command <- paste0("mv ", file.path(folder, file), " ", file.path(folder, filename.new))
    print(command)
    system(command = command, wait = TRUE)
}

## Export chromatin states
folder <- 'chromatin-states'
if (!file.exists(folder)) { dir.create(folder) }
exportCombinations(model = new.model, filename = file.path(folder,""), trackname = paste0("mode-", model$mode, ", ", opt$changeWhat, " = ", opt$cutoff))
files <- list.files(file.path(folder,""), pattern='_combinations_')
for (file in files) {
    filename.new <- sub("_combinations_", "combinations_", sub(".gz", "", file))
    command <- paste0("mv ", file.path(folder, file), " ", file.path(folder, filename.new))
    print(command)
    system(command = command, wait = TRUE)
}

total.time <- proc.time() - ptm.start; message("Total elapsed time: ",round(total.time[3],2),"s")

