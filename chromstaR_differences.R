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
  make_option(c("-m", "--minWidth"), type = 'numeric', default = 300,
              help="Minimum width in base-pairs for differential regions. Default is %default."),
  make_option(c("-s", "--differentialScore"), type = 'numeric', default = 0.9999,
              help="Minimum differential score to detect differences. Default is %default."),
  make_option(c("-p", "--differentialPosterior"), type = 'numeric', default = 0.9999,
              help="Minimum differential posterior to detect pairwise differences. Default is %default.")
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
marks <- unique(model$info$mark)
conditions <- unique(model$info$condition)

## Get chromatin state differences
folder <- 'chromatin-diffs'
if (!file.exists(folder)) { dir.create(folder) }
s <- model$segments
for (icond1 in 1:(length(conditions)-1)) {
    for (icond2 in (icond1+1):length(conditions)) {
        dmask <- s$differential.score >= opt$differentialScore
        diffs <- s[dmask]
        diffs <- diffs[width(diffs) >= opt$minWidth]
        diffs$name <- paste0(conditions[icond1], ":", as.character(mcols(diffs)[,paste0("combination.", conditions[icond1])]), ", ", conditions[icond2], ":", as.character(mcols(diffs)[,paste0("combination.", conditions[icond2])]))
        trackname <- paste0("chromatin state differences between ", conditions[icond1], " and ", conditions[icond2], ", mode-", model$mode, ", differentialScore=", opt$differentialScore)
        filename <- paste0("differential_chromatin-states_", conditions[icond1], "-", conditions[icond2])
        exportGRangesAsBedFile(diffs, namecol = 'name', trackname = trackname, filename = file.path(folder, filename))
        command <- paste0("mv ", file.path(folder, paste0(filename, ".bed.gz")), " ", file.path(folder, paste0(filename, ".bed")))
        print(command)
        system(command = command, wait = TRUE)
    }
}

## Get chromatin state differences
folder <- 'pairwise-diffs'
if (!file.exists(folder)) { dir.create(folder) }
s <- model$segments
for (imark in 1:length(marks)) {
    for (icond1 in 1:(length(conditions)-1)) {
        id1 <- paste0(marks[imark], "-", conditions[icond1])
        maxpost1 <- s$maxPostInPeak[,grep(paste0('\\<', id1, '\\>'), colnames(s$maxPostInPeak))[1]]
        for (icond2 in (icond1+1):length(conditions)) {
            id2 <- paste0(marks[imark], "-", conditions[icond2])
            maxpost2 <- s$maxPostInPeak[,grep(paste0('\\<', id2, '\\>'), colnames(s$maxPostInPeak))[1]]
            
            ## Stuff in conditionX but not in conditionY ##
            for (i1 in 1:2) {
                if (i1==1) {
                    maxpostX <- maxpost1
                    maxpostY <- maxpost2
                    icondX <- icond1
                    icondY <- icond2
                } else if (i1==2) {
                    maxpostX <- maxpost2
                    maxpostY <- maxpost1
                    icondX <- icond2
                    icondY <- icond1
                }
                dmask <- (maxpostX - maxpostY) >= opt$differentialPosterior # select by maxPostInPeak
                diffs <- s[dmask]
                diffs <- reduce(diffs)
                diffs <- diffs[width(diffs) >= opt$minWidth]
                
                trackname <- paste0(marks[imark], ", peaks in ", conditions[icond1], " and not in ", conditions[icond2], ", mode-", model$mode, ", differentialPosterior=", opt$differentialPosterior)
                filename <- paste0("differential_peaks_", marks[imark], "_", conditions[icondX], "-not-", conditions[icondY])
                exportGRangesAsBedFile(diffs, trackname = trackname, filename = file.path(folder, filename))
                command <- paste0("mv ", file.path(folder, paste0(filename, ".bed.gz")), " ", file.path(folder, paste0(filename, ".bed")))
                print(command)
                system(command = command, wait = TRUE)
            }
        }
    }
}

total.time <- proc.time() - ptm.start; message("Total elapsed time: ",round(total.time[3],2),"s")

