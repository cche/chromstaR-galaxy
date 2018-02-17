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
  make_option(c("-i", "--inputfolder"), type = 'character',
              help="Folder with input files"),
  make_option(c("-F", "--format"), type = 'character',
              help="Format of the input files. Either 'bam' or 'bed'."),
  make_option(c("-x", "--experiment.table"),
              help="Experiment table"),
  make_option(c("-o", "--outputfolder"), type = 'character',
              help="Folder for the output files"),
  make_option(c("-f", "--configfile"), type = 'character', default = NULL,
              help="chromstaR configuration file. Default is %default."),
  make_option(c("-n", "--numCPU"), type = 'numeric', default = 1, 
              help="Number of CPUs to use. Default is %default."),
  make_option(c("-B", "--binsize"), type = 'numeric', default = 1000,
              help="Bin size. Default is %default."),
  make_option(c("-S", "--stepsize"), type = 'numeric', default = 500,
              help="Step size for bin offset. Default is %default."),
  make_option(c("-a", "--assembly"), type = 'character', default = NULL,
              help="Genome assembly. Default is %default."),
  make_option(c("-c", "--chromosomes"), type = 'character',
              help="Chromosomes to use in the analysis"),
  make_option(c("-D", "--remove.duplicate.reads"), type = 'logical', default = TRUE, 
              help="Whether or not to remove duplicate reads. Default is %default."),
  make_option(c("-Q", "--min.mapq"), type = 'numeric', default = 10,
              help="Minimum mapping quality. Default is %default."),
  make_option(c("-P", "--prefit.on.chr"), type = 'character', default = NULL,
              help="Prefit on chromosome. Default is %default."),
  make_option(c("-U", "--eps.univariate"), type = 'numeric', default = 0.1,
              help="Univariate epsilon for Baum-Welch. Default is %default."),
  make_option(c("-T", "--max.time"), type = 'numeric', default = NULL,
              help="Maximum time for Baum-Welch. Default is %default."),
  make_option(c("-I", "--max.iter"), type = 'numeric', default = 5000,
              help="Maximum iterations for Baum-Welch. Default is %default."),
  make_option(c("-R", "--read.cutoff.absolute"), type = 'numeric', default = 500,
              help="Read count cutoff. Default is %default."),
  make_option(c("-k", "--keep.posteriors"),type = 'logical', default = TRUE, 
              help="Whether or not to keep posteriors. Default is %default."),
  make_option(c("-m", "--mode"), type = 'character', default = 'differential',
              help="Mode of analysis. One of 'differential', 'combinatorial', 'full', 'separate'. Default is %default."),
  make_option(c("-X", "--max.states"), type = 'numeric', default = 128,
              help="Maximum number of states. Default is %default."),
  make_option(c("-p", "--per.chrom"),type = 'logical', default = TRUE, 
              help="Whether or not to do the multivariate part per chromosome. Default is %default."),
  make_option(c("-M", "--eps.multivariate"), type = 'numeric', default = 0.01,
              help="Multivariate epsilon for Baum-Welch. Default is %default."),
  make_option(c("-e", "--exclusive.table"), type = 'character', default = NULL,
              help="Exclusive table. Default is %default.")
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

Chromstar(inputfolder = opt$inputfolder, experiment.table = opt$experiment.table, outputfolder = opt$outputfolder, configfile = opt$configfile, numCPU = opt$numCPU, binsize = opt$binsize, stepsize = opt$stepsize, assembly = opt$assembly, chromosomes = opt$chromosomes, remove.duplicate.reads = opt$remove.duplicate.reads, min.mapq = opt$min.mapq, format = opt$format, prefit.on.chr = opt$prefit.on.chr, eps.univariate = opt$eps.univariate, max.time = opt$max.time, max.iter = opt$max.iter, read.cutoff.absolute = opt$read.cutoff.absolute, keep.posteriors = opt$keep.posteriors, mode = opt$mode, max.states = opt$max.states, per.chrom = opt$per.chrom, eps.multivariate = opt$eps.multivariate, exclusive.table = opt$exclusive.table)

## Print information about produced files
print("Files produced:")
print(list.files(opt$outputfolder, recursive=TRUE))
## Move main file to top level folder
command <- paste0("mv ", file.path(opt$outputfolder, 'combined', paste0('combined_mode-', opt$mode, '_binsize', opt$binsize, '_stepsize', opt$stepsize, '.RData')), " ", file.path('chromstaR-result.RData'))
print(command)
system(command = command, wait = TRUE)
## Move peak calls to separate folder
folder.new <- 'peak-calls'
if (!file.exists(folder.new)) { dir.create(folder.new) }
files <- list.files(file.path(opt$outputfolder, 'BROWSERFILES'), pattern='_peaks_')
for (file in files) {
    filename.new <- sub(".*_peaks_", "peaks_", sub(".gz", "", file))
    command <- paste0("mv ", file.path(opt$outputfolder, 'BROWSERFILES', file), " ", file.path(folder.new, filename.new))
    print(command)
    system(command = command, wait = TRUE)
}
## Move chromatin states to separate folder
folder.new <- 'chromatin-states'
if (!file.exists(folder.new)) { dir.create(folder.new) }
files <- list.files(file.path(opt$outputfolder, 'BROWSERFILES'), pattern='_combinations_')
for (file in files) {
    filename.new <- sub(".*_combinations_", "combinations_", sub(".gz", "", file))
    command <- paste0("mv ", file.path(opt$outputfolder, 'BROWSERFILES', file), " ", file.path(folder.new, filename.new))
    print(command)
    system(command = command, wait = TRUE)
}
## Move univariate plots to separate folder
folder.new <- 'univariate-fits'
if (!file.exists(folder.new)) { dir.create(folder.new) }
files <- list.files(file.path(opt$outputfolder, 'PLOTS', 'univariate-distributions'), pattern='\\.png')
for (file in files) {
    filename.new <- sub("_binsize.*png", ".png", file)
    command <- paste0("mv ", file.path(opt$outputfolder, 'PLOTS', 'univariate-distributions', file), " ", file.path(folder.new, filename.new))
    print(command)
    system(command = command, wait = TRUE)
}

total.time <- proc.time() - ptm.start; message("Total elapsed time: ",round(total.time[3],2),"s")

