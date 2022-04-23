#! /usr/bin/env Rscript
#

# Load a few packages.
library(Matrix)
library(optparse)
library(flashier)

source("../code/do_fit.R")

# Process the command-line arguments.
parser <- OptionParser()
parser <- add_option(parser,"--datfile",type="character",default="dat.rds")
parser <- add_option(parser,"--method",type = "character",default = "nmf")
parser <- add_option(parser,c("--k","-k"),type = "integer",default = 3)
parser <- add_option(parser,c("--selgenes"),type = "logical",default = TRUE)
parser <- add_option(parser,c("--outfile","-o"),type="character",default="out.rds")
opts   <- parse_args(parser)

rm(parser)

res <- do_fit(opts$datfile, opts$method, opts$k, opts$selgenes, opts$outfile)

print(sessionInfo())
