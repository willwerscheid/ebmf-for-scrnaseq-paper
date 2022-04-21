# Pulse-seq dataset from Montoro et al. Run from project home directory.

library(Matrix)

zz <- download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE103354&format=file&file=GSE103354%5FPulseSeq%5FUMI%5Fcounts%2Erds%2Egz",
                    destfile = "./data/pulseseq.rds.gz")
R.utils::gunzip("./data/pulseseq.rds.gz")

pulseseq <- readRDS("./data/pulseseq.rds")

source("./code/preprocess.R")
pulseseq.pp <- preprocess(pulseseq)
saveRDS(pulseseq.pp, "./data/pulseseq_pp.rds")

zz <- file.remove("./data/pulseseq.rds")
