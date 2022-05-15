# Pulse-seq dataset from Montoro et al.

library(Matrix)

zz <- download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE103354&format=file&file=GSE103354%5FPulseSeq%5FUMI%5Fcounts%2Erds%2Egz",
                    destfile = "../data/pulseseq.rds.gz")
R.utils::gunzip("../data/pulseseq.rds.gz")

dat <- readRDS("../data/pulseseq.rds")
file.remove("../data/pulseseq.rds")

source("../code/preprocess.R")
pp.dat <- preprocess(dat)
saveRDS(pp.dat, "../data/montoro.rds")

mean.expr <- rowSums(pp.dat$counts) / ncol(pp.dat$counts)
var.gene.mean.expr <- mean.expr[pp.dat$var.genes]
saveRDS(
  list(mean.expr = mean.expr, var.gene.mean.expr = var.gene.mean.expr),
  "../data/montoro-mean-expr.rds"
)
