# From Deng et al.

library(Matrix)
library(singleCellRNASeqMouseDeng2014)

dat <- Matrix(exprs(Deng2014MouseESC))
meta <- pData(Deng2014MouseESC)
colnames(dat) <- paste0(colnames(dat), "_", meta$cell_type, "_", meta$embryo_id)

source("../code/preprocess.R")
pp.dat <- preprocess(dat)
saveRDS(pp.dat, "../data/deng.rds")

mean.expr <- rowSums(pp.dat$counts) / ncol(pp.dat$counts)
var.gene.mean.expr <- mean.expr[pp.dat$var.genes]
saveRDS(
  list(mean.expr = mean.expr, var.gene.mean.expr = var.gene.mean.expr),
  "../data/deng-mean-expr.rds"
)
