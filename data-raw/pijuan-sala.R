# Chimera-tal1 dataset from Pijuan-Sala et al.

library(tidyverse)
source("../code/preprocess.R")

dat <- Matrix::readMM("../../EmbryoTimecourse2018/download/chimera-tal1/raw_counts.mtx")

# Genes:
rownames(dat) <- read_tsv(
  "../../EmbryoTimecourse2018/download/chimera-tal1/genes.tsv",
  col_names = FALSE
)$X2

# Cells:
cells <- read_tsv("../../EmbryoTimecourse2018/download/chimera-tal1/meta.tab") %>%
  unite(CellName, 2:7) %>%
  pull(CellName)
colnames(dat) <- cells

pp.dat <- preprocess(dat)
saveRDS(pp.dat, "../data/pijuan-sala.rds")

mean.expr <- Matrix::rowSums(pp.dat$counts) / ncol(pp.dat$counts)
var.gene.mean.expr <- mean.expr[pp.dat$var.genes]
saveRDS(
  list(mean.expr = mean.expr, var.gene.mean.expr = var.gene.mean.expr),
  "../data/pijuan-sala-mean-expr.rds"
)
