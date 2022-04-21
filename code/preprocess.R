preprocess <- function(dat, min.nzcts = 10) {
  size.factors <- colSums(dat)
  size.factors <- size.factors / mean(size.factors)
  gene_cts <- rowSums(dat > 0)
  dat <- dat[gene_cts >= min.nzcts, ]

  lunpc <- max(1 / min(size.factors) - 1 / max(size.factors), 1)
  fl.dat <- log1p(t(t(dat) / size.factors) / lunpc)

  return(list(
    fl.dat = fl.dat,
    size.factors = size.factors,
    pc = lunpc,
    excluded.genes = gene_cts < min.nzcts)
  )
}
