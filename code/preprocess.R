preprocess <- function(dat, genelist = NULL, min.nzcts = 10) {
  size.factors <- Matrix::colSums(dat)
  size.factors <- size.factors / mean(size.factors)
  lunpc <- max(1 / min(size.factors) - 1 / max(size.factors), 1)
  gene_cts <- Matrix::rowSums(dat > 0)

  if (!is.null(genelist)) {
    rownames(dat) <- genelist
  }
  dat <- dat[gene_cts >= min.nzcts, ]

  liger.dat <- rliger::createLiger(list(dat = dat))
  liger.dat <- rliger::normalize(liger.dat)
  liger.dat <- rliger::selectGenes(liger.dat)

  return(list(
    counts = dat,
    sf = size.factors,
    pc = lunpc,
    var.genes = liger.dat@var.genes
  ))
}
