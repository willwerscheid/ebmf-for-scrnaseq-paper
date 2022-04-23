preprocess <- function(dat, min.nzcts = 10) {
  size.factors <- Matrix::colSums(dat)
  size.factors <- size.factors / mean(size.factors)
  lunpc <- max(1 / min(size.factors) - 1 / max(size.factors), 1)

  # Use liger for most variable gene selection.
  liger.dat <- rliger::createLiger(list(dat = dat))
  liger.dat <- rliger::normalize(liger.dat)
  liger.dat <- rliger::selectGenes(liger.dat)

  gene_cts <- Matrix::rowSums(dat > 0)
  dat <- dat[gene_cts >= min.nzcts, ]

  return(list(
    counts = dat,
    sf = size.factors,
    pc = lunpc,
    var.genes = intersect(liger.dat@var.genes, rownames(dat))
  ))
}
