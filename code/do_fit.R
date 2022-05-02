# datfile must be list including counts, pc, sf, and var.genes
# packages needed: fastTopics, RcppML

do_fit <- function(datfile, method, K, select.genes, outfile) {
  if (method == "fasttopics") {
    return(fit_fasttopics(datfile, K, select.genes, outfile))
  } else if (method == "nmf") {
    return(fit_nmf(datfile, K, select.genes, outfile))
  } else if (method == "snn-ebmf") {
    return(fit_snn_ebmf(datfile, K, select.genes, outfile))
  } else if (method == "nn-ebmf") {
    return(fit_nn_ebmf(datfile, K, select.genes, outfile))
  }
}

fit_fasttopics <- function(datfile, K, select.genes, outfile) {
  pp.dat <- readRDS(datfile)

  if (select.genes) {
    dat <- pp.dat$counts
  } else {
    dat <- pp.dat$counts[pp.dat$var.genes, ]
  }

  rm(pp.dat)
  dat <- as(dat, "dgCMatrix")

  t0 <- Sys.time()
  fit <- fastTopics::fit_poisson_nmf(dat, K, numiter = 100)
  t1 <- Sys.time()

  saveRDS(list(t = t1 - t0, fit = fit), outfile)
}

fit_nmf <- function(datfile, K, select.genes, outfile) {
  pp.dat <- readRDS(datfile)

  dat <- log1p(t(t(pp.dat$counts) / pp.dat$sf))
  if (select.genes) {
    dat <- dat[pp.dat$var.genes, ]
  }

  rm(pp.dat)

  t0 <- Sys.time()
  ntrials <- 30
  best_obj <- Inf
  for (i in 1:ntrials) {
    fit <- RcppML::nmf(dat, K, maxit = 100, seed = i)
    obj <- sum((dat - fit$w %*% fit$h)^2)
    if (obj < best_obj) {
      best_obj <- obj
      best_fit <- fit
    }
  }
  t1 <- Sys.time()

  saveRDS(list(t = t1 - t0, fit = best_fit), outfile)
}

fit_snn_ebmf <- function(datfile, K, select.genes, outfile) {
  fit_ebmf(datfile, K, select.genes, outfile, nonnegative = FALSE)
}

fit_nn_ebmf <- function(datfile, K, select.genes, outfile) {
  fit_ebmf(datfile, K, select.genes, outfile, nonnegative = TRUE)
}

fit_ebmf <- function(datfile, K, select.genes, outfile, nonnegative) {
  pp.dat <- readRDS(datfile)

  dat <- log1p(t(t(pp.dat$counts) / pp.dat$sf))

  if (select.genes) {
    dat <- dat[pp.dat$var.genes, ]
  }

  set.seed(666)
  n <- ncol(dat)
  min.sd <- sd(log1p(rpois(1e7, 1 / n) / median(pp.dat$sf)))

  rm(datfile)

  if (nonnegative) {
    ebnm.fn = ebnm::ebnm_point_exponential
    init.fn = function(f) init.fn.default(f, dim.signs = c(1, 1))
  } else {
    ebnm.fn = c(ebnm::ebnm_point_laplace, ebnm::ebnm_point_exponential)
    init.fn = function(f) init.fn.default(f, dim.signs = c(0, 1))
  }

  t0 <- Sys.time()

  intercept <- list(
    matrix(rowMeans(dat), ncol = 1),
    matrix(rep(1, ncol(dat)), ncol = 1)
  )
  init.fl <- flash.init(dat, S = min.sd, var.type = 1) %>%
    flash.init.factors(
      intercept,
      ebnm.fn = as.ebnm.fn(prior_family = "point_normal", mode = "estimate")
    ) %>%
    flash.fix.factors(kset = 1, mode = 2L) %>%
    flash.backfit(verbose = 0)

  fl <- init.fl %>%
    flash.add.greedy(
      Kmax = min(K - 1, 10),
      ebnm.fn = ebnm.fn,
      init.fn = init.fn
    )

  while (fl$n.factors < K) {
    fl <- fl %>%
      flash.backfit(maxiter = 10, verbose = 3) %>%
      flash.add.greedy(
        Kmax = min(K - fl$n.factors, 10),
        ebnm.fn = ebnm.fn,
        init.fn = init.fn,
        verbose = 1
      )
  }

  fl <- fl %>%
    flash.backfit(maxiter = 100, verbose = 3)

  t1 <- Sys.time()

  fl <- flash.reorder.factors(fl, c(1, order(fl$pve[-1], decreasing = TRUE) + 1))

  saveRDS(list(t = t1 - t0, fit = fl), outfile)
}
