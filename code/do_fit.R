# datfile must be list including counts, pc, sf, and var.genes
# packages needed: RcppML, fastTopics, pCMF

do_fit <- function(datfile, method, K, select.genes, outfile) {
  if (method == "nmf-log") {
    return(fit_nmf(datfile, K, select.genes, outfile, link = "log"))
  } else if (method == "nmf-identity") {
    return(fit_nmf(datfile, K, select.genes, outfile, link = "identity"))
  } else if (method == "fasttopics") {
    return(fit_fasttopics(datfile, K, select.genes, outfile))
  } else if (method == "pcmf") {
    return(fit_pcmf(datfile, K, select.genes, outfile))
  } else if (method == "ebmf-log") {
    return(fit_ebmf(datfile, K, select.genes, outfile, link = "log"))
  } else if (method == "ebmf-identity") {
    return(fit_ebmf(datfile, K, select.genes, outfile, link = "identity"))
  }
}

get_data <- function(datfile, select.genes, link) {
  pp.dat <- readRDS(datfile)

  dat <- t(t(pp.dat$counts) / pp.dat$sf)

  if (select.genes) {
    dat <- dat[pp.dat$var.genes, ]
  }

  if (link == "log") {
    dat <- log1p(dat)
  } else if (link == "identity") {
    dat <- dat / apply(dat, 1, sd)
  }

  return(dat)
}

fit_nmf <- function(datfile, K, select.genes, outfile, link) {
  dat <- get_data(datfile, select.genes, link)

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

fit_pcmf <- function(datfile, K, select.genes, outfile) {
  dat <- get_data(datfile, select.genes, link = "identity")

  t0 <- Sys.time()
  fit <- pCMF::pCMF(dat, K, zero_inflation = FALSE, iter_max = 100)
  t1 <- Sys.time()

  saveRDS(list(t = t1 - t0, fit = fit), outfile)

}

fit_ebmf <- function(datfile, K, select.genes, outfile, link) {
  dat <- get_data(datfile, select.genes, link)

  pp.dat <- readRDS(datfile)
  set.seed(666)
  min.sd <- sd(log1p(rpois(1e7, 1 / ncol(dat)) / median(pp.dat$sf)))
  rm(datfile)

  t0 <- Sys.time()

  ebnm.fn = ebnm::ebnm_point_exponential
  init.fn = function(f) init.fn.default(f, dim.signs = c(1, 1))

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
