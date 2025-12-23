# datfile must be list including counts, pc, sf, and var.genes
# packages needed: RcppML, fastTopics

maxiter <- 200

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
  } else if (method == "glmpca") {
    return(fit_glmpca(datfile, K, select.genes, outfile))
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
    fit <- RcppML::nmf(dat, K, maxit = maxiter, seed = i)
    obj <- sum((dat - fit$w %*% diag(fit$d) %*% fit$h)^2)
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
    dat <- pp.dat$counts[pp.dat$var.genes, ]
  } else {
    dat <- pp.dat$counts
  }

  rm(pp.dat)
  dat <- as(dat, "CsparseMatrix")

  t0 <- Sys.time()
  fit <- fastTopics::fit_poisson_nmf(
    dat,
    K,
    control = list(extrapolate = TRUE),
    numiter = maxiter
  )
  t1 <- Sys.time()

  saveRDS(list(t = t1 - t0, fit = fit), outfile)
}

fit_glmpca <- function(datfile, K, select.genes, outfile) {
  pp.dat <- readRDS(datfile)

  if (select.genes) {
    dat <- pp.dat$counts[pp.dat$var.genes, ]
  } else {
    dat <- pp.dat$counts
  }

  rm(pp.dat)
  dat <- as(dat, "CsparseMatrix")

  t0 <- Sys.time()
  fit0 <- fastglmpca::init_glmpca_pois(dat, K = K)
  fit <- fastglmpca::fit_glmpca_pois(dat, fit0 = fit0)
  t1 <- Sys.time()

  saveRDS(list(t = t1 - t0, fit = fit), outfile)
}

fit_pcmf <- function(datfile, K, select.genes, outfile) {
  pp.dat <- readRDS(datfile)

  if (select.genes) {
    dat <- pp.dat$counts
  } else {
    dat <- pp.dat$counts[pp.dat$var.genes, ]
  }

  rm(pp.dat)
  dat <- as.matrix(dat)

  t0 <- Sys.time()
  fit <- pCMF::pCMF(dat, K, zero_inflation = FALSE, iter_max = maxiter)
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

  fl <- flash_init(dat, S = min.sd, var_type = 1) |>
    flash_greedy(
      Kmax = min(K, 10),
      ebnm_fn = ebnm::ebnm_point_exponential
    )

  while (fl$n_factors < K) {
    fl <- fl |>
      flash_backfit(maxiter = 10, verbose = 3) |>
      flash_greedy(
        Kmax = min(K - fl$n_factors, 10),
        ebnm_fn = ebnm::ebnm_point_exponential,
        verbose = 1
      )
  }

  fl <- fl |>
    flash_backfit(maxiter = maxiter, verbose = 3)

  t1 <- Sys.time()

  fl <- flash_factors_reorder(fl, c(1, order(fl$pve[-1], decreasing = TRUE) + 1))
  fl$sampler <- NULL
  fl$flash.fit <- NULL

  saveRDS(list(t = t1 - t0, fit = fl), outfile)
}
