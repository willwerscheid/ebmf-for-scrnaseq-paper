library(Matrix)
library(flashier)

do_ebmf_fits <- function(pp.dat, K, folder.name, maxiter = 500) {
  set.seed(666)
  n <- ncol(pp.dat$fl.dat)
  min.sd <- sd(log1p(rpois(1e7, 1/n) / pp.dat$pc / median(pp.dat$size.factors)))
  disp.min.sd <- function(new, old, k) {
    return(sqrt(min(1 / ff.tau(new))))
  }

  # Intercept factor.
  intercept <- list(
    matrix(rowMeans(pp.dat$fl.dat), ncol = 1),
    matrix(rep(1, ncol(pp.dat$fl.dat)), ncol = 1)
  )
  init.fl <- flash.init(pp.dat$fl.dat, S = min.sd, var.type = 1) %>%
    flash.init.factors(
      intercept,
      ebnm.fn = as.ebnm.fn(prior_family = "point_normal", mode = "estimate")
    ) %>%
    flash.fix.factors(kset = 1, mode = 2L) %>%
    flash.backfit(verbose = 0)

  # SNMF fit.
  t0 <- Sys.time()
  snmf.fl <- init.fl %>%
    flash.add.greedy(
      Kmax = min(K - 1, 10),
      ebnm.fn = c(ebnm::ebnm_point_laplace, ebnm::ebnm_point_exponential),
      init.fn = function(f) init.fn.default(f, dim.signs = c(0, 1))
    )

  while (snmf.fl$n.factors < K) {
    snmf.fl <- snmf.fl %>%
      flash.backfit(maxiter = 10, verbose = 3) %>%
      flash.add.greedy(
        Kmax = min(K - snmf.fl$n.factors, 10),
        ebnm.fn = c(ebnm::ebnm_point_laplace, ebnm::ebnm_point_exponential),
        init.fn = function(f) init.fn.default(f, dim.signs = c(0, 1)),
        verbose = 1
      )
  }

  snmf.fl <- snmfl.fl %>%
    flash.backfit(maxiter = maxiter, verbose = 3)

  t1 <- Sys.time()
  snmf.t <- t1 - t0

  saveRDS(list(fl = snmf.fl, t = snmf.t), paste0("./output/", folder.name, "/snmf.rds"))
  rm(snmf.fl)

  # NMF fit.
  t0 <- Sys.time()
  nmf.fl <- init.fl %>%
    flash.add.greedy(
      Kmax = min(K - 1, 10),
      ebnm.fn = ebnm::ebnm_point_exponential,
      init.fn = function(f) init.fn.default(f, dim.signs = c(1, 1))
    )

  while (nmf.fl$n.factors < K) {
    nmf.fl <- nmf.fl %>%
      flash.backfit(maxiter = 10, verbose = 3) %>%
      flash.add.greedy(
        Kmax = min(K - nmf.fl$n.factors, 10),
        ebnm.fn = ebnm::ebnm_point_exponential,
        init.fn = function(f) init.fn.default(f, dim.signs = c(1, 1)),
        verbose = 1
      )
  }

  nmf.fl <- nmfl.fl %>%
    flash.backfit(maxiter = maxiter, verbose = 3)

  t1 <- Sys.time()
  nmf.t <- t1 - t0

  saveRDS(list(fl = nmf.fl, t = nmf.t), paste0("./output/", folder.name, "/nmf.rds"))
  rm(nmf.fl)
}
