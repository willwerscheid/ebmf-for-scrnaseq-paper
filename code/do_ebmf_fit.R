do_ebmf_fit <- function(pp.dat, K, nonnegative, file.name, maxiter = 500) {
  set.seed(666)
  n <- ncol(pp.dat$fl.dat)
  min.sd <- sd(log1p(rpois(1e7, 1/n) / pp.dat$pc / median(pp.dat$size.factors)))
  disp.min.sd <- function(new, old, k) {
    return(sqrt(min(1 / ff.tau(new))))
  }

  if (nonnegative) {
    ebnm.fn = ebnm::ebnm_point_exponential
    init.fn = function(f) init.fn.default(f, dim.signs = c(1, 1))
  } else {
    ebnm.fn = c(ebnm::ebnm_point_laplace, ebnm::ebnm_point_exponential)
    init.fn = function(f) init.fn.default(f, dim.signs = c(0, 1))
  }

  t0 <- Sys.time()

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
    flash.backfit(maxiter = maxiter, verbose = 3)

  t1 <- Sys.time()

  saveRDS(list(fl = fl, t = t1 - t0), paste0("./output/", file.name))
  return(TRUE)
}
