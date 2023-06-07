library(flashier)
library(tidyverse)

sim_data <- function(n = 100,
                     p = 200,
                     K = 6,
                     L.nn = 10,
                     F.nn = 20,
                     se = 0.1,
                     K.dense = 1,
                     seed = 666) {
  set.seed(seed)

  LL <- matrix(rexp(n * K), nrow = n, ncol = K)
  FF <- matrix(rexp(p * K), nrow = p, ncol = K)

  # "Mean" factor.
  LL[, 1] <- 3

  # Additional sparse nonnegative factors.
  for (k in (K.dense + 1):K) {
    L.nn.idx <- seq((k - K.dense - 1) * L.nn + 1, (k - K.dense) * L.nn)
    F.nn.idx <- seq((k - K.dense - 1) * (F.nn / 2) + 1, (k - K.dense + 1) * (F.nn / 2))
    LL[setdiff(1:n, L.nn.idx), k] <- 0
    FF[setdiff(1:p, F.nn.idx), k] <- 0
  }

  # Add normal noise.
  Y <- LL %*% t(FF) + rnorm(n * p, sd = se)

  # Add a constant (which can be absorbed by mean factor) to ensure nonnegativity.
  Y <- Y - min(Y)

  return(list(LL = LL, FF = FF, Y = Y))
}

set.seed(666)
simdat <- sim_data(K = 6)

# Use true value of K.
t0 <- Sys.time()
nnlm_res1 <- RcppML::nmf(simdat$Y, k = 6, verbose = FALSE, tol = 1e-6, maxit = 10000L)
nnlm_t1 <- Sys.time() - t0

niter <- 30
t0 <- Sys.time()
nnlm_res2 <- nnlm_res1
best_obj <- mean((simdat$Y - nnlm_res1$w %*% diag(nnlm_res1$d) %*% nnlm_res1$h)^2)
for (i in 2:niter) {
  nnlm_res <- RcppML::nmf(simdat$Y, k = 6, verbose = FALSE, tol = 1e-6, maxit = 10000L)
  obj <- mean((simdat$Y - nnlm_res$w %*% diag(nnlm_res$d) %*% nnlm_res$h)^2)
  if (obj < best_obj) {
    nnlm_res2 <- nnlm_res
    best_obj <- obj
  }
}
nnlm_t2 <- (Sys.time() - t0) + nnlm_t1

t0 <- Sys.time()
ebnmf_res <- flash.init(simdat$Y) %>%
  flash.set.verbose(0) %>%
  flash.add.greedy(
    Kmax = 6,
    ebnm.fn = ebnm::ebnm_point_exponential,
    init.fn = function(fl) init.fn.default(fl, dim.signs = c(1, 1))
  ) %>%
  flash.backfit() %>%
  flash.nullcheck()
ebnmf_t <- Sys.time() - t0

nnlm_idx1 <- c(1, apply(cor(simdat$LL[, -1], nnlm_res1$W[, -1]), 1, which.max) + 1)
nnlm_idx2 <- c(1, apply(cor(simdat$LL[, -1], nnlm_res2$W[, -1]), 1, which.max) + 1)
ebnmf_idx <- c(1, apply(cor(simdat$LL[, -1], ebnmf_res$L.pm[, -1]), 1, which.max) + 1)

to_tibble <- function(mat, type, dim, col_idx) {
  mat <- scale(mat, center = FALSE, scale = apply(mat, 2, function(x) max(abs(x))))
  mat <- mat[, col_idx]
  return(
    as_tibble(mat, .name_repair = "unique") %>%
      mutate(row = row_number()) %>%
      pivot_longer(!row, names_to = "k", values_to = "value") %>%
      add_column(type = type, dim = dim)
  )
}

# Suppress annoying "new names" messages.
suppressMessages({
  tib <- to_tibble(simdat$LL, "True", "L", 1:6) %>%
    bind_rows(to_tibble(simdat$FF, "True", "F", 1:6)) %>%
    bind_rows(to_tibble(ebnmf_res$L.pm, "EBNMF", "L", ebnmf_idx)) %>%
    bind_rows(to_tibble(ebnmf_res$F.pm, "EBNMF", "F", ebnmf_idx))
  tib <- tib %>%
    bind_rows(to_tibble(nnlm_res1$w, "NMF, one run", "L", nnlm_idx1)) %>%
    bind_rows(to_tibble(t(nnlm_res1$h), "NMF, one run", "F", nnlm_idx1))
  tib <- tib %>%
    bind_rows(to_tibble(nnlm_res2$w, "NMF, best run", "L", nnlm_idx2)) %>%
    bind_rows(to_tibble(t(nnlm_res2$h), "NMF, best run", "F", nnlm_idx2))
  tib <- tib %>%
    mutate(k = as.numeric(str_remove_all(k, "\\."))) %>%
    mutate(
      type = factor(type, levels = c(
        "True", "EBNMF", "NMF, one run", "NMF, best run"
      )),
      dim = factor(dim, levels = c("L", "F"))
    )
})

ggplot(tib, aes(x = k, y = row, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(high = "black") +
  facet_grid(rows = vars(dim), cols = vars(type), scales = "free", switch = "y") +
  theme_void() +
  theme(
    strip.text.x = element_text(size = 8, margin = margin(2, 2, 2, 2), family = "serif"),
    strip.text.y.left = element_text(size = 12, margin = margin(4, 4, 4, 4), family = "serif", angle = 0),
    legend.position = "none"
  )

ggsave("../figs/sim_res.png", width = 100, height = 60, units = "mm")

tib_poster <- tib %>%
  filter(type != "NMF, one run") %>%
  mutate(type = fct_drop(fct_recode(type, `Vanilla NMF` = "NMF, best run")))
ggplot(tib_poster, aes(x = k, y = row, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(high = "black") +
  facet_grid(rows = vars(dim), cols = vars(type), scales = "free", switch = "y") +
  theme_void() +
  theme(
    strip.text.x = element_text(size = 14, margin = margin(2, 2, 2), family = "serif"),
    strip.text.y.left = element_text(size = 14, margin = margin(4, 4, 4), family = "serif", angle = 0),
    legend.position = "none"
  )

ggsave("./figs/sim_res_poster.png", width = 100, height = 60, units = "mm")

cat("Time to run EBNMF:", format(ebnmf_t, units= "auto"), "\n")
cat("Time to do", niter, "NMF runs:", format(nnlm_t2, units= "auto"), "\n")
