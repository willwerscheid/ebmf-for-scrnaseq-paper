library(Matrix)
library(flashier)

source("../code/do_fit.R")

do_ps_fit <- function(method, K, outfile) {
  do_fit(datfile = "../data/pijuan-sala.rds", method = method, K = K,
         select.genes = TRUE, outfile = outfile)
}

do_montoro_fit <- function(method, K, outfile) {
  do_fit(datfile = "../data/montoro.rds", method = method, K = K,
         select.genes = TRUE, outfile = outfile)
}

do_ps_fit(method = "nmf-log",    K = 40, outfile = "../output/ps-nmf.rds")
do_ps_fit(method = "ebmf-log",   K = 20, outfile = "../output/ps-ebmf-20.rds")
do_ps_fit(method = "ebmf-log",   K = 30, outfile = "../output/ps-ebmf-30.rds")
do_ps_fit(method = "ebmf-log",   K = 40, outfile = "../output/ps-ebmf-40.rds")
do_ps_fit(method = "fasttopics", K = 20, outfile = "../output/ps-topics-20.rds")
do_ps_fit(method = "fasttopics", K = 30, outfile = "../output/ps-topics-30.rds")
do_ps_fit(method = "fasttopics", K = 40, outfile = "../output/ps-topics-40.rds")

do_montoro_fit(method = "nmf-log",    K = 30, outfile = "../output/montoro-nmf.rds")
do_montoro_fit(method = "ebmf-log",   K = 30, outfile = "../output/montoro-ebmf.rds")
do_montoro_fit(method = "fasttopics", K = 30, outfile = "../output/montoro-topics.rds")
