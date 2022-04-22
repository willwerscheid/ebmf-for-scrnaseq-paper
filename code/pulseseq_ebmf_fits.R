library(Matrix)
library(flashier)

source("./code/do_ebmf_fits.R")

pulseseq <- readRDS("./data/pulseseq_pp.rds")

do_ebmf_fit(pulseseq, 2, FALSE, "montoro/snmf.rds", maxiter = 5)
do_ebmf_fit(pulseseq, 2, TRUE, "montoro/nmf.rds", maxiter = 5)
