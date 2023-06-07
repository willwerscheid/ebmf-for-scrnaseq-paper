library(Matrix)
library(flashier)

source("./code/do_ebmf_fit.R")

pulseseq <- readRDS("./data/pulseseq_pp.rds")

do_ebmf_fit(pulseseq, 25, FALSE, "montoro/snmf.rds", maxiter = 250)
do_ebmf_fit(pulseseq, 25, TRUE, "montoro/nmf.rds", maxiter = 250)
