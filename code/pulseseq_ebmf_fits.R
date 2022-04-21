source("./code/do_ebmf_fits.R")

pulseseq <- readRDS("./data/pulseseq_pp.rds")

do_ebmf_fits(pulseseq, 25, "montoro", maxiter = 200)
