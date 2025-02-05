library(tidyverse)

go <- c("geneontology_Biological_Process",
        "geneontology_Cellular_Component",
        "geneontology_Molecular_Function")

run_GSEA <- function(LL, incl_thresh = 0.01) {
  all_res <- tibble()
  for (k in 1:ncol(LL)) {
    dat <- data.frame(
      Gene = names(LL[, k]),
      Loading = LL[, k]
    ) |>
      filter(Loading > max(Loading) * incl_thresh)
    cat("Running GSEA on component", k,
        "with", nrow(dat), "genes...\n")
    wgres <- WebGestaltR::WebGestaltR(
      organism = "mmusculus",
      enrichMethod = "GSEA",
      enrichDatabase = go,
      interestGene = dat, interestGeneType = "genesymbol",
      minNum = 10, maxNum = 500,
      topThr = 10, sigMethod = "top",
      perNum = 100,
      isOutput = FALSE
    )
    wgres <- wgres |>
      mutate(k = k, .before = geneSet)
    all_res <- all_res |>
      bind_rows(wgres)
  }
  return(all_res)
}

ebnmf <- readRDS("output/montoro-ebmf.rds")
ebnmf_LL <- ebnmf$fit$L.pm
ebnmf_GSEA <- run_GSEA(ebnmf_LL)
saveRDS(ebnmf_GSEA, "output/montoro-ebnmf-gsea.rds")

nmf <- readRDS("output/montoro-nmf.rds")
nmf_LL <- nmf$fit$w
rownames(nmf_LL) <- rownames(ebnmf_LL)
nmf_GSEA <- run_GSEA(nmf_LL)
saveRDS(nmf_GSEA, "output/montoro-nmf-gsea.rds")

tm <- readRDS("output/montoro-topics.rds")
tm_LL <- tm$fit$L
rownames(tm_LL) <- rownames(ebnmf_LL)
tm_GSEA <- run_GSEA(tm_LL)
saveRDS(tm_GSEA, "output/montoro-tm-gsea.rds")
