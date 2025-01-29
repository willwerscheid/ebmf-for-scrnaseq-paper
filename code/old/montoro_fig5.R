library(clusterProfiler)
library(org.Mm.eg.db)
library(ggrepel)
library(gt)
library(tidyverse)
library(cowplot)

# exprmean <- log10(readRDS("data/montoro-mean-expr.rds")$var.gene.mean.expr)

ebnmf <- readRDS("output/montoro-ebmf.rds")
ebnmf_FF <- ebnmf$fit$F.pm

calc_factor_sizes <- function(FF, LL, D = 1) {
  l1_norm <- apply(FF, 2, max)
  LL <- t(t(LL) * l1_norm * D)
  return(apply(LL, 2, function(x) sqrt(sum(x^2))))
}

ebnmf_widths <- calc_factor_sizes(ebnmf_FF, ebnmf$fit$L.pm)

external_info <- tibble(
  CellType = sapply(strsplit(rownames(ebnmf_FF), "_"), `[[`, 5),
  TimePoint = sapply(strsplit(rownames(ebnmf_FF), "_"), `[[`, 2),
  Replicate = sapply(strsplit(rownames(ebnmf_FF), "_"), `[[`, 3),
  Fluor = sapply(strsplit(rownames(ebnmf_FF), "_"), `[[`, 4)
) |>
  mutate(
    CellTypeCoarse = case_when(
      CellType %in% c("Basal", "Proliferating") ~ "Basal/Prolif",
      str_detect(CellType, "Club|Goblet") ~ "Club/Goblet",
      CellType == "Ciliated" ~ "Ciliated",
      TRUE ~ "Tuft/NEC/Iono"
    )
  ) |>
  mutate(
    CellType = factor(CellType, levels = c(
      "Basal",
      "Proliferating",
      "Club",
      "Club (hillock-associated)",
      "Goblet.progenitor",
      "Goblet.1",
      "Goblet.2",
      "Tuft.progenitor",
      "Tuft.1",
      "Tuft.2",
      "Neuroendocrine",
      "Ionocyte",
      "Ciliated"
    )),
    CellTypeCoarse = factor(CellTypeCoarse, levels = c(
      "Basal/Prolif",
      "Club/Goblet",
      "Tuft/NEC/Iono",
      "Ciliated"
    )),
    TimePoint = factor(TimePoint),
    Replicate = factor(Replicate),
    Fluor = factor(Fluor)
  ) |>
  mutate(CTxTP = fct_cross(TimePoint, CellTypeCoarse),
         CTxFlo = fct_cross(Fluor, CellTypeCoarse),
         TPxRep = fct_cross(Replicate, TimePoint))

make_tib <- function(FF, widths, info_col, ksort, rescale = TRUE) {
  info <- external_info[[info_col]]

  # Downsample the number of cells.
  set.seed(666)
  cell_idx <- numeric(0)
  for (k in levels(info)) {
    which_idx <- which(info == k)
    # Downsample common cell types. Duplicate rare ones.
    which_idx <- sample(which_idx, 200, replace = (length(which_idx) < 200))
    cell_idx <- c(cell_idx, which_idx)
  }

  FF <- cbind(rowSums(FF[cell_idx, -ksort]), FF[cell_idx, ksort])
  info <- info[cell_idx]

  tsne_res <- Rtsne::Rtsne(
    FF,
    dims = 1,
    pca = FALSE,
    normalize = FALSE,
    perplexity = 100,
    theta = 0.1,
    max_iter = 1000,
    eta = 200,
    check_duplicates = FALSE
  )$Y[, 1]

  sort_idx <- order(info, tsne_res)
  FF <- FF[sort_idx, ]
  colnames(FF) <- c("Other", ksort)

  if (rescale) {
    FF <- scale(FF, center = FALSE, scale = apply(FF, 2, max))
  }

  tib <- as_tibble(FF) |>
    mutate(
      CellIdx = row_number(),
      Info = info
    ) |>
    pivot_longer(
      -c(CellIdx, Info),
      names_to = "Component",
      values_to = "Loading",
      values_drop_na = TRUE
    ) |>
    mutate(Component = factor(Component, levels = colnames(FF)))

  width_tib <- tibble(
    Component = c(as.character(ksort), "Other"),
    Width = c(widths[ksort], ifelse(rescale, sum(widths[-ksort]), 1))
  ) |>
    mutate(x_end = cumsum(Width),
           x_pos = x_end - Width / 2)

  tib <- tib |>
    left_join(width_tib, by = "Component")

  return(tib)
}

get_type <- function(k) {
  lmod <- lm(tm_FF[, k] ~ CellType + TimePoint + Replicate + Fluor,
             data = external_info)
  pvals <- coef(summary(lmod))[-1, 4]
  which(pvals == min(pvals))
}
all_types <- lapply(1:30, get_type)
ct_ksort <- which(sapply(all_types, function(k) any(1:12 %in% k))) + 1
tp_ksort <- which(sapply(all_types, function(k) any(13:14 %in% k))) + 1
rep_ksort <- which(sapply(all_types, function(k) any(15:16 %in% k))) + 1
flo_ksort <- which(sapply(all_types, function(k) any(17 %in% k))) + 1

ct_ksort <- c(4, 5, 9, 21, 23, 17, 3, 7, 26, 15, 25, 28, 29,
              18, 19, 14, 20, 6, 13, 22)
ct_tib <- make_tib(ebnmf_FF, ebnmf_widths, "CellType", ct_ksort)

flo_ksort <- c(9, 15, 3)
flo_tib <- make_tib(ebnmf_FF, ebnmf_widths, "Fluor", flo_ksort)

tp_ksort <- c(2, 11, 12, 24) # 12/24 also rep (30)
tp_tib <- make_tib(ebnmf_FF, ebnmf_widths, "TPxRep", tp_ksort)

dense_ksort <- c(1, 8, 10, 16) # 2/16 also tp
dense_tib <- make_tib(ebnmf_FF, ebnmf_widths, "CellType", dense_ksort)

other_ksort <- c(27, 30)
other_tib <- make_tib(ebnmf_FF, ebnmf_widths, "CellType", other_ksort)

all_ksort <- unique(c(dense_ksort, tp_ksort, ct_ksort, other_ksort))
all_tib <- make_tib(ebnmf_FF, ebnmf_widths, "CellType", all_ksort)

nnlm <- readRDS("output/montoro-nmf.rds")
nnlm_FF <- t(nnlm$fit$h)
nnlm_widths <- calc_factor_sizes(nnlm_FF, nnlm$fit$w, nnlm$fit$d)

nnlm_tib <- make_tib(nnlm_FF, nnlm_widths, "CellType", 1:30)

nnlm_ct_ksort <- c(1, 7, 28, 21, 26, # dense
                   3, 10, 14, 13,
                   2, 11, 12, 6, 5, 4, 24, 19, # basal/prolif
                   9, 17, # club
                   18, 16, 25, # hillock
                   23, # goblet
                   29, # tuft
                   27, # NEC
                   30, 8, # ionocyte
                   20, 15, 22) # ciliated
nnlm_all_ksort <- c(1, 26, 28, 21, # dense
                    7, 3, 6, 5, 14, 10, 13, 19,
                    25, 16,
                    8, 30, # ionocyte
                    20,
                    2, 12, 11, 4, 24, # basal/prolif
                    9, 17, # club
                    23, # goblet
                    18, # hillock
                    29, # tuft
                    27, # NEC
                    22, 15) # ciliated
nnlm_ct_ksort <- c(1, # dense
                   6,
                   25, 16,
                   30, # ionocyte
                   2, 12, 4, 24, # basal/prolif
                   9, 17, # club
                   23, # goblet
                   18, # hillock
                   29, # tuft
                   27, # NEC
                   22, 15) # ciliated
nnlm_tib <- make_tib(nnlm_FF, nnlm_widths, "CellType", nnlm_ct_ksort)

tm <- readRDS("output/montoro-topics.rds")
tm_fit <- tm$fit
tm_fit$F <- tm$fit$L
tm_fit$L <- tm$fit$F
tm_fit$Fn <- tm$fit$Ln
tm_fit$Ln <- tm$fit$Fn
tm_fit$Fy <- tm$fit$Ly
tm_fit$Ly <- tm$fit$Fy
tm_fit <- fastTopics::poisson2multinom(tm_fit)
tm_FF <- tm_fit$L

tm_all_ksort <- c(
  30, 9, 29, 15, 21, 13, 26, 11, 27,
  17, 12, 24, 8, 10, 6, 4, #basal/prolif
  23, 3, 19, 18, 5, 28, 16, 7, 20, 22, #club/gob
  2, 1, 25
)
tm_ct_ksort <- c(
  30, 9, 29, 15, 13,
  12, 24, 8, 10, 6, #basal/prolif
  23, 3, 20, 16, 22, 7, #club/gob
  2, 1, 25
)
tm_tib <- make_tib(tm_FF, rep(1, 30), "CellType", tm_ct_ksort, rescale = FALSE)

next_tib <- tm_tib
struct_df <- next_tib |>
  mutate(Loading = Loading * Width) |>
  select(CellIdx:Loading) |>
  pivot_wider(names_from = Component, values_from = Loading)
Lmat <- struct_df |>
  select(-(1:2)) |>
  as.matrix()
group <- struct_df$Info
fastTopics::structure_plot(
  Lmat, topics = colnames(Lmat),
  grouping = group,
  loadings_order = 1:nrow(Lmat),
  gap = 40) +
  labs(y = "")

tib <- next_tib |>
  group_by(Info, CellIdx) |>
  summarize()
info_breaks <- c(1, which(tib$Info[-1] != tib$Info[-nrow(tib)]))
label_pos <- info_breaks / 2 + c(info_breaks[-1], nrow(tib)) / 2
max_idx <- max(tib$CellIdx)

heatmap_tib <- next_tib |> mutate(Method = "EBNMF")
p1 <- ggplot(heatmap_tib,
              aes(x = x_pos, y = -CellIdx, fill = Loading)) +
  geom_tile(aes(width = Width)) +
  scale_fill_gradient2(high = "black") +
  labs(x = "", y = "") +
  scale_y_continuous(breaks = -label_pos,
                     minor_breaks = NULL,
                     labels = levels(heatmap_tib$Info)) +
  theme_minimal() +
  geom_hline(yintercept = -info_breaks, linewidth = 0.1) +
  facet_wrap(~Method, scales = "free", ncol = 1) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12, margin = margin(2, 2, 2, 2), family = "serif"),
    axis.text.y = element_text(size = 8, family = "serif"),
    axis.text.x = element_blank(),
    axis.title = element_text(size = 8, margin = margin(2, 2, 2, 2), family = "serif"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  geom_text(
    data = filter(heatmap_tib, CellIdx == 1),
    aes(x = x_pos, y = -1.05 * max_idx, label = Component),
    size = 2
  )

p1





# plt <- ggplot(heatmap.tib, aes(x = SortedComponent, y = -CellIdx, fill = Loading)) +
#   geom_tile() +
#   scale_fill_gradient2(high = "black") +
#   labs(y = "") +
#   scale_y_continuous(breaks = -label_pos,
#                      minor_breaks = NULL,
#                      labels = levels(cell_type)) +
#   theme_minimal() +
#   geom_hline(yintercept = -cell_type_breaks, linewidth = 0.1) +
#   facet_wrap(~Method, scales = "free", ncol = 1) +
#   theme(
#     legend.position = "none",
#     strip.text = element_text(size = 12, margin = margin(2, 2, 2, 2), family = "serif"),
#     axis.text.y = element_text(size = 8, family = "serif"),
#     axis.text.x = element_text(size = 6, family = "serif"),
#     axis.title = element_text(size = 8, margin = margin(2, 2, 2, 2), family = "serif")
#   )
#
# ggsave("figs/montoro_heatmap.png", width = 163, height = 190, units = "mm")

# poster.tib.ebnmf <- heatmap.tib |>
#   filter(Method == "EBNMF") |>
#   mutate(Loading = ifelse(Component == "20", -Loading, Loading)) |>
#   mutate(CellType = fct_recode(CellType, `Club (hillock)` = "Club (hillock-associated)"))
# plt <- ggplot(poster.tib.ebnmf, aes(x = Component, y = -CellIdx, fill = Loading)) +
#   geom_tile() +
#   scale_fill_gradient2(low = "red", mid = "white", high = "black") +
#   labs(y = "", title = "EBNMF") +
#   scale_y_continuous(breaks = -label_pos,
#                      minor_breaks = NULL,
#                      labels = levels(poster.tib.ebnmf$CellType)) +
#   theme_minimal() +
#   geom_hline(yintercept = -cell_type_breaks, linewidth = 0.1) +
#   theme(
#     legend.position = "none",
#     plot.title = element_text(size = 16, family = "serif", hjust = 0.5),
#     axis.text.y = element_text(size = 13, family = "serif"),
#     axis.text.x = element_text(size = 9, family = "serif"),
#     axis.title = element_text(size = 14, margin = margin(2, 2, 2, 2), family = "serif")
#   )
#
# ggsave("figs/montoro_ebnmf_heatmap_poster.png", width = 163, height = 95, units = "mm")
#
# poster.tib.nmf <- heatmap.tib |>
#   filter(Method == "NMF") |>
#   mutate(Loading = ifelse(Component == "30", -Loading, Loading)) |>
#   mutate(CellType = fct_recode(CellType, `Club (hillock)` = "Club (hillock-associated)"))
# plt <- ggplot(poster.tib.nmf, aes(x = Component, y = -CellIdx, fill = Loading)) +
#   geom_tile() +
#   scale_fill_gradient2(low = "red", mid = "white", high = "black") +
#   labs(x = "", y = "", title = "NMF") +
#   scale_y_continuous(breaks = -label_pos,
#                      minor_breaks = NULL,
#                      labels = NULL) +
#   theme_minimal() +
#   geom_hline(yintercept = -cell_type_breaks, linewidth = 0.1) +
#   theme(
#     legend.position = "none",
#     plot.title = element_text(size = 16, family = "serif", hjust = 0.5),
#     axis.text.x = element_blank()
#   )
#
# ggsave("figs/montoro_nmf_heatmap_poster.png", width = 163 / 2, height = 95 / 2, units = "mm")

top.genes <- apply(ebnmf$fit$L.pm, 2, function(k) {
  rownames(ebnmf$fit$L.pm)[order(k, decreasing = TRUE)[1:12]]
})
colnames(top.genes) <- paste0("k=", formatC(1:ncol(top.genes), width = 2, flag = "0"))
genes.tib <- as_tibble(top.genes) |>
  pivot_longer(cols = everything(), names_to = "Component", values_to = "SYMBOL") |>
  mutate(Component = as.numeric(str_remove_all(Component, "k=")))
genes.tib <- genes.tib |>
  group_by(Component) |>
  summarize(TopGenes = paste(SYMBOL, collapse = ", "))

all.gsea.res <- character(30)
for (i in 1:30) {
  cat("Factor", i, "\n")
  gene.list <- ebnmf$fit$L.pm[, i]
  names(gene.list) <- rownames(ebnmf$fit$L.pm)
  gene.list <- sort(gene.list, decreasing = TRUE)[1:100]
  gene.list <- names(gene.list)
  gsea.res <- enrichGO(
    gene.list,
    universe = rownames(ebnmf$fit$L.pm),
    ont = "ALL",
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    readable = TRUE
  )
  top.sets <- gsea.res@result$Description[1]
  all.gsea.res[i] <- top.sets
}
genes.tib <- genes.tib |>
  add_column(GOTerms = all.gsea.res)

tbl <- genes.tib |>
  filter(!(Component %in% c(27, 30))) |>
  dplyr::rename(`Top Genes` = TopGenes,
                `Top GO Term` = GOTerms) |>
  gt() |>
  cols_align("left", columns = Component) |>
  cols_align("left", columns = `Top Genes`) |>
  cols_align("left", columns = `Top GO Term`) |>
  cols_width(
    Component ~ pct(10),
    `Top Genes` ~ pct(55),
    `Top GO Term` ~ pct(35)
  ) |>
  opt_row_striping() |>
  tab_options(table.font.size = 10)

gtsave(tbl, paste0("figs/montoro_topgenes.png"))


k <- 20
tib <- tibble(
  pm = ebnmf$fit$L.pm[, k] / max(abs(ebnmf$fit$L.pm[, k])),
  z = abs(ebnmf$fit$L.pm[, k]) / pmax(sqrt(.Machine$double.eps), ebnmf$fit$L.psd[, k]),
  SYMBOL = rownames(ebnmf$fit$L.pm),
  exprmean = exprmean
) |>
  mutate(expr_bin = factor(ceiling(3 * exprmean))) |>
  group_by(expr_bin) |>
  mutate(bin_quantile = rank(pm) / n()) |>
  # mutate(SYMBOL = ifelse(pm > .125 & bin_quantile > .994, SYMBOL, ""))
  mutate(SYMBOL = ifelse(pm > .125 & 8 * pm - exprmean > 3.45, SYMBOL, ""))

plt <- ggplot(tib, aes(x = exprmean, y = pm, label = SYMBOL)) +
  geom_point() +
  geom_text_repel(color = "darkgray",size = 2.25, fontface = "italic",
                  segment.color = "darkgray", segment.size = 0.25,
                  min.segment.length = 0, na.rm = TRUE, max.overlaps = 20) +
  theme_minimal() +
  labs(
    x = "Mean Expression (log10)",
    y = "Gene Loading (posterior mean)"
  ) +
  theme(
    axis.text.y = element_text(size = 8, family = "serif"),
    axis.text.x = element_text(size = 8, family = "serif"),
    axis.title = element_text(size = 10, margin = margin(2, 2, 2, 2), family = "serif")
  )

ggsave("figs/montoro_ionocyte.png", width = 133, height = 100, units = "mm")

# notable.genes <- c("Cftr", "Foxi1", "Ascl3", "Atp6v0d2")
# plt <- ggplot(tib, aes(x = exprmean, y = pm, label = SYMBOL)) +
#   geom_point() +
#   geom_text_repel(aes(color = SYMBOL %in% notable.genes), size = 4, fontface = "italic",
#                   segment.color = "darkgray", segment.size = 0.25,
#                   min.segment.length = 0, na.rm = TRUE, max.overlaps = 20) +
#   theme_minimal() +
#   labs(
#     x = "Mean Expression (log10)",
#     y = "Gene Loading \n",
#     title = "Component 20 (Ionocytes)"
#   ) +
#   scale_color_manual(values = c("darkgray", "red")) +
#   theme(
#     plot.title = element_text(size = 14, family = "serif", hjust = 0.5),
#     axis.text.y = element_text(size = 10, family = "serif"),
#     axis.text.x = element_text(size = 10, family = "serif"),
#     axis.title = element_text(size = 14, margin = margin(2, 2, 2, 2), family = "serif"),
#     legend.position = "none"
#   )
#
# ggsave("figs/montoro_ionocyte_poster.png", width = 133, height = 100, units = "mm")

nnlm_tib |>
  filter(CellType == "Ionocyte") |>
  ggplot(aes(x = Loading)) +
  geom_histogram(bins = 10) +
  facet_wrap(~Component)

