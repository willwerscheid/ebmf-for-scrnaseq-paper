library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggrepel)
library(gt)

exprmean <- log10(readRDS("data/montoro-mean-expr.rds")$var.gene.mean.expr)

ebnmf <- readRDS("output/montoro-ebmf.rds")
nnlm <- readRDS("output/montoro-nmf.rds")
tm <- readRDS("output/montoro-topics.rds")
tm.col_order <- order(
  apply(tm$fit$F, 2, function(x) sum(x^2)) *
    apply(tm$fit$L, 2, function(x) sum(x^2)),
  decreasing = TRUE
)

fl <- ebnmf$fit
FF <- fl$F.pm
FF <- scale(FF, center = FALSE)

cell.type <- sapply(strsplit(rownames(FF), "_"), `[[`, 5)
time.point <- sapply(strsplit(rownames(FF), "_"), `[[`, 2)

cell.type <- ifelse(
  cell.type == "Basal", paste0("Basal (", time.point, ")"), cell.type
)

cell.type <- factor(cell.type, levels = c(
  "Basal (Tp0)",
  "Basal (Tp30)",
  "Basal (Tp60)",
  "Club",
  "Club (hillock-associated)",
  "Ciliated",
  "Goblet.progenitor",
  "Goblet.1",
  "Goblet.2",
  "Tuft.progenitor",
  "Tuft.1",
  "Tuft.2",
  "Neuroendocrine",
  "Ionocyte",
  "Proliferating"
))

# Downsample the number of cells.
set.seed(666)
cell.idx <- numeric(0)
for (k in levels(cell.type)) {
  which.idx <- which(cell.type == k)
  # Downsample common cell types. Duplicate rare ones.
  which.idx <- sample(which.idx, 200, replace = (length(which.idx) < 200))
  cell.idx <- c(cell.idx, which.idx)
}

cell.type <- cell.type[cell.idx]
ebnmf.FF <- fl$F.pm[cell.idx, ]
nnlm.FF <- (t(nnlm$fit$h))[cell.idx, ]
tm.FF <- tm$fit$F[cell.idx, tm.col_order]

make.heatmap.tib <- function(FF) {
  # Sort cells using t-SNE.
  for (k in levels(cell.type)) {
    which.idx <- which(cell.type == k)
    tsne.res <- Rtsne::Rtsne(
      FF[which.idx, ],
      dims = 1,
      pca = FALSE,
      normalize = FALSE,
      perplexity = min(100, floor((length(which.idx) - 1) / 3) - 1),
      theta = 0.1,
      max_iter = 1000,
      eta = 200,
      check_duplicates = FALSE
    )$Y[, 1]
    FF[which.idx, ] <- FF[which.idx[order(tsne.res)], ]
  }

  colnames(FF) <- 1:ncol(FF)
  tib <- as_tibble(scale(FF, center = FALSE, scale = apply(FF, 2, max))) %>%
    mutate(
      Cell.idx = row_number(),
      Cell.type = cell.type
    )

  tib <- tib %>%
    pivot_longer(
      -c(Cell.idx, Cell.type),
      names_to = "Component",
      values_to = "Loading",
      values_drop_na = TRUE
    ) %>%
    mutate(Component = factor(Component, levels = colnames(FF)))

  return(tib)
}

ebnmf.tib <- make.heatmap.tib(ebnmf.FF)
nnlm.tib <- make.heatmap.tib(nnlm.FF)
tm.tib <- make.heatmap.tib(tm.FF)

heatmap.tib <- ebnmf.tib %>% mutate(Method = "EBNMF") %>%
  bind_rows(nnlm.tib %>% mutate(Method = "NMF")) %>%
  bind_rows(tm.tib %>% mutate(Method = "Topics")) %>%
  mutate(Method = factor(Method, levels = c("EBNMF", "NMF", "Topics")))

tib <- heatmap.tib %>%
  group_by(Cell.type, Cell.idx) %>%
  summarize()
cell_type_breaks <- c(1, which(tib$Cell.type[-1] != tib$Cell.type[-nrow(tib)]))
label_pos <- cell_type_breaks / 2 + c(cell_type_breaks[-1], nrow(tib)) / 2

plt <- ggplot(heatmap.tib, aes(x = Component, y = -Cell.idx, fill = Loading)) +
  geom_tile() +
  scale_fill_gradient2(high = "black") +
  labs(y = "") +
  scale_y_continuous(breaks = -label_pos,
                     minor_breaks = NULL,
                     labels = levels(cell.type)) +
  theme_minimal() +
  geom_hline(yintercept = -cell_type_breaks, linewidth = 0.1) +
  facet_wrap(~Method, scales = "free", ncol = 1) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12, margin = margin(2, 2, 2, 2), family = "serif"),
    axis.text.y = element_text(size = 8, family = "serif"),
    axis.text.x = element_text(size = 6, family = "serif"),
    axis.title = element_text(size = 8, margin = margin(2, 2, 2, 2), family = "serif")
  )

ggsave("figs/montoro_heatmap.png", width = 163, height = 190, units = "mm")

# poster.tib.ebnmf <- heatmap.tib %>%
#   filter(Method == "EBNMF") %>%
#   mutate(Loading = ifelse(Component == "20", -Loading, Loading)) %>%
#   mutate(Cell.type = fct_recode(Cell.type, `Club (hillock)` = "Club (hillock-associated)"))
# plt <- ggplot(poster.tib.ebnmf, aes(x = Component, y = -Cell.idx, fill = Loading)) +
#   geom_tile() +
#   scale_fill_gradient2(low = "red", mid = "white", high = "black") +
#   labs(y = "", title = "EBNMF") +
#   scale_y_continuous(breaks = -label_pos,
#                      minor_breaks = NULL,
#                      labels = levels(poster.tib.ebnmf$Cell.type)) +
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
# poster.tib.nmf <- heatmap.tib %>%
#   filter(Method == "NMF") %>%
#   mutate(Loading = ifelse(Component == "30", -Loading, Loading)) %>%
#   mutate(Cell.type = fct_recode(Cell.type, `Club (hillock)` = "Club (hillock-associated)"))
# plt <- ggplot(poster.tib.nmf, aes(x = Component, y = -Cell.idx, fill = Loading)) +
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
genes.tib <- as_tibble(top.genes) %>%
  pivot_longer(cols = everything(), names_to = "Component", values_to = "SYMBOL") %>%
  mutate(Component = as.numeric(str_remove_all(Component, "k=")))
genes.tib <- genes.tib %>%
  group_by(Component) %>%
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
genes.tib <- genes.tib %>%
  add_column(GOTerms = all.gsea.res)

tbl <- genes.tib %>%
  filter(!(Component %in% c(27, 30))) %>%
  dplyr::rename(`Top Genes` = TopGenes,
                `Top GO Term` = GOTerms) %>%
  gt() %>%
  cols_align("left", columns = Component) %>%
  cols_align("left", columns = `Top Genes`) %>%
  cols_align("left", columns = `Top GO Term`) %>%
  cols_width(
    Component ~ pct(10),
    `Top Genes` ~ pct(55),
    `Top GO Term` ~ pct(35)
  ) %>%
  opt_row_striping() %>%
  tab_options(table.font.size = 10)

gtsave(tbl, paste0("figs/montoro_topgenes.png"))


k <- 20
tib <- tibble(
  pm = ebnmf$fit$L.pm[, k] / max(abs(ebnmf$fit$L.pm[, k])),
  z = abs(ebnmf$fit$L.pm[, k]) / pmax(sqrt(.Machine$double.eps), ebnmf$fit$L.psd[, k]),
  SYMBOL = rownames(ebnmf$fit$L.pm),
  exprmean = exprmean
) %>%
  mutate(expr_bin = factor(ceiling(3 * exprmean))) %>%
  group_by(expr_bin) %>%
  mutate(bin_quantile = rank(pm) / n()) %>%
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

nnlm.tib |>
  filter(Cell.type == "Ionocyte") |>
  ggplot(aes(x = Loading)) +
  geom_histogram(bins = 10) +
  facet_wrap(~Component)

