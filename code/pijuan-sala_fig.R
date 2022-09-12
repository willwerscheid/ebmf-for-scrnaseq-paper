library(tidyverse)
library(gt)
library(Polychrome)
library(ggrepel)

exprmean <- log10(readRDS("data/pijuan-sala-mean-expr.rds")$var.gene.mean.expr)

ebnmf <- readRDS("output/pijuan-sel-ebmf-log-k=40.rds")
nnlm <- readRDS("output/pijuan-sel-nmf-log-k=40.rds")

fl <- ebnmf$fit
LL <- scale(fl$L.pm, center = FALSE, scale = apply(fl$L.pm, 2, max))
FF <- scale(fl$F.pm, center = FALSE, scale = apply(fl$F.pm, 2, max))

cell.type <- sapply(strsplit(rownames(FF), "_"), `[[`, 6)
wild.type <- ifelse(
  sapply(strsplit(rownames(FF), "_"), `[[`, 4),
  "Tal1-knockout",
  "Wild type"
)

cell.type <- factor(cell.type, levels = c(
  "Caudal epiblast",
  "PGC",
  "Notochord",
  "Def. endoderm",
  "Gut",
  "Intermediate mesoderm",
  "Caudal Mesoderm",
  "Paraxial mesoderm",
  "Somitic mesoderm",
  "Pharyngeal mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE mesoderm",
  "Mesenchyme",
  "Haematoendothelial progenitors",
  "Blood progenitors 1",
  "Blood progenitors 2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "Endothelium",
  "NMP",
  "Rostral neurectoderm",
  "Neural crest",
  "Forebrain/Midbrain/Hindbrain",
  "Spinal cord",
  "Surface ectoderm",
  "ExE endoderm",
  "ExE ectoderm",
  "Parietal endoderm",
  "Visceral endoderm",
  "Doublet",
  "Stripped"
))
wild.type <- factor(wild.type, levels = c(
  "Wild type",
  "Tal1-knockout"
))


# Comparison of wt and ko cells for components 12 and 22.

endo1 <- which(cell.type %in% c(
  "Endothelium", "Haematoendothelial progenitors"
))
endo2 <- which(cell.type %in% c(
  "Endothelium", "Blood progenitors 1"
))
endo.tib <- tibble(
  k = "Component 12",
  Loading = FF[endo1, 12],
  Lineage = cell.type[endo1],
  Type = wild.type[endo1]
) %>%
  bind_rows(tibble(
    k = "Component 22",
    Loading = FF[endo2, 22],
    Lineage = cell.type[endo2],
    Type = wild.type[endo2]
  ))

ggplot(endo.tib, aes(x = Loading, color = Lineage)) +
  geom_density(aes(linetype = Type)) +
  labs(
    x = "Cell Loading (posterior mean)",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8, family = "serif"),
    axis.text.x = element_text(size = 8, family = "serif"),
    axis.title = element_text(size = 10, margin = margin(2, 2, 2, 2), family = "serif"),
    strip.text.x = element_text(size = 10, margin = margin(2, 2, 2, 2), family = "serif")
  ) +
  facet_wrap(~k, nrow = 1)

ggsave("figs/pijuan-sala_endothelium.png", width = 100, height = 70, units = "mm")

endo <- which(cell.type %in% c(
  "Endothelium", "Haematoendothelial progenitors", "Blood progenitors 1"
))
endo.tib2 <- tibble(
  Comp12 = FF[endo, 12],
  Comp22 = FF[endo, 22],
  Lineage = cell.type[endo],
  Type = wild.type[endo]
)
ggplot(endo.tib2, aes(x = Comp12, y = Comp22, color = Type, shape = Lineage)) +
  geom_point(size = 1)


# Downsample the number of cells. Sort them using t-SNE.
set.seed(666)
cell.idx <- numeric(0)
for (k in levels(cell.type)) {
  which.idx <- which(cell.type == k)
  # Downsample common cell types. Duplicate rare ones.
  which.idx <- sample(which.idx, 200, replace = (length(which.idx) < 200))
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
  cell.idx <- c(cell.idx, which.idx[order(tsne.res)])
}

cell.type <- cell.type[cell.idx]
wild.type <- wild.type[cell.idx]
ebnmf.FF <- fl$F.pm[cell.idx, ]
nnlm.FF <- (t(nnlm$fit$h))[cell.idx, ]
# tm.FF <- tm$fit$F[cell.idx, ]

make.heatmap.tib <- function(FF) {
  colnames(FF) <- 1:ncol(FF)
  tib <- as_tibble(scale(FF, center = FALSE, scale = apply(FF, 2, max))) %>%
    mutate(
      Cell.type = cell.type,
      Wild.type = wild.type
    ) %>%
    arrange(Cell.type, Wild.type) %>%
    mutate(Cell.idx = row_number())

  tib <- tib %>%
    pivot_longer(
      -c(Cell.idx, Cell.type, Wild.type),
      names_to = "Component",
      values_to = "Loading",
      values_drop_na = TRUE
    ) %>%
    mutate(Component = factor(Component, levels = colnames(FF)))

  return(tib)
}

ebnmf.tib <- make.heatmap.tib(ebnmf.FF)
nnlm.tib <- make.heatmap.tib(nnlm.FF)
# tm.tib <- make.heatmap.tib(tm.FF)

heatmap.tib <- ebnmf.tib %>% mutate(Method = "EBNMF") %>%
  bind_rows(nnlm.tib %>% mutate(Method = "NMF")) %>%
  # bind_rows(tm.tib %>% mutate(Method = "Topics")) %>%
  # mutate(Method = factor(Method, levels = c("EBNMF", "NMF", "Topics")))
  mutate(Method = factor(Method, levels = c("EBNMF", "NMF")))

tib <- heatmap.tib %>%
  group_by(Cell.type, Wild.type, Cell.idx) %>%
  summarize()
cell_type_breaks <- c(1, which(tib$Cell.type[-1] != tib$Cell.type[-nrow(tib)]))
label_pos <- cell_type_breaks / 2 + c(cell_type_breaks[-1], nrow(tib)) / 2
wt_breaks <- which(tib$Wild.type[-1] != tib$Wild.type[-nrow(tib)])
wt_breaks <- setdiff(wt_breaks, cell_type_breaks)

ggplot(heatmap.tib, aes(x = Component, y = -Cell.idx, fill = Loading)) +
  geom_tile() +
  scale_fill_gradient2(high = "black") +
  labs(y = "") +
  scale_y_continuous(breaks = -label_pos,
                     minor_breaks = NULL,
                     labels = levels(cell.type)) +
  theme_minimal() +
  geom_hline(yintercept = -cell_type_breaks, size = 0.1) +
  geom_hline(yintercept = -wt_breaks, size = 0.1, linetype = "dashed") +
  facet_wrap(~Method, scales = "free", ncol = 1) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12, margin = margin(2, 2, 2, 2), family = "serif"),
    axis.text.y = element_text(size = 8, family = "serif"),
    axis.text.x = element_text(size = 6, family = "serif"),
    axis.title = element_text(size = 8, margin = margin(2, 2, 2, 2), family = "serif")
  )

ggsave("figs/pijuan-sala_heatmap.png", width = 175, height = 190, units = "mm")


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

tbl <- genes.tib %>%
  filter(!(Component %in% c(27, 31, 34, 37, 39))) %>%
  rename(`Top Genes` = TopGenes) %>%
  gt() %>%
  cols_align("left", columns = Component) %>%
  cols_align("left", columns = `Top Genes`) %>%
  cols_width(
    Component ~ pct(12),
    `Top Genes` ~ pct(88)
  ) %>%
  opt_row_striping()

gtsave(tbl, paste0("figs/pijuan-sala_topgenes.png"))





make_gene_tib <- function(k, lbl_slope, lbl_intercept) {
  tib <- tibble(
    comp = paste("Component", k),
    pm = ebnmf$fit$L.pm[, k] / max(abs(ebnmf$fit$L.pm[, k])),
    SYMBOL = rownames(ebnmf$fit$L.pm),
    exprmean = exprmean
  ) %>%
    mutate(disp_SYMBOL = ifelse(pm > .125 & exprmean < lbl_slope * pm + lbl_intercept, SYMBOL, ""))

  return(tib)
}
tib <- make_gene_tib(12, 10, -3.5) %>%
  bind_rows(make_gene_tib(22, 9, -3.6))

diff_genes <- tib %>%
  select(-disp_SYMBOL) %>%
  pivot_wider(names_from = comp, values_from = c(pm)) %>%
  mutate(diff = `Component 12`/`Component 22`) %>%
  filter(diff > 4 | diff < 0.25) %>%
  pull(SYMBOL)
tib <- tib %>%
  mutate(lbl_color = ifelse(disp_SYMBOL %in% diff_genes, "red", "darkgray"))

plt <- ggplot(tib, aes(x = pm, y = exprmean, label = disp_SYMBOL)) +
  geom_point() +
  geom_text_repel(color = tib$lbl_color, size = 2, fontface = "italic",
                  segment.color = "darkgray", segment.size = 0.25,
                  min.segment.length = 0, na.rm = TRUE, max.overlaps = 20) +
  theme_minimal() +
  labs(
    x = "Gene Loading (posterior mean)",
    y = "Mean Expression (log10)",
  ) +
  theme(
    axis.text.y = element_text(size = 8, family = "serif"),
    axis.text.x = element_text(size = 8, family = "serif"),
    axis.title = element_text(size = 10, margin = margin(2, 2, 2, 2), family = "serif")
  ) +
  facet_wrap(~comp, ncol = 1)

ggsave("figs/pijuan-sala_comp.png", width = 133, height = 200, units = "mm")

