library(clusterProfiler)
library(org.Mm.eg.db)
library(ggrepel)
library(gt)
library(tidyverse)
library(cowplot)

# exprmean <- log10(readRDS("data/montoro-mean-expr.rds")$var.gene.mean.expr)

ebnmf <- readRDS("output/montoro-ebmf.rds")
nnlm <- readRDS("output/montoro-nmf.rds")
tm <- readRDS("output/montoro-topics.rds")

ebnmf_FF <- ebnmf$fit$F.pm
nnlm_FF <- t(nnlm$fit$h)
tm_FF <- tm$fit$F

calc_factor_sizes <- function(FF, LL, D = 1) {
  l1_norm <- apply(FF, 2, max)
  LL <- t(t(LL) * l1_norm * D)
  return(apply(LL, 2, function(x) sqrt(sum(x^2))))
}

ebnmf_widths <- calc_factor_sizes(ebnmf_FF, ebnmf$fit$L.pm)
nnlm_widths <- calc_factor_sizes(nnlm_FF, nnlm$fit$w, nnlm$fit$d)
# Topic models are fit to untransformed counts:
tm_widths <- log10(calc_factor_sizes(tm_FF, tm$fit$L))

external_info <- tibble(
  CellType = sapply(strsplit(rownames(ebnmf_FF), "_"), `[[`, 5),
  TimePoint = sapply(strsplit(rownames(ebnmf_FF), "_"), `[[`, 2),
  Replicate = sapply(strsplit(rownames(ebnmf_FF), "_"), `[[`, 3),
  Fluor = sapply(strsplit(rownames(ebnmf_FF), "_"), `[[`, 4)
)
cell_type <- external_info$CellType

# big_types <- c("Basal", "Club")
big_types <- "Basal"
cell_type <- ifelse(
  cell_type %in% big_types, paste0(
    cell_type, " ", external_info$TimePoint, " (", external_info$Replicate, ")"
  ),
  cell_type
)
lvls <- unique(cell_type)
cell_type <- factor(cell_type, levels = c(
  lvls[str_detect(lvls, "Basal Tp")],
  "Proliferating",
  # lvls[str_detect(lvls, "Club Tp")],
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
  "Ionocyte"
))

# Downsample the number of cells.
set.seed(666)
cell_idx <- numeric(0)
for (k in levels(cell_type)) {
  which_idx <- which(cell_type == k)
  # Downsample common cell types. Duplicate rare ones.
  # cat(k, " ", length(which_idx), "\n")
  which_idx <- sample(which_idx, 200, replace = (length(which_idx) < 200))
  cell_idx <- c(cell_idx, which_idx)
}
external_info <- slice(external_info, cell_idx)
cell_type <- cell_type[cell_idx]
ebnmf_FF <- ebnmf_FF[cell_idx, ]
nnlm_FF <- nnlm_FF[cell_idx, ]
tm_FF <- tm_FF[cell_idx, ]

make_heatmap_tib <- function(FF, widths, type_thresh = 0.2, r2_thresh = 0.15, intercept = TRUE) {
  # Sort cells using t-SNE.
  for (k in levels(cell_type)) {
    which_idx <- which(cell_type == k)
    tsne_res <- Rtsne::Rtsne(
      FF[which_idx, ],
      dims = 1,
      pca = FALSE,
      normalize = FALSE,
      perplexity = min(100, floor((length(which_idx) - 1) / 3) - 1),
      theta = 0.1,
      max_iter = 1000,
      eta = 200,
      check_duplicates = FALSE
    )$Y[, 1]
    FF[which_idx, ] <- FF[which_idx[order(tsne_res)], ]
  }

  colnames(FF) <- 1:ncol(FF)
  tib <- as_tibble(scale(FF, center = FALSE, scale = apply(FF, 2, max))) |>
    mutate(
      CellIdx = row_number(),
      CellType = cell_type
    )

  tib <- tib |>
    pivot_longer(
      -c(CellIdx, CellType),
      names_to = "Component",
      values_to = "Loading",
      values_drop_na = TRUE
    ) |>
    mutate(Component = factor(Component, levels = colnames(FF)))

  all_coef <- NULL
  all_r2 <- NULL
  for (k in 1:ncol(FF)) {
    y <- FF[, k] / max(FF[, k])
    lmod <- lm(FF[, k] ~ TimePoint + Replicate + Fluor + CellType,
               data = external_info)
    all_coef <- cbind(all_coef, coef(lmod)[-1])
    all_r2 <- c(all_r2, broom::glance(lmod)$r.squared)
  }
  all_coef <- as_tibble(all_coef) |>
    mutate(Type = str_extract(rownames(all_coef), "TimePoint|Replicate|Fluor|CellType"))
  coef_sum <- all_coef |>
    group_by(Type) |>
    summarize(across(everything(), ~max(abs(.x)))) |>
    mutate(across(-1, ~ round(.x / max(.x), 2)))
  coef_sum <- coef_sum |>
    pivot_longer(-Type, names_to = "Component", values_to = "TypeWt")
  coef_sum <- coef_sum |>
    bind_rows(tibble(
      Type = "None", Component = paste0("V", 1:ncol(FF)), TypeWt = 0
    ))
  component_types <- coef_sum |>
    mutate(Type = factor(Type, levels = c("CellType", "TimePoint", "Replicate", "Fluor", "None")),
           Component = factor(parse_number(Component)),
           r2 = rep(all_r2, times = 5)) |>
    filter((r2 > r2_thresh & TypeWt > type_thresh) | (r2 <= r2_thresh & Type == "None")) |>
    mutate(Width = widths[Component])

  get_subtype <- function(type, k) {
    tib |>
      filter(Component == k) |>
      mutate(Val = case_when(
        type == "TimePoint" ~ str_extract(CellType, "Tp[0-9]+"),
        type == "Replicate" ~ str_extract(CellType, "R[0-9]"),
        type == "CellType" ~ str_remove(CellType, " Tp.*"),
        TRUE ~ "None"
      )) |>
      filter(!is.na(Val)) |>
      group_by(Component, Val) |>
      summarize(Loading = mean(Loading), .groups = "drop_last") |>
      mutate(MeanLoading = mean(Loading)) |>
      filter(Loading == max(Loading)) |>
      pull(Val)
  }

  component_types <- component_types |>
    mutate(SubType = mapply(get_subtype, Type, Component)) |>
    mutate(SubType = factor(SubType, levels = c(
      unique(str_remove(levels(tib$CellType), " Tp.*")),
      unique(str_extract(levels(tib$CellType), "Tp[0-9]+")),
      unique(str_extract(levels(tib$CellType), "R[0-9]")),
      "None"
    )))

  if (intercept) {
    sorted_components <- component_types |>
      arrange(desc(Component == 1), SubType, desc(Width))
  } else {
    sorted_components <- component_types |>
      arrange(SubType, desc(Width))
  }

  sorted_components <- sorted_components |>
    mutate(SortedComponent = as_factor(1:n()),
           x_end = cumsum(Width),
           x_pos = x_end - Width / 2) |>
    select(Component, SortedComponent, Width, x_end, x_pos)

  tib <- tib |>
    left_join(sorted_components, by = "Component")

  return(tib)
}

ebnmf_tib <- make_heatmap_tib(ebnmf_FF, ebnmf_widths)
nnlm_tib <- make_heatmap_tib(nnlm_FF, nnlm_widths)
tm_tib <- make_heatmap_tib(tm_FF, tm_widths, intercept = FALSE)

heatmap_tib <- ebnmf_tib |> mutate(Method = "EBNMF") |>
  bind_rows(nnlm_tib |> mutate(Method = "NMF")) |>
  bind_rows(tm_tib |> mutate(Method = "Topics")) |>
  mutate(Method = factor(Method, levels = c("EBNMF", "NMF", "Topics")))

tib <- heatmap_tib |>
  group_by(CellType, CellIdx) |>
  summarize()
cell_type_breaks <- c(1, which(tib$CellType[-1] != tib$CellType[-nrow(tib)]))
label_pos <- cell_type_breaks / 2 + c(cell_type_breaks[-1], nrow(tib)) / 2

max_idx <- max(heatmap_tib$CellIdx)

p1 <- ggplot(heatmap_tib,
              aes(x = x_pos, y = -CellIdx, fill = Loading)) +
  geom_tile(aes(width = Width)) +
  scale_fill_gradient2(high = "black") +
  labs(x = "", y = "") +
  scale_y_continuous(breaks = -label_pos,
                     minor_breaks = NULL,
                     labels = levels(cell_type)) +
  theme_minimal() +
  geom_hline(yintercept = -cell_type_breaks, linewidth = 0.1) +
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

