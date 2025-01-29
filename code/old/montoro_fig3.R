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

external_info <- rownames(ebnmf_FF)
external_info <- str_split(external_info, "_")
external_info <- tibble(
  TimePoint = factor(sapply(external_info, `[[`, 2)),
  Replicate = factor(sapply(external_info, `[[`, 3)),
  Fluorescence = factor(sapply(external_info, `[[`, 4)),
  CellType = factor(sapply(external_info, `[[`, 5))
)
external_info <- external_info |>
  mutate(CellType = factor(CellType, levels = c(
    "Basal",
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
  )))

calc_factor_sizes <- function(FF, LL, D = 1) {
  l1_norm <- apply(FF, 2, max)
  LL <- t(t(LL) * l1_norm * D)
  return(apply(LL, 2, function(x) sqrt(sum(x^2))))
}

ebnmf_widths <- calc_factor_sizes(ebnmf_FF, ebnmf$fit$L.pm)
nnlm_widths <- calc_factor_sizes(nnlm_FF, nnlm$fit$w, nnlm$fit$d)
# Topic models are fit to untransformed counts:
tm_widths <- log10(calc_factor_sizes(tm_FF, tm$fit$L))

# Downsample the number of cells.
set.seed(666)
cell_idx <- numeric(0)
cell_type <- pull(external_info, CellType)
for (k in levels(cell_type)) {
  which_idx <- which(cell_type == k)
  # Downsample common cell types. Duplicate rare ones.
  which_idx <- sample(which_idx, 200, replace = (length(which_idx) < 200))
  cell_idx <- c(cell_idx, which_idx)
}
external_info <- slice(external_info, cell_idx)
ebnmf_FF <- ebnmf_FF[cell_idx, ]
nnlm_FF <- nnlm_FF[cell_idx, ]
tm_FF <- tm_FF[cell_idx, ]

make_heatmap_tib <- function(FF, widths, ftype_thresh = 0.5, intercept = TRUE) {
  next_coef <- function(k) {
    lmod <- lm(FF[, k] ~ TimePoint + Replicate + Fluorescence + CellType,
               data = external_info)
    return(coef(lmod)[-1])
  }
  all_coef <- sapply(1:30, next_coef)
  all_coef <- as_tibble(all_coef) |>
    mutate(Type = str_extract(rownames(all_coef), "TimePoint|Replicate|Fluorescence|CellType"))
  coef_sum <- all_coef |>
    group_by(Type) |>
    summarize(across(everything(), ~max(abs(.x)))) |>
    mutate(across(-1, ~ round(.x / max(.x), 2)))

  coef_mat <- coef_sum |>
    column_to_rownames("Type") |>
    as.matrix()

  ftype_list <- apply(coef_mat, 1, function(x) which(x > ftype_thresh))
  all_tibs <- NULL
  for (ftype in names(ftype_list)) {
    info <- pull(external_info[, ftype])
    which_k <- ftype_list[[ftype]]

    if (length(which_k) > 0) {
      # Sort cells using t-SNE.
      sorted_FF <- NULL
      all_idx <- NULL
      for (i in levels(info)) {
        which_idx <- which(info == i)
        all_idx <- c(all_idx, which_idx)
        tsne_res <- Rtsne::Rtsne(
          FF[which_idx, which_k, drop = FALSE],
          dims = 1,
          pca = FALSE,
          normalize = FALSE,
          perplexity = min(100, floor((length(which_idx) - 1) / 3) - 1),
          theta = 0.1,
          max_iter = 1000,
          eta = 200,
          check_duplicates = FALSE
        )$Y[, 1]
        sorted_FF <- rbind(sorted_FF, FF[which_idx[order(tsne_res)], which_k, drop = FALSE])
      }

      colnames(sorted_FF) <- which_k
      next_tib <- as_tibble(scale(sorted_FF, center = FALSE, scale = apply(sorted_FF, 2, max))) |>
        mutate(
          CellIdx = row_number(),
          Info = info[all_idx]
        )

      next_tib <- next_tib |>
        pivot_longer(
          -c(CellIdx, Info),
          names_to = "Component",
          values_to = "Loading",
          values_drop_na = TRUE
        ) |>
        mutate(Component = factor(Component, levels = which_k))

      sorted_components <- next_tib |>
        group_by(Component, Info) |>
        summarize(Loading = mean(Loading), .groups = "drop_last") |>
        mutate(MeanLoading = mean(Loading)) |>
        filter(Loading == max(Loading)) |>
        ungroup() |>
        mutate(Width = widths[which_k])

      if (intercept) {
        sorted_components <- sorted_components |>
          arrange(desc(Component == "1"), Info, desc(Loading))
      } else {
        sorted_components <- sorted_components |>
          arrange(Info, desc(Loading))
      }

      sorted_components <- sorted_components |>
        mutate(SortedComponent = as_factor(1:n()),
               x_end = cumsum(Width),
               x_pos = x_end - Width / 2) |>
        select(Component, SortedComponent, Width, x_end, x_pos)

      next_tib <- next_tib |>
        left_join(sorted_components, by = "Component")

      all_tibs <- all_tibs |>
        bind_rows(all_tibs, next_tib |> mutate(InfoType = ftype))
    }
  }

  return(all_tibs)
}

ebnmf_tib <- make_heatmap_tib(ebnmf_FF, ebnmf_widths, ftype_thresh = 0.5)
nnlm_tib <- make_heatmap_tib(nnlm_FF, nnlm_widths, ftype_thresh = 0.5)
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
    aes(x = x_pos, y = -1.05 * max_idx, label = SortedComponent),
    size = 2
  )

p1

make_plot <- function(tib, ftype) {
  next_tib <- tib |>
    filter(InfoType == ftype) |>
    mutate(Info = fct_drop(Info))

  next_cells <- next_tib |>
    group_by(Info, CellIdx) |>
    summarize()
  info_breaks <- c(1, which(next_cells$Info[-1] != next_cells$Info[-nrow(next_cells)]))
  label_pos <- info_breaks / 2 + c(info_breaks[-1], nrow(next_cells)) / 2
  max_idx <- max(next_tib$CellIdx)

  p <- ggplot(next_tib,
               aes(x = x_pos, y = -CellIdx, fill = Loading)) +
    geom_tile(aes(width = Width)) +
    scale_fill_gradient2(high = "black") +
    labs(x = "", y = "") +
    scale_y_continuous(breaks = -label_pos,
                       minor_breaks = NULL,
                       labels = levels(next_tib$Info)) +
    theme_minimal() +
    geom_hline(yintercept = -info_breaks, linewidth = 0.1) +
    #facet_wrap(~Method, scales = "free", ncol = 1) +
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
      data = filter(next_tib, CellIdx == 1),
      aes(x = x_pos, y = -1.05 * max_idx, label = Component),
      size = 2
    )
}
