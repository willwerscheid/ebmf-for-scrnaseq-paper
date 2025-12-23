library(tidyverse)
library(cowplot)
library(ggrepel)
library(reactable)
library(htmltools)

scale_FF <- function(FF, LL, D = 1) {
  LL_norms <- apply(LL, 2, function(x) sqrt(sum(x^2)))
  return(t(t(FF) * D * LL_norms))
}

ebnmf <- readRDS("output/montoro-ebmf.rds")
ebnmf_FF <- scale_FF(ebnmf$fit$F_pm, ebnmf$fit$L_pm)

cell_types <- c(
  "Ciliated",
  "Proliferating",
  "Basal",
  "Club",
  "Club (hillock-associated)",
  "Goblet.progenitor",
  "Goblet.1",
  "Goblet.2",
  "Tuft.progenitor",
  "Tuft.1",
  "Tuft.2",
  "Neuroendocrine",
  "Ionocyte"
)

ct_abbr <- cell_types
ct_abbr <- str_replace(ct_abbr, "\\.", "\\-")
ct_abbr[2] <- "Prolif."
ct_abbr[5] <- "Hillock"
ct_abbr[12] <- "NEC"
ct_abbr <- str_replace(ct_abbr, "progenitor", "Pr.")

external_info <- tibble(
  CellType = sapply(strsplit(rownames(ebnmf_FF), "_"), `[[`, 5),
  TimePoint = sapply(strsplit(rownames(ebnmf_FF), "_"), `[[`, 2),
  Replicate = sapply(strsplit(rownames(ebnmf_FF), "_"), `[[`, 3)
) |>
  left_join(tibble(CellType = cell_types, CtAbbr = ct_abbr), by = "CellType") |>
  mutate(
    CellType = factor(CellType, levels = cell_types),
    CtAbbr = factor(CtAbbr, levels = ct_abbr),
    TimePoint = factor(TimePoint),
    Replicate = factor(Replicate),
  ) |>
  mutate(
    TPxRep = fct_cross(TimePoint, Replicate),
    CommonType = fct_collapse(CtAbbr, Other = ct_abbr[6:13])
  )

# Make data frames. -----

make_tib <- function(FF, info_col, ksort, ncells = 50, duprate = 1) {
  info <- external_info[[info_col]]

  # Downsample the number of cells.
  set.seed(666)
  cell_idx <- numeric(0)
  for (k in levels(info)) {
    which_idx <- which(info == k)
    # Downsample common cell types. Duplicate rare ones.
    if (ncells < length(which_idx)) {
      which_idx <- sample(which_idx, ncells, replace = FALSE)
    } else if (ncells >= length(which_idx) * duprate) {
      which_idx <- rep(which_idx, duprate)
    }
    cell_idx <- c(cell_idx, which_idx)
  }

  FF <- FF[cell_idx, ]
  info <- info[cell_idx]
  colnames(FF) <- 1:ncol(FF)

  lvls <- as.numeric(ksort)
  if (any(setdiff(1:ncol(FF), ksort))) {
    lvls <- c("Other", lvls)
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
    mutate(
      Component = ifelse(as.numeric(Component) %in% ksort, Component, "Other")
    ) |>
    group_by(CellIdx, Info, Component) |>
    summarize(Loading = sum(Loading), .groups = "drop") |>
    mutate(
      Component = factor(Component, levels = lvls)
    ) |>
    arrange(Info, Component)

  levels(tib$Component) <- c("Other", names(ksort))

  return(tib)
}


# # As a first pass:
# all_sums <- list()
# for (k in 1:30) {
#   all_sums[[k]] <- summary(
#     lm(ebnmf_FF[, k] ~ external_info$CellType + external_info$TimePoint + external_info$Replicate)
#   )
# }


## EBNMF ------

ebnmf_ksort <- c(
  26, 27, 30, # Sparse
  1, 8, 14, 21, # Dense
  6, 11, 2, 23, # Time points
  5, 18, 25, # Ciliated
  9, 20, 22, # Proliferating
  4, 7, # Basal
  3, 13, 10, 16, # 26, Club-Hillock
  24, 28, 29, # Goblet
  15, 17, 12, 19 # Tuft-NEC-Ionocyte
)
names(ebnmf_ksort) <- c(
  "Sprs-a", "Sprs-b", "Sprs-c",
  "Dns-a", "Dns-b", "Dns-c", "Dns-d",
  "Tp0", "Tp30", "Tp60", "Tp30xR1",
  "Cil-a", "Cil-b", "Cil-c",
  "Prl-a", "Prl-b", "Prl-c",
  "Bas", "B/Cl",
  "Cl/H/G1", "Cl/H/G2", "Hil-a", "Hil-b",
  "GobP/1", "Gob1", "Gob2",
  "Tuf1", "Tuf2", "Nec", "Ion"
)

ebnmf_sparse_k <- 1:3
ebnmf_dense_k <- 4:7
ebnmf_tp_k <- 8:11
ebnmf_common_k <- c(12:14, 18:19, 22:23)
ebnmf_rare_k <- c(15:17, 20:21, 24:30)

ebnmf_sparse <- make_tib(ebnmf_FF, "CommonType", ebnmf_ksort[ebnmf_sparse_k], Inf, 1)
ebnmf_dense <- make_tib(ebnmf_FF, "CommonType", ebnmf_ksort[ebnmf_dense_k], Inf, 1)
ebnmf_tp <- make_tib(ebnmf_FF, "TPxRep", ebnmf_ksort[ebnmf_tp_k], Inf, 1)
ebnmf_common <- make_tib(ebnmf_FF, "CommonType", ebnmf_ksort[ebnmf_common_k], Inf, 1)
ebnmf_rare <- make_tib(ebnmf_FF, "CtAbbr", ebnmf_ksort[ebnmf_rare_k], 100, 3)


## NMF -----

nmf <- readRDS("output/montoro-nmf.rds")
nmf_FF <- scale_FF(t(nmf$fit$h), nmf$fit$w, nmf$fit$d)

nmf_ksort <- c(
  1, 2, 24, 27, 30, # dense
  6, 15, 3, 13, 19, # time points
  14, 23, # ciliated
  9, 18, 28, # prolif
  4, 5, 7, 8, 10, 16, # basal
  12, 20, 17, 22, # club
  11, 25, # club/goblet
  21, # hillock
  29, # tuft
  26 # NEC/ionocyte
)
names(nmf_ksort) <- c(
  "Dns-a", "Dns-b", "Dns-c", "Dns-d", "Dns-e",
  "Tp0", "Tp0/30", "Tp30", "Tp60-a", "Tp60-b",
  "Cil-a", "Cil-b",
  "Prl-a", "Prl-b", "Prl-c",
  "Bas-a", "Bas-b", "Bas-c", "Bas-d", "Bas-e", "Bas-f",
  "B/Cl", "B/H", "Cl/H", "B/Cl/H",
  "Cl/G-a", "Cl/G-b",
  "Hil",
  "Tuf",
  "N/I"
)

nmf_dense_k <- 1:5
nmf_tp_k <- 6:10
nmf_basal_k <- 16:21
nmf_common_k <- c(11:12, 22:25, 28)
nmf_rare_k <- c(13:15, 26:27, 29:30)

nmf_dense <- make_tib(nmf_FF, "CommonType", nmf_ksort[nmf_dense_k], Inf, 1)
nmf_tp <- make_tib(nmf_FF, "TPxRep", nmf_ksort[nmf_tp_k], Inf, 1)
nmf_basal <- make_tib(nmf_FF, "CommonType", nmf_ksort[nmf_basal_k], Inf, 1)
nmf_common <- make_tib(nmf_FF, "CommonType", nmf_ksort[nmf_common_k], Inf, 1)
nmf_rare <- make_tib(nmf_FF, "CtAbbr", nmf_ksort[nmf_rare_k], 100, 3)


## Topic modeling. -----

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

tm_ksort <- c(
  7, 9, 11, # Dense
  14, 15, 5, 6, 12, 20, 8, # Time points
  1, # Ciliated
  21, 30, # Proliferating
  22, 28, 18, 27, # Basal, Club
  16, 19, # Club/Goblet
  23, 29, 3, 4, # Hillock
  10, 17, 26, # Goblet
  2, 24, 13, 15 # Tuft-NEC-Ion
)
names(tm_ksort) <- c(
  "Dns-a", "Dns-b", "Dns-c",
  "Tp0-a", "Tp0-b", "Tp0/60", "Tp30", "Tp60-a", "Tp60-b", "Tp30xR1",
  "Cil",
  "Prl-a", "Prl-b",
  "Bas-a", "Bas-b", "Cl-a", "Cl-b",
  "Cl/H/G-a", "Cl/H/G-b",
  "Hil", "B/Cl/H", "H/G", "Hil",
  "GobP/1-a", "GobP/1-b", "GobP/1-c",
  "T/N-a", "T/N-b", "Ion", "P/B/I"
)

tm_dense_k <- 1:3
tm_tp_k <- 4:10
tm_common_k <- c(11, 14:17, 20:21, 23)
tm_rare_k <- c(12:13, 18:19, 22, 24:30)

tm_dense <- make_tib(tm_FF, "CommonType", tm_ksort[tm_dense_k], Inf, 1)
tm_tp <- make_tib(tm_FF, "TPxRep", tm_ksort[tm_tp_k], Inf, 1)
tm_common <- make_tib(tm_FF, "CommonType", tm_ksort[tm_common_k], Inf, 1)
tm_rare <- make_tib(tm_FF, "CtAbbr", tm_ksort[tm_rare_k], 100, 3)


# Make plots. -----

do_plot <- function(tib, colors, ptitle, ntiles = 400) {
  struct_df <- tib |>
    pivot_wider(names_from = Component, values_from = Loading)
  Lmat <- struct_df |>
    select(-(1:2)) |>
    as.matrix()
  group <- struct_df$Info

  if ("Other" %in% tib$Component) {
    colors <- c("grey90", colors)
  }

  p <- fastTopics::structure_plot(
    Lmat, topics = colnames(Lmat),
    grouping = group,
    colors = colors,
    # embed_method = fastTopics::umap_from_topics,
    # dims = 1,
    gap = 10,
    verbose = FALSE
  ) + theme(
    axis.text.y = element_blank()
  ) + labs(
    y = "",
    title = ptitle
  )

  return(p)
}

sparse_colors <- fastTopics:::glasbey()[2:4]
dense_colors <- fastTopics:::glasbey()[5:13]
tp_colors <- fastTopics:::glasbey()[14:20]
common_colors <- fastTopics:::glasbey()[21:33]
rare_colors <- fastTopics:::glasbey()[34:49]

# Colors: 1-3 proliferating; 4-6 club/hillock; 7-9 goblet 1; 10 goblet;
#   11-12 tuft; 13-14 NEC; 15-16 ionocyte
ebnmf_rare_p <- do_plot(ebnmf_rare, rare_colors[c(1:3, 4:5, 7:8, 10, 11:12, 13, 15)], "EBNMF")
nmf_rare_p <- do_plot(nmf_rare, rare_colors[c(1:3, 4:5, 11, 13)], "Vanilla NMF")
tm_rare_p <- do_plot(tm_rare, rare_colors[c(1:2, 4:6, 7:9, 13:14, 15:16)], "Topic Model")
plot_grid(ebnmf_rare_p, nmf_rare_p, tm_rare_p, nrow = 3, ncol = 1)

# sparse_title <- "Very sparse factors"
# dense_title <- "Dense factors"
# tp_title <- "Time point factors"
# common_title <- "Prevalent cell type factors"
# rare_title <- "Rare cell type factors (common cell types are downsampled)"

# sparse_p <- do_plot(ebnmf_sparse, glasbey(sparse_colors), sparse_title)
# dense_p <- do_plot(ebnmf_dense, glasbey(dense_colors), dense_title)
# tp_p <- do_plot(ebnmf_tp, glasbey(tp_colors), tp_title)
# common_p <- do_plot(ebnmf_common, glasbey(common_colors), common_title)
# rare_p <- do_plot(ebnmf_rare, glasbey(rare_colors), rare_title)
# plot_grid(common_p, rare_p, tp_p, dense_p, sparse_p, nrow = 5, ncol = 1)

# nmf_dense_p <- do_plot(nmf_tib, glasbey(dense_colors), dense_title)
# nmf_tp_p <- do_plot(nmf_tib, glasbey(tp_colors), tp_title)
# nmf_basal_p <- do_plot(nmf_tib, glasbey(common_colors[1:6]), "Basal cell factors")
# nmf_common_p <- do_plot(nmf_tib, glasbey(common_colors[-(1:6)]), paste("Other", tolower(common_title)))
# nmf_rare_p <- do_plot(nmf_tib, glasbey(rare_colors), rare_title)
# plot_grid(nmf_basal_p, nmf_common_p, nmf_rare_p, nmf_tp_p, nmf_dense_p, nrow = 5, ncol = 1)

# tm_dense_p <- do_plot(tm_tib, glasbey(dense_colors), dense_title)
# tm_tp_p <- do_plot(tm_tib, glasbey(tp_colors), tp_title)
# tm_common_p <- do_plot(tm_tib, glasbey(common_colors), common_title)
# tm_rare_p <- do_plot(tm_tib, glasbey(rare_colors), rare_title)
# plot_grid(tm_common_p, tm_rare_p, tm_tp_p, tm_dense_p, nrow = 4, ncol = 1)

# Gene sets. -----

make_gene_set_tbl <- function(GSEA, ksort) {
  klabs <- tibble(
    korig = sort(ksort),
    ksort = ksort,
    label = factor(names(ksort), levels = names(ksort))
  )

  keep_res <- GSEA |>
    filter(enrichmentScore > 0) |>
    left_join(klabs, by = c("k" = "ksort")) |>
    group_by(k) |>
    filter(FDR < 0.05 | rank(FDR) == min(rank(FDR))) |>
    ungroup()

  shared_genesets <- keep_res |>
    left_join(keep_res |> select(geneSet, shared_by = label),
              by = "geneSet",
              relationship = "many-to-many") |>
    filter(shared_by != label) |>
    arrange(label) |>
    group_by(label, geneSet) |>
    summarize(shared_by = paste(shared_by, collapse = ", "))

  keep_res <- keep_res |>
    select(label, geneSet, description, FDR) |>
    left_join(shared_genesets, by = c("label", "geneSet")) |>
    arrange(label, FDR, description)

  keep_res <- keep_res |>
    mutate(FDR = ifelse(FDR < 0.01, "<.01", format(round(FDR, 2), nsmall = 2)))

  react_tbl <- reactable(
    keep_res,
    columns = list(
      label = colDef(name = "Component"),
      geneSet = colDef(name = "Gene Set"),
      description = colDef(name = "Description"),
      shared_by = colDef(name = "Shared By")
    ),
    groupBy = "label"
  )

  return(react_tbl)
}

bar_chart <- function(width = "100%", height = "1rem", fill = "#00bfc4", background = NULL) {
  bar <- div(style = list(background = fill, width = width, height = height))
  chart <- div(style = list(flexGrow = 1, marginLeft = "0.5rem", background = background), bar)
  div(style = list(display = "flex", alignItems = "center"), chart)
}

bar_chart_pos_neg <- function(label, value, max_value = 1, height = "1rem",
                              pos_fill = "#005ab5", neg_fill = "#dc3220") {
  neg_chart <- div(style = list(flex = "1 1 0"))
  pos_chart <- div(style = list(flex = "1 1 0"))
  width <- paste0(abs(value / max_value) * 100, "%")

  if (value < 0) {
    bar <- div(style = list(marginLeft = "0.5rem", background = neg_fill, width = width, height = height))
    chart <- div(
      style = list(display = "flex", alignItems = "center", justifyContent = "flex-end"),
      label,
      bar
    )
    neg_chart <- tagAppendChild(neg_chart, chart)
  } else {
    bar <- div(style = list(marginRight = "0.5rem", background = pos_fill, width = width, height = height))
    chart <- div(style = list(display = "flex", alignItems = "center"), bar, label)
    pos_chart <- tagAppendChild(pos_chart, chart)
  }

  div(style = list(display = "flex"), neg_chart, pos_chart)
}

make_top_genes_tbl <- function(LL, ksort, gene_names, top_n = 30) {
  LL <- t(t(LL) / apply(LL, 2, max)) # Scale so maximum loading for each k is 1
  LL <- LL[, ksort]
  colnames(LL) <- names(ksort)
  LL_tib <- as_tibble(LL) |>
    mutate(Gene = gene_names) |>
    pivot_longer(-Gene, names_to = "Component", values_to = "Loading") |>
    mutate(Component = factor(Component, levels = names(ksort))) |>
    group_by(Gene) |>
    mutate(Uniqueness = log10(Loading) - log10(sort(Loading, decreasing = TRUE)[2]))

  shared_genes <- LL_tib |>
    left_join(LL_tib |> select(Gene, OtherLoading = Loading, SharedBy = Component),
              by = join_by(Gene == Gene, Loading < OtherLoading),
              relationship = "many-to-many") |>
    filter(SharedBy != Component) |>
    arrange(Component) |>
    group_by(Component, Gene) |>
    summarize(SharedBy = paste(SharedBy, collapse = ", "))

  top_genes <- LL_tib |>
    group_by(Component) |>
    slice_max(Loading, n = top_n) |>
    left_join(shared_genes, by = c("Component", "Gene")) |>
    ungroup()

  react_tbl <- reactable(
    top_genes,
    columns = list(
      Loading = colDef(name = "Loading", align = "right", cell = function(value) {
        width <- paste0(value * 100, "%")
        bar_chart(width = width, background = "#e1e1e1")
      }),
      Uniqueness = colDef(name = "Uniqueness", align = "center", cell = function(value) {
        label <- round(value, 2)
        bar_chart_pos_neg(label, value)
      }),
      SharedBy = colDef(name = "Larger In")
    ),
    groupBy = "Component"
  )

  return(react_tbl)
}
