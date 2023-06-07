library(tidyverse)

ebnmf <- readRDS("output/montoro-sel-ebmf-log-k=30.rds")
nnlm <- readRDS("output/montoro-sel-nmf-log-k=30.rds")
tm <- readRDS("output/montoro-sel-topics-k=30.rds")

fl <- ebnmf$fit
FF <- fl$F.pm
FF <- scale(FF, center = FALSE)

cell.type <- sapply(strsplit(rownames(FF), "_"), `[[`, 5)
cell.type <- factor(cell.type, levels = sort(unique(cell.type))[c(
  1, 3, 4, 2, 7, 5, 6, 13, 11, 12, 9, 8, 10
)])

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

ebnmf.FF <- fl$F.pm[cell.idx, ]
nnlm.FF <- (t(nnlm$fit$h))[cell.idx, ]
tm.FF <- tm$fit$F[cell.idx, ]

make.heatmap.tib <- function(FF) {
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

ggplot(heatmap.tib, aes(x = Component, y = -Cell.idx, fill = Loading)) +
  geom_tile() +
  scale_fill_gradient2(high = "black") +
  labs(y = "") +
  scale_y_continuous(breaks = -label_pos,
                     minor_breaks = NULL,
                     labels = levels(cell.type)) +
  theme_minimal() +
  geom_hline(yintercept = -cell_type_breaks, size = 0.1) +
  facet_wrap(~Method, scales = "free", ncol = 1) +
  theme(legend.position = "none",
        strip.text = element_text(size = 12, margin = margin(2, 2, 2, 2), family = "serif"),
        axis.text.y = element_text(size = 8, family = "serif"),
        axis.text.x = element_text(size = 6, family = "serif"),
        axis.title = element_text(size = 8, margin = margin(2, 2, 2, 2), family = "serif"))

ggsave("figs/montoro_heatmap.png", width = 133, height = 150, units = "mm")
