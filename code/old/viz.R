
tmp <- Matrix::readMM("~/Github/EmbryoTimecourse2018/download/chimera-tal1/raw_counts.mtx")
Hbb <- tmp[13409, ]
libsizes <- Matrix::colSums(tmp)
rm(tmp)

tib <- tibble(
  Hbb = Hbb,
  libsize = libsizes,
  cell.type = sapply(strsplit(rownames(FF), "_"), `[[`, 6),
  tomato = sapply(strsplit(rownames(FF), "_"), `[[`, 4)
)
tib <- tib %>%
  group_by(cell.type, tomato) %>%
  filter(n() > 100) %>%
  summarize(Hbb = mean(Hbb / libsize))
tib <- tib %>%
  pivot_wider(names_from = tomato, values_from = Hbb, names_prefix = "Tom_")
tib <- tib %>%
  filter(!is.na(Tom_TRUE) & !is.na(Tom_FALSE)) %>%
  filter(!(cell.type %in% c("Stripped", "Doublet")))

ggplot(tib, aes(x = Tom_TRUE, y = Tom_FALSE)) +
  geom_point()

tib <- make.heatmap.tib(FF)
tib <- tib %>%
  group_by(Cell.type, Tomato, Factor) %>%
  summarize(Loading = mean(Loading)) %>%
  ungroup() %>%
  mutate(Tomato = as.character(Tomato)) %>%
  pivot_wider(names_from = Tomato, values_from = Loading, names_prefix = "Tom_") %>%
  mutate(Tom_diff = Tom_TRUE - Tom_FALSE)

ebmf.LL <- fl$L.pm[, order(-fl$pve)]
