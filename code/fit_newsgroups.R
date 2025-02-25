library("tidyverse")
library("flashier")
library("Matrix")

library("reticulate")
source_python("code/read_pickle.py")
newsgroups <- read_pickle_file("data/newsgroups.pkl")

topic_names <- c('alt.atheism',
                 'comp.graphics',
                 'comp.os.ms-windows.misc',
                 'comp.sys.ibm.pc.hardware',
                 'comp.sys.mac.hardware',
                 'comp.windows.x',
                 'misc.forsale',
                 'rec.autos',
                 'rec.motorcycles',
                 'rec.sport.baseball',
                 'rec.sport.hockey',
                 'sci.crypt',
                 'sci.electronics',
                 'sci.med',
                 'sci.space',
                 'soc.religion.christian',
                 'talk.politics.guns',
                 'talk.politics.mideast',
                 'talk.politics.misc',
                 'talk.religion.misc')

library("tidytext")
# Tokenize documents.
newsgroups_tokens <- newsgroups |>
  mutate(document = 1:n()) |>
  unnest_tokens(output = word, input = text)

# Only keep recognizable English words (dictionary includes names).
newsgroups_tokens <- newsgroups_tokens |>
  filter(word %in% qdapDictionaries::GradyAugmented)

# Stem words.
newsgroups_tokens <- newsgroups_tokens |>
  mutate(word = SnowballC::wordStem(word, )) |>
  count(document, word)

# Remove one-letter stems.
newsgroups_tokens <- newsgroups_tokens |>
  filter(str_length(word) > 1)

# Remove words appearing in fewer than ten documents.
rare_tokens <- newsgroups_tokens |>
  group_by(word) |>
  summarize(n_docs = n(), total_count = sum(n)) |>
  filter(n_docs < 10)
newsgroups_tokens <- newsgroups_tokens |>
  anti_join(rare_tokens, by = "word")

# Only keep documents with 50 words or more.
newsgroups_tokens <- newsgroups_tokens |>
  group_by(document) |>
  filter(sum(n) >= 50) |>
  ungroup()

# Create a sparse data matrix.
newsgroups_dtm <- newsgroups_tokens |>
  cast_dtm(document, word, n)
newsgroups_mat <- sparseMatrix(
  i = newsgroups_dtm$i,
  j = newsgroups_dtm$j,
  x = newsgroups_dtm$v,
  dims = c(newsgroups_dtm$nrow, newsgroups_dtm$ncol)
)
rownames(newsgroups_mat) <- rownames(newsgroups_dtm)
colnames(newsgroups_mat) <- colnames(newsgroups_dtm)
topics <- newsgroups$target[as.numeric(rownames(newsgroups_mat))]

# Normalize counts.
size_factors <- rowSums(newsgroups_mat)
size_factors <- size_factors / median(size_factors)

# Transform data.
pc <- 1
dat <- newsgroups_mat / size_factors
dat <- log1p(dat / pc)
saveRDS(dat, "./data/newsgroups_pp.rds")

K <- 30
maxiter <- 100

t0 <- Sys.time()

fl <- flash_init(dat, var_type = 0) |>
  flash_greedy(
    Kmax = min(K, 10),
    ebnm_fn = ebnm_point_exponential
  )

while (fl$n_factors < K) {
  fl <- fl |>
    flash_backfit(maxiter = 10, verbose = 3) |>
    flash_greedy(
      Kmax = min(K - fl$n_factors, 10),
      ebnm_fn = ebnm_point_exponential
    )
}

fl <- fl |>
  flash_backfit(maxiter = maxiter, verbose = 3)

t1 <- Sys.time()

fl <- flash_factors_reorder(fl, order(fl$pve, decreasing = TRUE))

# Select the 10 components that vary most across newsgroups topics.
topic_means <- apply(fl$L_pm, 2, tapply, topic_names[topics + 1], mean)
which_k <- order(-apply(topic_means, 2, max) / apply(topic_means, 2, min))[1:10]

flash_plot_structure(
  fl, pm_which = "loadings", pm_groups = topic_names[topics + 1], gap = 10, kset = which_k
)
F_pm <- fl$F_pm
colnames(F_pm) <- paste0("k", 1:ncol(F_pm))
fastTopics::annotation_heatmap(F_pm[, which_k], feature_sign = "positive", n = 4)
