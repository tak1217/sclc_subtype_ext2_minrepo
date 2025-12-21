source("R/io.R")
source("R/preprocess.R")
source("R/signatures.R")
source("R/nmf_cluster.R")
source("R/plots.R")

dir.create("out", showWarnings = FALSE)

ranks <- 2:6
nrun_rank <- 30
chosen_k <- 4
nrun_consensus <- 50
seed <- 0

expr <- read_expression("data/expression.csv")
expr <- filter_low_expression(expr, min_mean = 0)

expr_z <- row_zscore(expr)
scores <- score_signatures(expr_z)

genes <- select_variable_genes(expr, n = 300)
X <- samples_x_genes(expr[genes, , drop = FALSE])  # samples x genes
X_nn <- make_nonnegative(X, method = "minmax")

est <- estimate_rank(X_nn, ranks = ranks, nrun = nrun_rank, seed = seed, verbose = TRUE)
quality <- est$quality
write.csv(quality, "out/nmf_rank_quality.csv", quote = FALSE, row.names = FALSE)
plot_rank_quality(quality, out_pdf = "out/nmf_rank_quality_plots.pdf")

res <- run_nmf_consensus(X_nn, k = chosen_k, nrun = nrun_consensus, seed = seed)
labels <- res$labels
names(labels) <- rownames(X)

lab_res <- label_clusters_by_signatures(labels, scores)

write.csv(data.frame(sample=names(labels), nmf_cluster=unname(labels)),
          "out/consensus_nmf_cluster.csv", row.names = FALSE, quote = FALSE)
write.csv(data.frame(sample=names(lab_res$labeled), nmf_labeled_subtype=unname(lab_res$labeled)),
          "out/consensus_nmf_labeled_subtype.csv", row.names = FALSE, quote = FALSE)

cons <- res$consensus
write.csv(cons, "out/consensus_matrix.csv", quote = FALSE)
heat_res <- consensus_heatmap_pdf(cons, out_pdf = "out/consensus_heatmap.pdf", show_names = FALSE)
ord <- heat_res$order
write.csv(cons[ord, ord], "out/consensus_matrix_ordered.csv", quote = FALSE)

if (file.exists("data/metadata.csv")) {
  meta <- read_metadata("data/metadata.csv")
  merged <- merge(meta, data.frame(sample=names(lab_res$labeled), nmf_labeled_subtype=unname(lab_res$labeled)),
                  by="sample", all.x=TRUE)
  ct <- table(merged$true_subtype, merged$nmf_labeled_subtype)
  write.csv(as.data.frame.matrix(ct), "out/consensus_nmf_confusion_matrix.csv", quote = FALSE)
}

cat("Done: rank metrics + consensus NMF outputs in out/\n")
