# Load helpers for I/O, preprocessing, signatures, NMF, and plotting.
source("R/io.R")
source("R/preprocess.R")
source("R/signatures.R")
source("R/nmf_cluster.R")
source("R/plots.R")

dir.create("out", showWarnings = FALSE)

# NMF rank search and consensus settings.
# ranks: candidate k values for rank estimation.
# nrun_rank: repeated NMF runs per k for stability metrics.
# chosen_k: final k used for consensus NMF.
# nrun_consensus: repeated runs to build the consensus matrix.
# seed: reproducible random initialization for NMF.
ranks <- 2:6
nrun_rank <- 30
chosen_k <- 4
nrun_consensus <- 50
seed <- 0

# Read expression data (genes x samples) and basic filtering.
# min_mean=0 keeps all genes; adjust to drop lowly expressed genes.
expr <- read_expression("data/expression.csv")
expr <- filter_low_expression(expr, min_mean = 0)

# Signature scoring uses row-wise z-scores per gene.
# expr_z is genes x samples with per-gene mean=0, sd=1.
expr_z <- row_zscore(expr)
scores <- score_signatures(expr_z)

# Select variable genes and prepare nonnegative matrix for NMF.
# Variable genes are chosen by variance; NMF uses samples x genes.
genes <- select_variable_genes(expr, n = 300)
X <- samples_x_genes(expr[genes, , drop = FALSE])  # samples x genes
# Nonnegative transform required by NMF; "minmax" scales each gene to [0,1].
X_nn <- make_nonnegative(X, method = "minmax")

# Evaluate candidate ranks; save metrics and plots.
# quality includes cophenetic/RSS if available from NMF::summary().
est <- estimate_rank(X_nn, ranks = ranks, nrun = nrun_rank, seed = seed, verbose = TRUE)
quality <- est$quality
write.csv(quality, "out/nmf_rank_quality.csv", quote = FALSE, row.names = FALSE)
plot_rank_quality(quality, out_pdf = "out/nmf_rank_quality_plots.pdf")

# Run consensus NMF at chosen_k and label clusters via signature scores.
# labels are 0-based cluster IDs (NMF::basis max component per sample).
res <- run_nmf_consensus(X_nn, k = chosen_k, nrun = nrun_consensus, seed = seed)
labels <- res$labels
names(labels) <- rownames(X)

# Assign each NMF cluster a subtype label by mean signature score.
lab_res <- label_clusters_by_signatures(labels, scores)

# Export per-sample cluster labels and labeled subtypes.
write.csv(data.frame(sample=names(labels), nmf_cluster=unname(labels)),
          "out/consensus_nmf_cluster.csv", row.names = FALSE, quote = FALSE)
write.csv(data.frame(sample=names(lab_res$labeled), nmf_labeled_subtype=unname(lab_res$labeled)),
          "out/consensus_nmf_labeled_subtype.csv", row.names = FALSE, quote = FALSE)

# Save consensus matrix and a heatmap; also save ordered matrix by clustering.
# heat_res$order is the dendrogram order used for the heatmap.
cons <- res$consensus
write.csv(cons, "out/consensus_matrix.csv", quote = FALSE)
heat_res <- consensus_heatmap_pdf(cons, out_pdf = "out/consensus_heatmap.pdf", show_names = FALSE)
ord <- heat_res$order
write.csv(cons[ord, ord], "out/consensus_matrix_ordered.csv", quote = FALSE)

# Optional confusion matrix if ground-truth metadata is provided.
# Requires data/metadata.csv with columns: sample, true_subtype.
if (file.exists("data/metadata.csv")) {
  meta <- read_metadata("data/metadata.csv")
  merged <- merge(meta, data.frame(sample=names(lab_res$labeled), nmf_labeled_subtype=unname(lab_res$labeled)),
                  by="sample", all.x=TRUE)
  ct <- table(merged$true_subtype, merged$nmf_labeled_subtype)
  write.csv(as.data.frame.matrix(ct), "out/consensus_nmf_confusion_matrix.csv", quote = FALSE)
}

cat("Done: rank metrics + consensus NMF outputs in out/\n")
