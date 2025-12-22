source("R/io.R")
source("R/preprocess.R")
source("R/signatures.R")
source("R/centroid.R")

dir.create("out", showWarnings = FALSE)

expr <- read_expression("data/expression.csv")
expr <- filter_low_expression(expr, min_mean = 0)
expr_z <- row_zscore(expr)
Z <- samples_x_genes(expr_z) # samples x genes

scores <- score_signatures(expr_z = expr_z, gene_sets = 'data/gene_sets.tsv') # 2025.12.22 add gene_sets
train_labels <- assign_subtype(scores)

centroids <- train_centroids(Z, train_labels)
pred <- predict_by_centroid(Z, centroids)

write.csv(centroids, "out/centroids.csv", quote = FALSE)
write.csv(data.frame(sample=names(pred), centroid_pred_subtype=unname(pred)),
          "out/centroid_pred_subtype.csv", row.names = FALSE, quote = FALSE)

if (file.exists("data/metadata.csv")) {
  meta <- read_metadata("data/metadata.csv")
  merged <- merge(meta, data.frame(sample=names(pred), centroid_pred_subtype=unname(pred)),
                  by="sample", all.x=TRUE)
  ct <- table(merged$true_subtype, merged$centroid_pred_subtype)
  write.csv(as.data.frame.matrix(ct), "out/centroid_confusion_matrix.csv", quote = FALSE)
}

cat("Done: out/centroids.csv, out/centroid_pred_subtype.csv\n")
