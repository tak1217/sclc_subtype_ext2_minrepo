# Load helpers for I/O, preprocessing, and signature scoring.
source("R/io.R")
source("R/preprocess.R")
source("R/signatures.R")

dir.create("out", showWarnings = FALSE)

# Read expression data (genes x samples) and basic filtering.
expr <- read_expression("data/expression.csv")
expr <- filter_low_expression(expr, min_mean = 0)
# Per-gene z-score for signature scoring.
expr_z <- row_zscore(expr)

# Load external gene sets and compute signature scores.
gene_sets <- read_gene_sets("data/gene_sets.tsv")
scores <- score_signatures(expr_z, gene_sets = gene_sets)
# Assign SCLC-A/N/P/I based on signature scores.
pred <- assign_subtype(scores)

# Export scores and predicted subtype labels.
write.csv(scores, "out/signature_scores.csv", quote = FALSE)
write.csv(data.frame(sample=names(pred), pred_subtype=unname(pred)),
          "out/signature_pred_subtype.csv", row.names = FALSE, quote = FALSE)

# Optional confusion matrix if ground-truth metadata is provided.
if (file.exists("data/metadata.csv")) {
  meta <- read_metadata("data/metadata.csv")
  merged <- merge(meta, data.frame(sample=names(pred), pred_subtype=unname(pred)), by="sample", all.x=TRUE)
  ct <- table(merged$true_subtype, merged$pred_subtype)
  write.csv(as.data.frame.matrix(ct), "out/signature_confusion_matrix.csv", quote = FALSE)
}

cat("Done: out/signature_scores.csv, out/signature_pred_subtype.csv\n")
