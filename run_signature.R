source("R/io.R")
source("R/preprocess.R")
source("R/signatures.R")

dir.create("out", showWarnings = FALSE)

expr <- read_expression("data/expression.csv")
expr <- filter_low_expression(expr, min_mean = 0)
expr_z <- row_zscore(expr)

gene_sets <- read_gene_sets("data/gene_sets.tsv")
scores <- score_signatures(expr_z, gene_sets = gene_sets)
pred <- assign_subtype(scores)

write.csv(scores, "out/signature_scores.csv", quote = FALSE)
write.csv(data.frame(sample=names(pred), pred_subtype=unname(pred)),
          "out/signature_pred_subtype.csv", row.names = FALSE, quote = FALSE)

if (file.exists("data/metadata.csv")) {
  meta <- read_metadata("data/metadata.csv")
  merged <- merge(meta, data.frame(sample=names(pred), pred_subtype=unname(pred)), by="sample", all.x=TRUE)
  ct <- table(merged$true_subtype, merged$pred_subtype)
  write.csv(as.data.frame.matrix(ct), "out/signature_confusion_matrix.csv", quote = FALSE)
}

cat("Done: out/signature_scores.csv, out/signature_pred_subtype.csv\n")
