filter_low_expression <- function(expr, min_mean = 0) {
  keep <- rowMeans(expr, na.rm = TRUE) >= min_mean
  expr[keep, , drop = FALSE]
}

row_zscore <- function(expr) {
  mu <- rowMeans(expr, na.rm = TRUE)
  sdv <- apply(expr, 1, sd, na.rm = TRUE)
  sdv[sdv == 0 | is.na(sdv)] <- 1
  z <- (expr - mu) / sdv
  z[is.na(z)] <- 0
  z
}

samples_x_genes <- function(expr_genes_x_samples) {
  t(expr_genes_x_samples)
}
