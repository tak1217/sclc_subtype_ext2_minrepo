train_centroids <- function(expr_samples_x_genes, labels) {
  labs <- labels[rownames(expr_samples_x_genes)]
  if (any(is.na(labs))) stop("labels missing for some samples")
  df <- as.data.frame(expr_samples_x_genes)
  df$label <- labs
  centroids <- aggregate(. ~ label, data = df, FUN = mean)
  rownames(centroids) <- centroids$label
  centroids$label <- NULL
  as.matrix(centroids)
}

predict_by_centroid <- function(expr_samples_x_genes, centroids) {
  common <- intersect(colnames(expr_samples_x_genes), colnames(centroids))
  X <- as.matrix(expr_samples_x_genes[, common, drop=FALSE])
  C <- as.matrix(centroids[, common, drop=FALSE])

  cor_mat <- matrix(0, nrow=nrow(X), ncol=nrow(C))
  rownames(cor_mat) <- rownames(X)
  colnames(cor_mat) <- rownames(C)

  for (i in seq_len(nrow(X))) {
    for (j in seq_len(nrow(C))) {
      a <- X[i,]; b <- C[j,]
      if (sd(a) == 0 || sd(b) == 0) {
        cor_mat[i,j] <- 0
      } else {
        cor_mat[i,j] <- suppressWarnings(cor(a,b))
        if (is.na(cor_mat[i,j])) cor_mat[i,j] <- 0
      }
    }
  }
  pred <- colnames(cor_mat)[apply(cor_mat, 1, which.max)]
  names(pred) <- rownames(X)
  pred
}
