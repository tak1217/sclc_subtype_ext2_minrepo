select_variable_genes <- function(expr, n = 300) {
  v <- apply(expr, 1, var, na.rm = TRUE)
  names(sort(v, decreasing = TRUE))[seq_len(min(n, length(v)))]
}

make_nonnegative <- function(X_samples_x_genes, method = c("minmax","shift")) {
  method <- match.arg(method)
  X <- as.matrix(X_samples_x_genes)
  storage.mode(X) <- "numeric"
  if (method == "minmax") {
    mins <- apply(X, 2, min, na.rm = TRUE)
    maxs <- apply(X, 2, max, na.rm = TRUE)
    rng <- maxs - mins
    rng[rng == 0 | is.na(rng)] <- 1
    Xnn <- sweep(X, 2, mins, "-")
    Xnn <- sweep(Xnn, 2, rng, "/")
    Xnn[is.na(Xnn)] <- 0
    return(Xnn)
  } else {
    m <- min(X, na.rm = TRUE)
    Xnn <- X - m
    Xnn[Xnn < 0] <- 0
    Xnn[is.na(Xnn)] <- 0
    return(Xnn)
  }
}

nmf_require <- function() {
  if (!requireNamespace("NMF", quietly = TRUE)) {
    stop("Package 'NMF' is required. Install with: install.packages('NMF')")
  }
}

# rank estimation (your NMF requires 'range=')
estimate_rank <- function(X_nn, ranks = 2:6, nrun = 30, seed = 0, verbose = TRUE) {
  nmf_require()
  set.seed(seed)

  old_cores <- NMF::nmf.getOption("cores")
  NMF::nmf.options(cores = 1)
  on.exit(NMF::nmf.options(cores = old_cores), add = TRUE)

  est <- NMF::nmfEstimateRank(X_nn, range = ranks, method = "brunet",
                             nrun = nrun, verbose = verbose)
  q <- NMF::summary(est)
  return(list(estimate = est, quality = q))
}

run_nmf_consensus <- function(X_nn, k = 4, nrun = 50, seed = 0) {
  nmf_require()
  set.seed(seed)

  old_cores <- NMF::nmf.getOption("cores")
  NMF::nmf.options(cores = 1)
  on.exit(NMF::nmf.options(cores = old_cores), add = TRUE)

  fit <- NMF::nmf(X_nn, rank = k, method = "brunet", nrun = nrun)
  W <- NMF::basis(fit)   # samples x k
  H <- NMF::coef(fit)    # k x genes
  labels <- apply(W, 1, which.max) - 1
  cons <- NMF::consensus(fit)

  list(W = W, H = H, labels = labels, consensus = cons, fit = fit)
}

label_clusters_by_signatures <- function(labels, scores) {
  lab <- labels
  names(lab) <- rownames(scores)

  df <- cbind(scores, nmf_cluster = lab)
  cluster_means <- aggregate(. ~ nmf_cluster, data = df, FUN = mean)
  rownames(cluster_means) <- cluster_means$nmf_cluster
  cluster_means$nmf_cluster <- NULL

  top_sig <- apply(as.matrix(cluster_means), 1, function(x) colnames(cluster_means)[which.max(x)])
  sig_to_sub <- c(A="SCLC-A", N="SCLC-N", P="SCLC-P", I="SCLC-I")
  cluster_to_sub <- sig_to_sub[top_sig]

  labeled <- cluster_to_sub[as.character(lab)]
  names(labeled) <- names(lab)

  list(labeled = labeled,
       cluster_means = cluster_means,
       top_signature = top_sig,
       cluster_to_subtype = cluster_to_sub)
}
