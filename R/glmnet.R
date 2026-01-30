# glmnet marker selection (multinomial), following the description you quoted
# - Fit multinomial GLM with penalized likelihood (glmnet)
# - Repeat 20 times; each time do 10-fold CV to pick best alpha and features
# - Keep features that appear in >=80% of repeats
# - Refit on those features with 10-fold CV to pick final alpha and features

suppressPackageStartupMessages({
  library(glmnet)
})

run_glmnet_markers <- function(expr, y,
                              n_repeats = 20,
                              n_folds   = 10,
                              alpha_grid = c(0, 0.25, 0.5, 0.75, 1),
                              lambda_choice = c("lambda.1se", "lambda.min"),
                              feature_keep_frac = 0.80,
                              seed = 1) {
  lambda_choice <- match.arg(lambda_choice)

  # expr: genes x samples OR samples x genes; we convert to samples x genes for glmnet
  if (ncol(expr) == length(y)) {
    # assume genes x samples
    x <- t(expr)
  } else if (nrow(expr) == length(y)) {
    # assume samples x genes
    x <- expr
  } else {
    stop("Dimensions don't match: y length must equal nrow(expr) or ncol(expr).")
  }

  y <- as.factor(y)
  if (nlevels(y) < 2) stop("y must have at least 2 classes.")

  # optional: ensure numeric matrix
  x <- as.matrix(x)
  storage.mode(x) <- "double"

  # ---------- repeat CV to select alpha + features ----------
  set.seed(seed)
  selected_by_repeat <- vector("list", n_repeats)
  best_alpha_each    <- numeric(n_repeats)

  for (r in seq_len(n_repeats)) {
    # new CV splits each repeat
    foldid <- sample(rep(seq_len(n_folds), length.out = nrow(x)))

    best_cvm   <- Inf
    best_alpha <- NA_real_
    best_fit   <- NULL

    for (a in alpha_grid) {
      fit <- cv.glmnet(
        x = x, y = y,
        family = "multinomial",
        alpha = a,
        foldid = foldid,
        type.measure = "class",   # misclassification error (works well for multinomial)
        standardize = TRUE
      )

      # choose best alpha by minimum CV error
      min_cvm <- min(fit$cvm, na.rm = TRUE)
      if (min_cvm < best_cvm) {
        best_cvm   <- min_cvm
        best_alpha <- a
        best_fit   <- fit
      }
    }

    best_alpha_each[r] <- best_alpha

    # extract non-zero features at chosen lambda
    lam <- if (lambda_choice == "lambda.1se") best_fit$lambda.1se else best_fit$lambda.min
    beta_list <- coef(best_fit, s = lam)  # list of class-specific coefficient matrices

    # union of non-zero genes across classes (excluding intercept)
    genes <- colnames(x)
    nz <- logical(length(genes)); names(nz) <- genes

    for (k in names(beta_list)) {
      b <- beta_list[[k]]                 # dgCMatrix: (Intercept + genes) x 1
      idx <- which(b[-1, 1] != 0)         # drop intercept
      if (length(idx)) nz[idx] <- TRUE
    }
    selected_by_repeat[[r]] <- names(nz)[nz]
  }

  # ---------- keep features that appear in >=80% repeats ----------
  all_sel <- unlist(selected_by_repeat)
  freq <- sort(table(all_sel), decreasing = TRUE)
  keep_n <- ceiling(feature_keep_frac * n_repeats)
  stable_features <- names(freq)[freq >= keep_n]

  # ---------- refit using stable features to pick final alpha + features ----------
  if (length(stable_features) == 0) {
    warning("No stable features found at the given threshold; returning repeat results only.")
    return(list(
      stable_features = character(0),
      best_alpha_each = best_alpha_each,
      selected_by_repeat = selected_by_repeat,
      final = NULL
    ))
  }

  x2 <- x[, stable_features, drop = FALSE]
  set.seed(seed + 999)
  foldid2 <- sample(rep(seq_len(n_folds), length.out = nrow(x2)))

  best_cvm2   <- Inf
  best_alpha2 <- NA_real_
  best_fit2   <- NULL

  for (a in alpha_grid) {
    fit2 <- cv.glmnet(
      x = x2, y = y,
      family = "multinomial",
      alpha = a,
      foldid = foldid2,
      type.measure = "class",
      standardize = TRUE
    )
    min_cvm2 <- min(fit2$cvm, na.rm = TRUE)
    if (min_cvm2 < best_cvm2) {
      best_cvm2   <- min_cvm2
      best_alpha2 <- a
      best_fit2   <- fit2
    }
  }

  lam2 <- if (lambda_choice == "lambda.1se") best_fit2$lambda.1se else best_fit2$lambda.min
  beta_list2 <- coef(best_fit2, s = lam2)

  genes2 <- colnames(x2)
  nz2 <- logical(length(genes2)); names(nz2) <- genes2
  for (k in names(beta_list2)) {
    b <- beta_list2[[k]]
    idx <- which(b[-1, 1] != 0)
    if (length(idx)) nz2[idx] <- TRUE
  }
  final_features <- names(nz2)[nz2]

  list(
    stable_features = stable_features,
    best_alpha_each = best_alpha_each,
    selected_by_repeat = selected_by_repeat,
    final = list(
      alpha = best_alpha2,
      lambda = lam2,
      features = final_features,
      cv_fit = best_fit2
    )
  )
}

# ------------------
# Example usage:
# expr: genes x samples (recommended) or samples x genes
# y: subtype labels (factor/character), length = number of samples

# res <- run_glmnet_markers(expr, y,
#                           n_repeats=20, n_folds=10,
#                           alpha_grid=c(0,0.25,0.5,0.75,1),
#                           lambda_choice="lambda.1se",
#                           feature_keep_frac=0.80,
#                           seed=1)
# res$final$alpha
# length(res$final$features)
# head(res$final$features)