# Signature gene sets are loaded from an external file.
# Default: data/gene_sets.tsv (columns: set, gene)

read_gene_sets_tsv <- function(path) {
  df <- read.delim(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("set","gene") %in% colnames(df))) {
    stop("TSV gene set file must have columns: 'set' and 'gene'")
  }
  df <- df[!is.na(df$set) & !is.na(df$gene) & df$set != "" & df$gene != "", , drop=FALSE]
  split(df$gene, df$set)
}

# Minimal GMT reader (each line: setName<TAB>description<TAB>gene1<TAB>gene2...)
read_gene_sets_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  sets <- list()
  for (ln in lines) {
    if (nchar(ln) == 0) next
    parts <- strsplit(ln, "\t")[[1]]
    if (length(parts) < 3) next
    set_name <- parts[[1]]
    genes <- parts[3:length(parts)]
    genes <- genes[genes != ""]
    sets[[set_name]] <- genes
  }
  if (length(sets) == 0) stop("No gene sets parsed from GMT file")
  sets
}

read_gene_sets <- function(path = "data/gene_sets.tsv") {
  if (!file.exists(path)) {
    stop(paste0("Gene set file not found: ", path, "\n",
                "Create it (e.g., data/gene_sets.tsv) or pass a valid path."))
  }
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv","txt","tab")) {
    return(read_gene_sets_tsv(path))
  } else if (ext %in% c("gmt")) {
    return(read_gene_sets_gmt(path))
  } else if (ext == "") {
    # assume TSV by default
    return(read_gene_sets_tsv(path))
  } else {
    stop(paste0("Unsupported gene set file extension: .", ext, " (supported: tsv/txt/tab/gmt)"))
  }
}

score_signature <- function(expr_z, genes) {
  genes_present <- intersect(genes, rownames(expr_z))
  if (length(genes_present) == 0) {
    return(setNames(rep(0, ncol(expr_z)), colnames(expr_z)))
  }
  colMeans(expr_z[genes_present, , drop = FALSE])
}

# Returns: samples x signatures (data.frame)
score_signatures <- function(expr_z, gene_sets) {
  if (missing(gene_sets) || is.null(gene_sets)) {
    stop("gene_sets is required. Use gene_sets = read_gene_sets('data/gene_sets.tsv').")
  }
  scores <- lapply(gene_sets, function(gs) score_signature(expr_z, gs))
  df <- as.data.frame(scores, check.names = FALSE)
  rownames(df) <- colnames(expr_z)
  return(df)
}

assign_subtype <- function(scores, tf_cols = c("A","N","P"),
                           tf_threshold = 0.10, i_threshold = 0.30) {
  required <- c(tf_cols, "I")
  miss <- setdiff(required, colnames(scores))
  if (length(miss) > 0) stop(paste("Scores missing required columns:", paste(miss, collapse=", ")))

  tf_mat <- as.matrix(scores[, tf_cols, drop = FALSE])
  tf_max <- apply(tf_mat, 1, max)
  tf_argmax <- tf_cols[apply(tf_mat, 1, which.max)]

  is_I <- (tf_max < tf_threshold) & (scores[, "I"] > i_threshold)

  pred <- tf_argmax
  pred[is_I] <- "I"

  mapping <- c(A="SCLC-A", N="SCLC-N", P="SCLC-P", I="SCLC-I")
  pred <- unname(mapping[pred])
  names(pred) <- rownames(scores)
  return(pred)
}
