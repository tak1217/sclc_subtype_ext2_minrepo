# Generic ID -> HGNC symbol conversion helpers.

read_id_map_tsv <- function(path) {
  df <- read.delim(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("source_id", "hgnc_symbol") %in% colnames(df))) {
    stop("Mapping file must have columns: source_id, hgnc_symbol")
  }
  df$source_id <- as.character(df$source_id)
  df$hgnc_symbol <- as.character(df$hgnc_symbol)
  df <- df[df$source_id != "" & df$hgnc_symbol != "", , drop = FALSE]
  setNames(df$hgnc_symbol, df$source_id)
}

map_ids_to_hgnc <- function(expr_genes_x_samples, id_to_symbol,
                            drop_missing = TRUE) {
  ids <- rownames(expr_genes_x_samples)
  if (is.null(ids)) stop("expr must have rownames (source IDs).")

  sym <- unname(id_to_symbol[as.character(ids)])
  if (drop_missing) {
    keep <- !is.na(sym) & sym != ""
    expr <- expr_genes_x_samples[keep, , drop = FALSE]
    sym <- sym[keep]
  } else {
    expr <- expr_genes_x_samples
    sym[is.na(sym) | sym == ""] <- ids[is.na(sym) | sym == ""]
  }

  sums <- rowsum(expr, group = sym, reorder = FALSE)
  counts <- as.integer(table(sym)[rownames(sums)])
  sweep(sums, 1, counts, "/")
}

fetch_hgnc_map_bioconductor <- function(ids, keytype) {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("org.Hs.eg.db is required for automatic mapping.")
  }
  AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = as.character(ids),
    column = "SYMBOL",
    keytype = keytype,
    multiVals = "first"
  )
}
