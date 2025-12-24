# GeneID -> HGNC symbol conversion helpers.

read_geneid_map_tsv <- function(path) {
  df <- read.delim(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("gene_id", "hgnc_symbol") %in% colnames(df))) {
    stop("Mapping file must have columns: gene_id, hgnc_symbol")
  }
  df$gene_id <- as.character(df$gene_id)
  df$hgnc_symbol <- as.character(df$hgnc_symbol)
  df <- df[df$gene_id != "" & df$hgnc_symbol != "", , drop = FALSE]
  setNames(df$hgnc_symbol, df$gene_id)
}

split_genenames_field <- function(x) {
  x <- as.character(x)
  if (is.na(x) || x == "") return(character())
  parts <- unlist(strsplit(x, "[,;]"))
  parts <- trimws(parts)
  parts[parts != ""]
}

read_genenames_to_id_map <- function(path,
                                     include_fields = c("Approved symbol",
                                                        "Previous symbols",
                                                        "Alias symbols",
                                                        "RefSeq IDs",
                                                        "NCBI Gene ID",
                                                        "Ensembl gene ID",
                                                        "UniProt ID(supplied by UniProt)"),
                                     status_filter = "Approved",
                                     include_approved_symbol = TRUE) {
  df <- read.delim(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if (!("Approved symbol" %in% colnames(df))) {
    stop("genenames.tsv must contain column: Approved symbol")
  }
  if (!is.null(status_filter) && "Status" %in% colnames(df)) {
    df <- df[df$Status == status_filter, , drop = FALSE]
  }

  rows <- list()
  for (i in seq_len(nrow(df))) {
    symbol <- trimws(as.character(df$`Approved symbol`[i]))
    if (symbol == "") next

    ids <- character()
    if (include_approved_symbol) ids <- c(ids, symbol)
    for (fld in include_fields) {
      if (!fld %in% colnames(df)) next
      ids <- c(ids, split_genenames_field(df[[fld]][i]))
    }
    if (length(ids) == 0) next

    rows[[length(rows) + 1]] <- data.frame(
      source_id = ids,
      hgnc_symbol = rep(symbol, length(ids)),
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)
  out <- out[!duplicated(out$source_id), , drop = FALSE]
  setNames(out$hgnc_symbol, out$source_id)
}

detect_genenames_id_matches <- function(path, ids,
                                        fields = c("Approved symbol",
                                                   "Previous symbols",
                                                   "Alias symbols",
                                                   "RefSeq IDs",
                                                   "NCBI Gene ID",
                                                   "Ensembl gene ID",
                                                   "UniProt ID(supplied by UniProt)"),
                                        status_filter = "Approved") {
  df <- read.delim(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if (!is.null(status_filter) && "Status" %in% colnames(df)) {
    df <- df[df$Status == status_filter, , drop = FALSE]
  }
  ids <- unique(as.character(ids))
  ids <- ids[!is.na(ids) & ids != ""]

  counts <- integer(0)
  for (fld in fields) {
    if (!fld %in% colnames(df)) next
    if (fld == "Approved symbol") {
      vals <- as.character(df[[fld]])
    } else {
      vals <- unlist(lapply(df[[fld]], split_genenames_field))
    }
    vals <- unique(vals[!is.na(vals) & vals != ""])
    counts[fld] <- sum(ids %in% vals)
  }
  counts
}

map_geneids_to_hgnc <- function(expr_genes_x_samples, geneid_to_symbol,
                                drop_missing = TRUE) {
  ids <- rownames(expr_genes_x_samples)
  if (is.null(ids)) stop("expr must have rownames (NCBI GeneID).")

  sym <- unname(geneid_to_symbol[as.character(ids)])
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
