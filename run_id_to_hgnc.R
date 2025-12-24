# Convert expression rownames from various ID types to HGNC symbols.
source("R/io.R")
source("R/id_map.R")
source("R/geneid_map.R")

# ---- settings ----
input_path <- "data/expression.csv"
output_path <- "data/expression_hgnc.csv"
genenames_path <- "data/genenames.tsv"

# One of the org.Hs.eg.db keytypes (e.g., "ENTREZID", "ENSEMBL", "REFSEQ", "UNIPROT").
keytype <- "ENTREZID"

drop_missing <- TRUE
# ------------------

expr <- read_expression(input_path)

if (file.exists(genenames_path)) {
  fields <- c("Approved symbol",
              "Previous symbols",
              "Alias symbols",
              "RefSeq IDs",
              "NCBI Gene ID",
              "Ensembl gene ID",
              "UniProt ID(supplied by UniProt)")
  counts <- detect_genenames_id_matches(genenames_path, rownames(expr), fields = fields)
  if (length(counts) > 0) {
    cat("ID match counts by genenames.tsv column:\n")
    for (nm in names(counts)) {
      cat("  ", nm, ": ", counts[[nm]], "\n", sep = "")
    }
    best_field <- names(counts)[which.max(counts)]
    cat("Using best-matching column: ", best_field, "\n", sep = "")
  } else {
    cat("No matching columns found in genenames.tsv.\n")
    best_field <- NULL
  }
  id_to_symbol <- read_genenames_to_id_map(
    genenames_path,
    include_fields = if (is.null(best_field)) character() else best_field,
    include_approved_symbol = identical(best_field, "Approved symbol")
  )
} else {
  id_to_symbol <- fetch_hgnc_map_bioconductor(rownames(expr), keytype = keytype)
}

expr_hgnc <- map_ids_to_hgnc(expr, id_to_symbol, drop_missing = drop_missing)
write.csv(expr_hgnc, output_path, quote = FALSE)

cat("Done: ", output_path, "\n", sep = "")
