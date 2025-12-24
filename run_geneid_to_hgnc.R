# Convert expression rownames from NCBI GeneID to HGNC symbol.
source("R/io.R")
source("R/geneid_map.R")
source("R/id_map.R")

expr <- read_expression("data/expression.csv")

map_path <- "data/geneid_to_hgnc.tsv"
if (file.exists(map_path)) {
  geneid_to_symbol <- read_geneid_map_tsv(map_path)
} else if (requireNamespace("AnnotationDbi", quietly = TRUE) &&
           requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  ids <- rownames(expr)
  geneid_to_symbol <- fetch_hgnc_map_bioconductor(ids, keytype = "ENTREZID")
} else {
  stop(paste0("No mapping file found at ", map_path, " and org.Hs.eg.db is not installed.\n",
              "Provide data/geneid_to_hgnc.tsv or install org.Hs.eg.db."))
}

expr_hgnc <- map_geneids_to_hgnc(expr, geneid_to_symbol, drop_missing = TRUE)
write.csv(expr_hgnc, "data/expression_hgnc.csv", quote = FALSE)

cat("Done: data/expression_hgnc.csv\n")
