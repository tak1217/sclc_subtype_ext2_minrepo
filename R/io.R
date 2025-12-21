read_expression <- function(path) {
  expr <- read.csv(path, row.names = 1, check.names = FALSE)
  expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"
  return(expr) # genes x samples
}

read_metadata <- function(path) {
  meta <- read.csv(path, stringsAsFactors = FALSE)
  if (!("sample" %in% colnames(meta))) stop("metadata must contain a 'sample' column")
  return(meta)
}
