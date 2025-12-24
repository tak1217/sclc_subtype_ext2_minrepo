# SCLC subtype minimal repo (ext2-min policy)

This repo is a **minimal, internally consistent** R scaffold focused on:
- Signature scoring (mean z-score) and SCLC-A/N/P/I assignment
- Consensus NMF with rank estimation (**nmfEstimateRank uses `range=`**)
- Optional centroid classifier

## Install
```r
# Before install NMF, you may need the following:
renv::install("BiocManager")
BiocManager::install("Biobase")
renv::install("NMF")

install.packages("NMF")
install.packages("pheatmap")   # optional (prettier heatmap)
```

## Data format
- `data/expression.csv` or `data/expression.tsv`: **genes x samples** (row names = HGNC symbols or other IDs; see ID mapping)
- `data/metadata.csv`: optional; columns `sample`, `true_subtype` (for validation only)

## Run (RStudio)
```r
source("run_signature.R")
source("run_consensus_nmf.R")
source("run_centroid.R")  # optional
```

## ID mapping (to HGNC symbols)
If your expression row names are not HGNC symbols, convert them first:
```r
source("run_id_to_hgnc.R")
```
This reads `data/genenames.tsv`, detects which column matches your IDs best,
and maps only that column to `Approved symbol`. Output is:
- `data/expression_hgnc.csv`

Then use `data/expression_hgnc.csv` as your expression input (update the run scripts
or overwrite `data/expression.csv`).

## What to edit for paper integration
- Replace gene sets in `data/gene_sets.tsv` (or a `.gmt`) with the paper-derived lists/modules.
- Optionally tune thresholds in `assign_subtype()` (especially for SCLC-I).

## Design choices
- No GSVA/ssGSEA to avoid version/argument mismatches and keep dependencies minimal.
- Heatmap annotation is disabled by default (pheatmap annotation can break if rownames don't align).


## Gene sets configuration
Gene sets are **not hard-coded**. Provide them in an external file:

- Default: `data/gene_sets.tsv` with columns `set` and `gene`
- Also supported: `.gmt`

Example TSV (already included as a template):
```
set	gene
A	ASCL1
...
I	CXCL9
```

To use your paper-derived sets, overwrite `data/gene_sets.tsv` (or provide a GMT and change the path in the run scripts).
