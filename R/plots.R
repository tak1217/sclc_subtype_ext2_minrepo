plot_rank_quality <- function(quality, out_pdf = "out/nmf_rank_quality_plots.pdf") {
  dir.create(dirname(out_pdf), showWarnings = FALSE, recursive = TRUE)
  pdf(out_pdf, width = 7, height = 5)
  op <- par(no.readonly = TRUE)
  on.exit({par(op); dev.off()}, add = TRUE)

  if ("cophenetic" %in% colnames(quality)) {
    plot(quality$rank, quality$cophenetic, type="b",
         xlab="k (rank)", ylab="Cophenetic correlation",
         main="NMF rank estimation: cophenetic")
  } else { plot.new(); text(0.5,0.5,"No 'cophenetic' column") }

  if ("rss" %in% colnames(quality)) {
    plot(quality$rank, quality$rss, type="b",
         xlab="k (rank)", ylab="RSS",
         main="NMF rank estimation: RSS")
  } else { plot.new(); text(0.5,0.5,"No 'rss' column") }
}

consensus_heatmap_pdf <- function(consensus_mat, out_pdf = "out/consensus_heatmap.pdf",
                                  show_names = FALSE) {
  dir.create(dirname(out_pdf), showWarnings = FALSE, recursive = TRUE)

  d <- as.dist(1 - consensus_mat)
  hc <- hclust(d, method = "average")
  ord <- hc$order
  cm_ord <- consensus_mat[ord, ord]

  if (requireNamespace("pheatmap", quietly = TRUE)) {
    pdf(out_pdf, width = 8, height = 8)
    pheatmap::pheatmap(cm_ord,
                       cluster_rows = FALSE, cluster_cols = FALSE,
                       show_rownames = show_names, show_colnames = show_names,
                       border_color = NA)
    dev.off()
  } else {
    pdf(out_pdf, width = 8, height = 8)
    op <- par(mar=c(2,2,2,2))
    image(1:nrow(cm_ord), 1:ncol(cm_ord), t(cm_ord[nrow(cm_ord):1,]),
          axes=FALSE, xlab="", ylab="", main="Consensus matrix (ordered)")
    par(op)
    dev.off()
  }
  list(order = ord, hc = hc)
}
