# -----------------------------------------------------------------------------
# Functional enrichment analysis
# -----------------------------------------------------------------------------
# GO Biological Process enrichment for top-selected proteins, using
# clusterProfiler with org.Hs.eg.db. The full set of measured proteins is
# used as the background universe, following the methodology in the report.
# -----------------------------------------------------------------------------


#' Map protein symbols to Entrez gene IDs
#'
#' @param symbols Character vector of gene symbols.
#' @return A data frame with SYMBOL and ENTREZID columns.
map_to_entrez <- function(symbols) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("Packages 'clusterProfiler' and 'org.Hs.eg.db' are required. ",
         "Install them via BiocManager::install().")
  }
  suppressMessages(
    clusterProfiler::bitr(
      symbols,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db::org.Hs.eg.db
    )
  )
}


#' Run GO Biological Process enrichment on a set of target proteins
#'
#' @param target_symbols Character vector of target gene symbols (e.g. top
#'   model predictors).
#' @param background_symbols Character vector of background universe symbols
#'   (typically all measured proteins).
#' @param p_cutoff Adjusted p-value threshold. Default 0.05.
#' @return An `enrichResult` object from clusterProfiler.
run_go_enrichment <- function(target_symbols,
                              background_symbols,
                              p_cutoff = 0.05) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package 'clusterProfiler' is required.")
  }

  target_entrez <- unique(map_to_entrez(target_symbols)$ENTREZID)
  bg_entrez <- unique(map_to_entrez(background_symbols)$ENTREZID)

  clusterProfiler::enrichGO(
    gene          = target_entrez,
    universe      = bg_entrez,
    OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = p_cutoff,
    readable      = TRUE
  )
}


#' Dot plot of enriched GO terms
#'
#' @param ego An `enrichResult` object.
#' @param show_n Number of top terms to display. Default 15.
#' @param title Plot title.
#' @return A ggplot object.
plot_go_dotplot <- function(ego, show_n = 15, title = "GO Pathway Enrichment") {
  if (!requireNamespace("enrichplot", quietly = TRUE)) {
    stop("Package 'enrichplot' is required.")
  }
  enrichplot::dotplot(ego, showCategory = show_n) +
    ggplot2::ggtitle(title)
}
