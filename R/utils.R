#' Convert gene id to GO term
#'
#' @param query character vector that can consist of mixed types of gene IDs (proteins, transcripts, microarray IDs, etc), SNP IDs, chromosomal intervals or term IDs.
#' @param organism organism name. Organism names are constructed by concatenating the first letter of the name and the family name. Example: human - 'hsapiens', mouse - 'mmusculus'.
#' @param ... Additional parameters that can be passed to \code{gconvert}
#'
#' @return The input contains multiple fields including the query id, GO term id and term name, converted gene name, description and namespace.
#' @export
#'
#' @examples
#' gene2goterm(c("ENSMUSG00000025981", "ENSMUSG00000057363"),
#' organism = 'mmusculus')
#'
#' @importFrom GO.db GOTERM
#' @importFrom AnnotationDbi Term Ontology
#' @importFrom gprofiler2 gconvert
gene2goterm <- function(query, organism = 'hsapiens', ...) {
  goterm_df <- data.frame(target=names(AnnotationDbi::Term(GO.db::GOTERM)),
                       term=AnnotationDbi::Term(GO.db::GOTERM),
                       ont=AnnotationDbi::Ontology(GO.db::GOTERM))
  geneGO <- gprofiler2::gconvert(query, organism = organism, target = 'GO', ...)
  geneGO <- merge(geneGO, goterm_df, by='target')
  geneGO <- geneGO[,c('input', 'name', 'description', 'target', 'term', 'ont')]
  return(geneGO)
}

#' Volcano plot for DEGs
#'
#' @param data Differential analysis result table.
#' @param lfc.col Column name of the log fold-change, default: "log2FoldChange".
#' @param p.col Column name of the adjusted p-value, default: "padj".
#' @param up.lfc.cutoff Log2 fold-change cutoff for up-regulated significant differential expression genes, default: 1.
#' @param down.lfc.cutoff Log2 fold-change cutoff for down-regulated significant differential expression genes, default: -1.
#' @param p.cutoff P-value (adjusted) cutoff for significant differential expression genes, default: 0.05.
#' @param title Title of the plot.
#' @param label.genes Vector of gene ids to label. Should be matched with the rowname of \code{data}
#' @param repel Whether to use \code{ggrepel} to avoid overlapping text labels or not, default TRUE.
#' @param pt.color Vector of colors for point
#'
#' @return ggplot2 object
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
ggVolcano <- function(data,
                      lfc.col = "log2FoldChange",
                      p.col = "padj",
                      up.lfc.cutoff = 1,
                      down.lfc.cutoff = -1,
                      p.cutoff = 0.05,
                      title = NULL,
                      label.genes = NULL,
                      repel = TRUE,
                      pt.color = c("#4F99B4","#808080","#CBC28D")) {
  data <- data %>%
    dplyr::mutate(stat = if_else(!!sym(p.col) < p.cutoff & !!sym(lfc.col) >= up.lfc.cutoff, "Up",
                                 if_else(!!sym(p.col) < p.cutoff & !!sym(lfc.col) <= down.lfc.cutoff, "Down", "NS", missing = "NS"), missing = "NS"))
  count.dat <- data %>% dplyr::count(stat) %>% dplyr::mutate(label = paste0(stat, ": ", n))
  data$logp <- -log10(data[,p.col])

  if (length(pt.color) != 3) {
    stop("Provided colors for point should be 3.")
  }

  p1 <- ggplot(data, aes_(x=as.name(lfc.col), y=quote(logp))) +
    geom_point(aes_(color=quote(stat))) +
    theme_classic() +
    theme(legend.position = "top",
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(color='black')) +
    scale_color_manual(values = pt.color, labels = count.dat$label) +
    geom_vline(xintercept = c(down.lfc.cutoff, up.lfc.cutoff), lty = "dashed") +
    geom_hline(yintercept = -log10(p.cutoff), lty = "dashed") +
    labs(x='log2 Fold Change', y='-log10 FDR', color='', title = title)
  # create gene column
  data$gene <- rownames(data)
  if (!is.null(label.genes)) {
    if (!all(unique(label.genes) %in% rownames(data))) {
      stop("All provided labels should be included in rownames of DE table.")
    } else {
      if (repel) {
        p1 <- p1 + ggrepel::geom_text_repel(aes(label=gene),
                                            max.overlaps = 20,
                                            color="black",
                                            min.segment.length = 0,
                                            data=subset(data, gene %in% label.genes))
      } else {
        p1 <- p1 + geom_text(aes(label=gene),
                             color="black",
                             data=subset(data, gene %in% label.genes))
      }
    }
  }

  return(p1)
}

# For adjusting no visible binding
## ggPCA
utils::globalVariables(c("PC1", "PC2", "group"))
## ggVolcano
utils::globalVariables(c("padj", "log2FoldChange","gene"))
