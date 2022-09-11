#' PCA plot from counts matrix
#'
#' @param object A count matrix.
#' @param group Vector of sample groups.
#' @param label Vector of sample names or labels.
#' @param vst.norm if TRUE perform vst transformation.
#' @param palette The color palette for different groups.
#'
#' @return ggplot2 object
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom DESeq2 vst
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats prcomp
#' @importFrom paintingr paint_palette
ggPCA <- function(object, group, label=NULL, vst.norm=FALSE, palette=NULL) {
  if (vst.norm) {
    counts_norm <- DESeq2::vst(as.matrix(object))
  } else {
    counts_norm <- object
  }

  pca <- prcomp(t(counts_norm))
  pc.var <- round(summary(pca)$importance[2,], 3)
  pca_dat <- as.data.frame(pca$x) %>%
    mutate(group = group)

  if (is.null(palette)) {
    palette <- paintingr::paint_palette("Spring", length(unique(pca_dat$group)), 'continuous')
  }

  p <- pca_dat %>%
    ggplot(aes(x=PC1, y=PC2)) +
    geom_point(aes(color=group), size=3) +
    geom_vline(xintercept=0, color='grey80', lty=2) +
    geom_hline(yintercept=0, color='grey80', lty=2) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = 'top',
          axis.text = element_text(color='black')) +
    scale_color_manual(values = palette) +
    labs(x=paste0('PC1: ', pc.var[1]*100, '%'),
         y=paste0('PC2: ', pc.var[2]*100, '%'))

  if (is.null(label)) {
    label <- colnames(object)
  }
  p <- p + ggrepel::geom_text_repel(label=label, max.overlaps = 20)
  return(p)
}

#' Box-violin plot comparing values between groups
#'
#' @param data A data frame (or a tibble).
#' @param x The grouping variable from the \code{data}.
#' @param y The value variable from the \code{data}.
#' @param color The color variable from the \code{data}.
#' @param palette The color palette for different groups.
#'
#' @return ggplot2 object
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @import rstatix
#' @importFrom paintingr paint_palette
#' @importFrom ggpubr stat_pvalue_manual
#' @importFrom stats as.formula
#' @importFrom rlang .data
BetweenStatPlot <- function(data, x, y, color, palette = NULL) {
  stat.formula <- as.formula(paste(y, "~", x))
  stat_dat <- data %>%
    rstatix::wilcox_test(stat.formula) %>%
    rstatix::adjust_pvalue() %>%
    rstatix::p_format(.data$p.adj, digits = 2, leading.zero = FALSE,
             trailing.zero = TRUE, add.p = TRUE, accuracy = 2e-16) %>%
    rstatix::add_xy_position(x = x, dodge=0.8, step.increase=0.2)

  x.labs <- paste0(unique(data[,x]), "\n(n=", tabulate(data[,x]),")")
  x.num <- length(unique(data[,x])) # number of x types
  if (is.null(palette)) palette <- paint_palette("Spring", x.num, 'continuous')

  data %>%
    ggplot(aes_string(x, y, color = color)) +
    geom_violin(width = 0.8) +
    geom_boxplot(width = 0.3, outlier.shape = NA) +
    theme_classic() +
    theme(legend.position = 'none',
          axis.text = element_text(color='black')) +
    scale_color_manual(values = palette) +
    scale_x_discrete(labels = x.labs) +
    ggpubr::stat_pvalue_manual(data = stat_dat, label = "p.adj", tip.length = 0.01, size = 3)+
    labs(x='')
}

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
#' @importFrom AnnotationDbi Term
#' @importFrom gprofiler2 gconvert
gene2goterm <- function(query, organism = 'hsapiens', ...) {
  goterm <- Term(GOTERM)
  goterm_df <- data.frame(target=names(goterm), term_name=goterm)
  geneGO <- gconvert(query, organism = organism, target = 'GO', ...)
  geneGO <- merge(geneGO, goterm_df, by='target')
  geneGO <- geneGO[,c('input_number', 'input', 'target_number', 'target', 'term_name', 'name', 'description', 'namespace')]
  return(geneGO)
}

# For adjusting no visible binding
## ggPCA
utils::globalVariables(c("PC1", "PC2", "group"))
