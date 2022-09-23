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

  if (!is.null(label)) {
    p <- p + ggrepel::geom_text_repel(label=label, max.overlaps = 20)
  }
  return(p)
}

#' Box-violin plot comparing values between groups
#'
#' @param data A data frame (or a tibble).
#' @param x The grouping variable from the \code{data}.
#' @param y The value variable from the \code{data}.
#' @param color The color variable from the \code{data}.
#' @param palette The color palette for different groups.
#' @param test Perform "wilcox.test" or "t.test" or not test.
#' @param step.increase numeric vector with the increase in fraction of total height for every additional comparison to minimize overlap.
#' @param comparisons	A list of length-2 vectors specifying the groups of interest to be compared. For example to compare groups "A" vs "B" and "B" vs "C", the argument is as follow: comparisons = list(c("A", "B"), c("B", "C"))
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
BetweenStatPlot <- function(data, x, y, color, palette = NULL,
                            test = c('wilcox.test', 't.test', 'none'),
                            comparisons = NULL,
                            step.increase=0.3) {
  stat.formula <- as.formula(paste(y, "~", x))

  test <- match.arg(test, choices = c('wilcox.test', 't.test', 'none'))
  if (test != 'none') {
    if (test == 'wilcox.test') {
      stat_dat <- data %>%
        wilcox_test(stat.formula, comparisons = comparisons)
    }
    if (test == 't.test') {
      stat_dat <- data %>%
        t_test(stat.formula, comparisons = comparisons)
    }
    stat_dat <- stat_dat %>%
      adjust_pvalue() %>%
      p_format(.data$p.adj, digits = 2, leading.zero = FALSE,
               trailing.zero = TRUE, add.p = TRUE, accuracy = 2e-16) %>%
      add_xy_position(x = x, dodge=0.8, step.increase=step.increase)
  }

  x.labs <- paste0(unique(data[,x]), "\n(n=", tabulate(as.factor(data[,x])),")")
  x.num <- length(unique(data[,color])) # number of x types
  if (is.null(palette)) palette <- paint_palette("Spring", x.num, 'continuous')

  p <- data %>%
    ggplot(aes_string(x, y, color = color)) +
    geom_violin(width = 0.8) +
    geom_boxplot(width = 0.3, outlier.shape = NA) +
    theme_classic() +
    theme(legend.position = 'none',
          axis.text = element_text(color='black')) +
    scale_color_manual(values = palette) +
    scale_x_discrete(labels = x.labs) +
    labs(x='')

  if (exists('stat_dat')) {
    p <- p + stat_pvalue_manual(data = stat_dat, label = "p.adj", tip.length = 0.01, size = 3)
  }

  return(p)
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

# volcano plot

#' Volcano plot for DEGs
#'
#' @param data DEseq2 result table.
#' @param title Title of the plot.
#'
#' @return ggplot2 object
#' @export
#'
#' @import dplyr
#' @import ggplot2
ggVolcano <- function(data, title=NULL) {
  data <- data %>%
    dplyr::mutate(stat = if_else(padj < 0.05 & log2FoldChange >= 1, "Up",
                          if_else(padj < 0.05 & log2FoldChange <= -1, "Down", "NS", missing = "NS"), missing = "NS"))
  count.dat <- data %>% dplyr::count(stat) %>% dplyr::mutate(label = paste0(stat, ": ", n))

  data %>%
    ggplot(aes(x=log2FoldChange, y = -log10(padj), color=stat)) +
    geom_point() +
    theme_minimal() +
    theme(legend.position = "top",
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(color='black')) +
    scale_color_manual(values = c("#4F99B4","#808080","#CBC28D"), labels = count.dat$label) +
    geom_vline(xintercept = c(-1, 1), lty = "dashed") +
    geom_hline(yintercept = -log10(0.05), lty = "dashed") +
    labs(color='', title = title)
}

# For adjusting no visible binding
## ggPCA
utils::globalVariables(c("PC1", "PC2", "group"))
## ggVolcano
utils::globalVariables(c("padj", "log2FoldChange"))
