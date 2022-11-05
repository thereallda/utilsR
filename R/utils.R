#' PCA plot from counts matrix
#'
#' @param object A count matrix.
#' @param use.pc Which two PCs to be used, default PC1 in x-axis and PC2 in y-axis.
#' @param color Vector indicates the color mapping of samples, default NULL.
#' @param label Vector of sample names or labels, default NULL.
#' @param shape Vector indicates the shape mapping of samples, default NULL.
#' @param vst.norm Whether to perform \code{vst} transformation, default FALSE.
#' @param palette The color palette for different groups.
#' @param repel Whether to use \code{ggrepel} to avoid overlapping text labels or not, default TRUE.
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
ggPCA <- function(object, use.pc=c(1,2),
                  color=NULL, label=NULL, shape=NULL,
                  vst.norm=FALSE, palette=NULL, repel=TRUE) {
  if (vst.norm) {
    counts_norm <- DESeq2::vst(as.matrix(object))
  } else {
    counts_norm <- object
  }

  # perform PCA
  pca <- prcomp(t(counts_norm))
  pc.var <- round(summary(pca)$importance[2,], 3)
  pca_dat <- as.data.frame(pca$x)

  # check if use.pc exceed the range of pcs
  use.pc <- paste0('PC', use.pc)
  if (!all(use.pc %in% colnames(pca_dat))) {
    stop(use.pc, "exceed the range of PCs.")
  }
  # mapping data
  var.ls <- list(color = color,
                 shape = shape
                 )
  var.length <- unlist(lapply(var.ls, length))
  var.ls <- var.ls[var.length == max(var.length)]
  map_df <- as.data.frame(Reduce(cbind, var.ls))
  colnames(map_df) <- names(var.ls)
  # combine with pca_dat if not empty
  if (!any(dim(map_df) == 0)) {
    pca_dat <- cbind(pca_dat, map_df)
  }

  # generate color palette
  if (is.null(palette)) {
    palette <- paintingr::paint_palette("Spring", length(unique(pca_dat$color)), 'continuous')
  }

  # create aes mapping
  map_ls <- list(x = use.pc[1],
                 y = use.pc[2],
                 color = "color",
                 shape = "shape")
  mapping <- do.call(ggplot2::aes_string, map_ls)

  p <- ggplot(pca_dat, mapping) +
    geom_point(size=3) +
    geom_vline(xintercept=0, color='grey80', lty=2) +
    geom_hline(yintercept=0, color='grey80', lty=2) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = 'top',
          axis.text = element_text(color='black')) +
    scale_color_manual(values = palette) +
    labs(x=paste0(use.pc[1], ': ', pc.var[1]*100, '%'),
         y=paste0(use.pc[2], ': ', pc.var[2]*100, '%'))

  # add text label
  if (!is.null(label)) {
    if (repel) {
      p <- p + ggrepel::geom_text_repel(label=label, max.overlaps = 20, color='black')
    } else {
      p <- p + geom_text(label=label, color='black')
    }
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
      add_xy_position(x=x, dodge=0.8, step.increase=step.increase)
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

#' Volcano plot for DEGs
#'
#' @param data Differential analysis result table.
#' @param lfc.col Column name of the log fold-change, default: "log2FoldChange".
#' @param p.col Column name of the adjusted p-value, default: "padj".
#' @param up.lfc.cutoff Log2 fold-change cutoff for up-regulated significant differential expression genes, default: 1.
#' @param down.lfc.cutoff Log2 fold-change cutoff for down-regulated significant differential expression genes, default: -1.
#' @param p.cutoff P-value (adjusted) cutoff for significant differential expression genes, default: 0.05.
#' @param title Title of the plot.
#'
#' @return ggplot2 object
#' @export
#'
#' @import dplyr
#' @import ggplot2
ggVolcano <- function(data,
                      lfc.col = "log2FoldChange",
                      p.col = "padj",
                      up.lfc.cutoff = 1,
                      down.lfc.cutoff = -1,
                      p.cutoff = 0.05,
                      title = NULL) {
  data <- data %>%
    dplyr::mutate(stat = if_else(!!sym(p.col) < p.cutoff & !!sym(lfc.col) >= up.lfc.cutoff, "Up",
                          if_else(!!sym(p.col) < p.cutoff & !!sym(lfc.col) <= down.lfc.cutoff, "Down", "NS", missing = "NS"), missing = "NS"))
  count.dat <- data %>% dplyr::count(stat) %>% dplyr::mutate(label = paste0(stat, ": ", n))
  data$logp <- -log10(data[,p.col])

  data %>%
    ggplot(aes_(x=as.name(lfc.col), y=quote(logp), color=quote(stat))) +
    geom_point() +
    theme_minimal() +
    theme(legend.position = "top",
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(color='black')) +
    scale_color_manual(values = c("#4F99B4","#808080","#CBC28D"), labels = count.dat$label) +
    geom_vline(xintercept = c(down.lfc.cutoff, up.lfc.cutoff), lty = "dashed") +
    geom_hline(yintercept = -log10(p.cutoff), lty = "dashed") +
    labs(x='log2 Fold Change', y='-log10 FDR', color='', title = title)
}

# For adjusting no visible binding
## ggPCA
utils::globalVariables(c("PC1", "PC2", "group"))
## ggVolcano
utils::globalVariables(c("padj", "log2FoldChange"))
