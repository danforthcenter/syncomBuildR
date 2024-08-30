#' Function to run a UMAP, plot and return the data with coordinates
#'
#' @param df Dataframe to ordinate. Generally UMAP works best on principal components,
#' so this is assumed to be output from \link{pcadf}.
#' @param cols columns to reduce dimensions of. Can be specified with names or positions.
#' Defaults to all column names starting with "pc" then numbers for use with \link{pcadf}.
#' @param color column name used to color points in the umap plot.
#' @param returnData Logical, should data be returned?
#' Defaults to TRUE where data and a ggplot are returned.
#' @param ... Additional arguments passed to \code{uwot::umap}. The n_neighbors, n_components,
#' and pca arguments are potentially commonly useful.
#' @keywords pca, umap
#'
#' @import ggplot2
#' @import uwot
#'
#' @examples
#' print(load("/home/jsumner/Desktop/stargate/SINC/sincUtils/syncomBuilder/cal_output.rdata"))
#' asv <- are_c[[1]]
#' zinbCalibrated <- are_c[[2]][are_c[[2]]$model == "ZINB", "asv"]
#' asv <- are_c[[1]][, c("tissue", "plot", "row", "genotype", "biomass", "sd", zinbCalibrated)]
#' pdf <- pcadf(df = asv, cols = NULL, color = c("tissue", "genotype"), returnData = T, ncp = NULL)
#' udf <- umapdf(df = pdf$data, cols = NULL, color = c("tissue", "genotype"), returnData = T, pca = 50)
#'
#' @export

umapdf <- function(df = NULL, cols = NULL, color = NULL, returnData = TRUE, ...) {
  if (is.null(cols)) {
    cols <- which(grepl("^pc[0-9]+", colnames(df)))
  } else if (is.character(cols) && length(cols) == 1) {
    cols <- which(grepl(cols, colnames(df)))
  } else if (!is.numeric(cols)) {
    cols <- which(colnames(df) %in% cols)
  }

  ump <- uwot::umap(df[, cols], ...)
  coords <- as.data.frame(ump)
  colnames(coords) <- c(paste0("umap", seq_len(ncol(coords))))
  umap.df <- cbind(df[, -cols], coords)
  if (length(color) > 1) {
    umap.df[[paste0(color, collapse = ".")]] <- interaction(umap.df[, c(color)])
    color <- paste0(color, collapse = ".")
  }
  if (is.null(color)) {
    umap.df$dummyVariableForColor <- 1
    color <- "dummyVariableForColor"
  }
  p <- ggplot2::ggplot(umap.df, ggplot2::aes(
    x = .data[["umap1"]], y = .data[["umap2"]],
    color = .data[[color]]
  )) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::labs(x = paste0("UMAP 1"), y = paste0("UMAP 2")) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1))) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(size = 14),
      axis.line.y.left = ggplot2::element_line(),
      axis.line.x.bottom = ggplot2::element_line()
    )
  if (color == "dummyVariableForColor") {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  if (returnData) {
    return(list("data" = umap.df, "plot" = p))
  } else {
    return(p)
  }
}
