#' Function to run a PCA, plot and return the data with PC coordinates
#'
#' @param df Dataframe to ordinate
#' @param cols columns to reduce dimensions of. Can be specified with names or positions. Defaults to
#' all column names containing "ASV".
#' @param color column name used to color points in the pca plot.
#' @param returnData Logical, should data be returned? Defaults to TRUE where data and a ggplot are
#' eturned.
#' @param ncp Optional, number of principal components to return attached to dataframe if data is
#' returned. Defaults to all.
#' @param umap Logical, should a UMAP also be performed? Defaults to FALSE. If TRUE then UAMP1 and 2
#' columns will be added and a list of 2 ggplots will be returned.
#' @param distance Distance to use. The default ("euclidean") will use \code{FactoMineR::PCA}, distances
#' supported by \code{vegan::vegdist} will use \code{ape::pcoa}.
#' @keywords pca umap
#' @return a ggplot, optionally a dataframe and FactoMineR::PCA output if returnData is TRUE.
#'
#' @import ggplot2
#' @importFrom FactoMineR PCA
#' @importFrom uwot umap
#' @importFrom ape pcoa
#'
#' @examples
#'
#' x <- pcadf(
#'   df = asv, cols = NULL,
#'   color = "tissue", returnData = TRUE, ncp = NULL
#' )
#' x$plot
#'
#' @export

pcadf <- function(df = NULL, cols = NULL, color = NULL, returnData = TRUE, ncp = NULL, umap = FALSE,
                  distance = "euclidean") {
  if (is.null(cols)) {
    cols <- which(grepl("ASV", colnames(df)))
  } else if (is.character(cols) && length(cols) == 1) {
    cols <- which(grepl(cols, colnames(df)))
  } else if (!is.numeric(cols)) {
    cols <- which(colnames(df) %in% cols)
  }
  if (is.null(ncp)) {
    ncp <- min(dim(df[, cols])) - 1
  }
  if (distance == "euclidean") {
    pca <- FactoMineR::PCA(df[, cols], ncp = ncp, graph = FALSE)
    pc1Var <- round(pca$eig[1, 2], 3)
    pc2Var <- round(pca$eig[2, 2], 3)
    coords <- as.data.frame(pca$ind)
    coords <- coords[, grepl("coord", colnames(coords))]
    colnames(coords) <- gsub("coord.Dim.", "pc", colnames(coords))
    pca.df <- cbind(df[, -cols], coords)
  } else {
    d <- vegan::vegdist(as.matrix(df[, cols]), method = distance)
    pca <- ape::pcoa(d, correction = "none")
    pc1Var <- round(pca$values$Relative_eig[1], 3) * 100
    pc2Var <- round(pca$values$Relative_eig[2], 3) * 100
    coords <- as.data.frame(pca$vectors[, seq_len(min(ncp, ncol(pca$vectors)))])
    colnames(coords) <- gsub("Axis.", "pc", colnames(coords))
    pca.df <- cbind(df[, -cols], coords)
  }
  if (length(color) > 1) {
    pca.df[[paste0(color, collapse = ".")]] <- interaction(pca.df[, c(color)])
    color <- paste0(color, collapse = ".")
  }
  if (is.null(color)) {
    pca.df$dummyVariableForColor <- 1
    color <- "dummyVariableForColor"
  }
  p <- ggplot2::ggplot(pca.df, ggplot2::aes(
    x = .data[["pc1"]],
    y = .data[["pc2"]],
    color = .data[[color]]
  )) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::labs(x = paste0("PC 1 (", pc1Var, "%)"), y = paste0("PC 2 (", pc2Var, "%)")) +
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
  if (umap) {
    pca.umap <- uwot::umap(pca.df[, grepl("pc", colnames(pca.df))])
    umap.df <- cbind(pca.df, as.data.frame(pca.umap))
    colnames(umap.df) <- c(colnames(pca.df), "UMAP1", "UMAP2")
    pca.df <- umap.df

    p2 <- ggplot2::ggplot(pca.df, ggplot2::aes(
      x = .data[["UMAP1"]],
      y = .data[["UMAP2"]], color = .data[[color]]
    )) +
      ggplot2::geom_point(alpha = 0.85) +
      ggplot2::labs(x = "UMAP1", y = "UMAP2") +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1))) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(size = 14),
        axis.line.y.left = ggplot2::element_line(),
        axis.line.x.bottom = ggplot2::element_line()
      )
    if (color == "dummyVariableForColor") {
      p2 <- p2 + ggplot2::theme(legend.position = "none")
    }
    p <- list("pca" = p, "umap" = p2)
  }

  if (returnData) {
    return(list("data" = pca.df, "plot" = p, "pca" = pca))
  } else {
    return(p)
  }
}
