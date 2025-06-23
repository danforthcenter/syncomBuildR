#' Function to plot results of a changepoint models as upset plots
#'
#' @param thresh Output from \code{\link{thresh}}
#' @param cutoff Alpha, the significance threshold. Defaults to 0.05.
#'
#' @keywords changepoint, threshold, regression, phenotype, upset
#'
#' @import patchwork
#' @importFrom ComplexUpset upset
#'
#' @return A ggplot or list of ggplots showing changepoint models against some set of phenotypes.
#'
#' @examples
#'
#' set.seed(123)
#' asv$biomass_z <- rnorm(nrow(asv))
#' asv$height_z <- rnorm(nrow(asv))
#' tm <- thresh(asv, c("biomass_z", "height_z"))
#'
#' threshUpset(tm, cutoff = 0.5)
#'
#' @export

threshUpset <- function(thresh, cutoff = 0.05) {
  sig <- unlist(lapply(unique(thresh$predictor), function(pred) {
    sub <- thresh[thresh$predictor == pred]
    if (sub$type == "chngptm") {
      sig_index <- which(sub$pval < cutoff)
    } else if (sub$type == "mcp") {
      sig_index <- which(sub$post_prob > cutoff)
    }
    if (length(sig_index) > 0) {
      return(pred)
    }
  }))
  if (is.null(sig)) {
    stop("No data with p past cutoff")
  }
  d <- as.data.frame(thresh[thresh$predictor %in% sig][c("predictor", "phenotype")])
  if (thresh$type == "chngptm") {
    d$sig <- thresh[thresh$predictor %in% sig]$pval < cutoff
  } else if (thresh$type == "mcp") {
    d$sig <- thresh[thresh$predictor %in% sig]$post_prob > cutoff
  }
  data.table::setDT(d)
  upsetPlotData <- as.data.frame(data.table::dcast(d,
    predictor ~ phenotype,
    value.var = "sig"
  ))
  p <- ComplexUpset::upset(upsetPlotData, intersect = colnames(upsetPlotData)[-1])
  p[[2]] <- p[[2]] +
    ggplot2::labs(
      title = paste0("ASV ~ Phenotype Correlations"),
      subtitle = paste0(
        ifelse(thresh$type == "chngptm", "With P.value < ", "With Post. Prob > "),
        cutoff
      )
    ) +
    ggplot2::theme(
      legend.title.align = NULL,
      legend.text.align = NULL
    )
  return(p)
}
