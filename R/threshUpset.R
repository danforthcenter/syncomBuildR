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
    sig_index <- which(sub$pval < cutoff)
    if (length(sig_index) > 0) {
      return(pred)
    }
  }))
  if (is.null(sig)) {
    stop("No data with p.values below cutoff")
  }
  d <- as.data.frame(thresh[thresh$predictor %in% sig][c("predictor", "phenotype")])
  d$sig <- thresh[thresh$predictor %in% sig]$pval
  d$sig <- ifelse(d$sig < cutoff, TRUE, FALSE)
  data.table::setDT(d)
  upsetPlotData <- as.data.frame(data.table::dcast(d,
    predictor ~ phenotype,
    value.var = "sig"
  ))
  p <- ComplexUpset::upset(upsetPlotData, intersect = colnames(upsetPlotData)[-1])
  p[[2]] <- p[[2]] +
    ggplot2::labs(
      title = paste0("ASV ~ Phenotype Correlations"),
      subtitle = paste0("With P.value < ", cutoff)
    ) +
    ggplot2::theme(
      legend.title.align = NULL,
      legend.text.align = NULL
    )
  return(p)
  return(p)
}
