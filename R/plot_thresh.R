#' Plot a \code{thresh} object.
#'
#' @aliases plot.thresh
#'
#' @param x An object of class \code{thresh}.
#' @param predictors A vector of ASV column names. Defaults to NULL in which case all columns containing
#' "ASV" are used and a list of ggplots is returned.
#' @param outcomes A vector of phenotype names in \code{thresh}. Defaults to NULL where all phenotypes
#' are used and a list of plots is returned per ASV.
#' @param ... further arguments, ignored.
#' @import ggplot2
#' @import patchwork
#' @import viridis
#' @examples
#'
#' asv$biomass_z <- rnorm(nrow(asv))
#' tm <- thresh(asv, "biomass_z")
#' plot(tm, "ASV9")
#'
#' @method plot thresh
#' @export

plot.thresh <- function(x, predictors = NULL, outcomes = NULL, ...) {
  if (is.null(predictors)) {
    predictors <- unique(x$predictor)
  }
  x <- x[x$predictor %in% predictors]
  d <- x$data
  outList <- lapply(predictors, function(microbe) {
    outcomes_iter <- outcomes
    x_iter <- x[x$predictor == microbe]
    if (is.null(outcomes_iter)) {
      outcomes_iter <- unique(x_iter$phenotype)
    }
    phenoPlots <- lapply(outcomes_iter, function(pheno) {
      postCptCol <- if (x_iter$slope[[1]]$padj < 0.05) {
        viridis::plasma(1, begin = 0.7)
      } else {
        "black"
      }
      p <- ggplot2::ggplot(d) +
        ggplot2::geom_point(
          data = d[d[[x_iter$predictor]] <= x_iter$changepoint[[1]], ],
          ggplot2::aes(x = .data[[x_iter$predictor]], y = .data[[x_iter$phenotype]]),
          color = "gray40", size = 2, alpha = 0.5
        ) +
        ggplot2::geom_point(
          data = d[d[[x_iter$predictor]] > x_iter$changepoint[[1]], ],
          ggplot2::aes(x = .data[[x_iter$predictor]], y = .data[[x_iter$phenotype]]),
          color = postCptCol, size = 2, alpha = 0.85
        ) +
        ggplot2::geom_vline(xintercept = x_iter$changepoint[[1]], linetype = 5, color = "black") +
        ggplot2::geom_segment(
          x = min(d[[x_iter$predictor]]),
          xend = x_iter$changepoint[[1]],
          y = x_iter$intercept[[1]],
          yend = x_iter$intercept[[1]],
          color = "black"
        ) +
        geom_segment(
          x = x_iter$changepoint[[1]],
          xend = max(d[[x_iter$predictor]]),
          y = x_iter$intercept[[1]],
          yend = x_iter$intercept[[1]] + x_iter$slope[[1]]$est *
            (max(d[[x_iter$predictor]]) - x_iter$changepoint[[1]]),
          color = postCptCol
        ) +
        ggplot2::labs(
          y = paste0(
            ifelse(is.null(x_iter$control$calibration), "", "Calibrated "), pheno
          ),
          x = microbe,
          title = x_iter$predictor,
          subtitle = paste0("P-value of slope: ", round(x_iter$slope[[1]]$padj, 3))
        ) +
        theme_light() +
        theme(
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 9)
        )
      return(p)
    })
    if (length(outcomes_iter) > 1) {
      phenoPlots <- Reduce(phenoPlots, `+`) +
        patchwork::plot_layout(
          guides = "collect",
          axes = "collect",
          axis_titles = "collect"
        )
    } else {
      phenoPlots <- phenoPlots[[1]]
    }
    return(phenoPlots)
  })

  if (length(predictors) == 1) {
    outList <- outList[[1]]
  }
  return(outList)
}
