#' Function to plot results of a changepoint models as upset plots
#'
#' @param thresh Output from \code{\link{thresh}}
#' @param mode Either "phenotype" (the default) or "pvalue". "phenotype" will show an upset plot
#' of significant changepoint models against the phenotypes present in thresh, "pvalue" will show
#' the results of various p-value adjustment methods.
#' @param pvalcol The column in thresh that has p-value data, defaults to "p.value." which reflects
#' standard output from \code{thresh}.
#' @param cutoff Alpha, the significance threshold. Defaults to 0.05.
#' @param cores Number of cores to use in parallel. This should not make much of a difference generally
#' but can be useful with very large tables. Defaults to 1 if the "mc.cores" option is not set.
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
    sig_index <- which(unlist(lapply(sub$slope, function(x) x$padj)) < cutoff)
    if (length(sig_index) > 0) {
      #d <- data.frame(predictor = asv, pheno = sub$phenotype[sig_index])
      return(pred)
    }
  }))
  if (is.null(sig)) {
    stop("No data with p.values below cutoff")
  }
  d <- as.data.frame(thresh[thresh$predictor %in% sig][c("predictor", "phenotype")])
  d$sig <- unlist(lapply(thresh[thresh$predictor %in% sig]$slope, function(i) {
    i$padj
  }))
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

