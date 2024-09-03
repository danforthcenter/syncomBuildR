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
#' tm <- thresh(asv, "biomass_z")
#' tm[11:20, "phenotype"] <- "not_biomass_z"
#' threshUpset(tm, cutoff = 0.5)
#'
#' @export

threshUpset <- function(thresh, mode = "phenotype",
                        pvalcol = "p.value.", cutoff = 0.05, cores = getOption("mc.cores", 1)) {
  if (mode == "phenotype") {
    p <- .pheno_upset_plot(thresh, pvalcol, cutoff, cores)
  } else {
    p <- .pvalue_upset_plot(thresh, pvalcol)
  }
  return(p)
}


.pheno_upset_plot <- function(thresh, pvalcol = "p.value.", cutoff = 0.05, cores = 1) {
  sig <- do.call(rbind, parallel::mclapply(unique(thresh$asv), function(asv){
    sub <- thresh[!grepl("Intercept", thresh$Source) & thresh$asv == asv, ]
    if(any(sub[[pvalcol]] < cutoff)){
      return(sub)
    } else{NULL}
  }, mc.cores=cores))
  if (is.null(sig)) {
    stop("No data with p.values below cutoff")
  }
  sub_upsetData <- sig[, c("phenotype", pvalcol, "asv")]
  sub_upsetData[[pvalcol]] <- ifelse(sub_upsetData[[pvalcol]] < cutoff, TRUE, FALSE)
  data.table::setDT(sub_upsetData)
  upsetPlotData <- as.data.frame(data.table::dcast(sub_upsetData, asv ~ phenotype, value.var = pvalcol))
  p <- ComplexUpset::upset(upsetPlotData, intersect = colnames(upsetPlotData)[-1])
  p[[2]] <- p[[2]] + 
    ggplot2::labs(title = paste0("ASV ~ Phenotype Correlations"),
                  subtitle = paste0("With ", pvalcol, " < ", cutoff)) +
    ggplot2::theme(legend.title.align = NULL,
                   legend.text.align = NULL)
  return(p)
}

#' @keywords internal
#' @noRd

.pvalue_upset_plot <- function(thresh, pvalcol = "p.value.") {
  df <- do.call(rbind, lapply(split(thresh, thresh$phenotype), function(d){
    d$p.adj.bonferroni <- p.adjust(d[[pvalcol]], method = "bonferroni")
    d$p.adj.holm <- p.adjust(d[[pvalcol]], method = "holm")
    d$p.adj.hochberg <- p.adjust(d[[pvalcol]], method = "hochberg")
    d$p.adj.fdr <- p.adjust(d[[pvalcol]], method = "fdr")
    d$p.adj.by <- p.adjust(d[[pvalcol]], method = "BY")
    d$p.adj.none <- p.adjust(d[[pvalcol]], method = "none")
    return(d)
  }))
  
  adj_results <- as.data.frame(do.call(cbind,
                                       lapply(df[, grepl("p.adj", colnames(df))], function(col){
                                         col<0.05
                                       })))
  p0 <- ComplexUpset::upset(adj_results, intersect = colnames(adj_results))
  p0[[2]] <- p0[[2]] +
    ggplot2::labs(title = paste0("P adjustment options")) +
    ggplot2::theme(legend.text.align = NULL)
  return(p0)
}