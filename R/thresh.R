#' Function to run changepoint models of individual microbes vs phenotypes or communities/clusters vs phenotypes.
#'
#'
#' @param asvTab An asv table or analogous dataframe with a row per observation and columns for traits.
#' @param phenoCols A character vector of column names for phenotypes to be used in changepoint models.
#' @param asvCols A character vector of columns representing microbes (predictor variables). Defaults NULL where all column names containing the string "ASV" will be used.
#' @param model Type of changepoint model in chngpt::chngptm labeling convention. Currently hinge, upperhinge, and segmented are supported. See Figure 2.1 of the \link{[chngpt vignette](https://cran.r-project.org/web/packages/chngpt/vignettes/chngpt-vignette.pdf)}
#' @param cores Number of cores to run in parallel.
#' @param calibratePheno An optional vector of column names to calibrate the phenotypes by. This should generally correspond to those used in `cal` if the ASV table has been calibrated or just represent confounders that you wish to remove effects from in the changepoint regression. /code{lm} is used to get residuals of the phenotype after these effects are regressed out.
#'
#' @keywords changepoint, threshold, regression, phenotype
#' @import chngpt
#' @return A dataframe summarizing changepoint models for individual ASVs vs phenotypes.
#'
#' @examples
#'
#' # a<-qc(); b<-cal(a); c<-thresh(b)
#' print(load("/home/jsumner/Desktop/stargate/SINC/sincUtils/syncomBuilder/cal_output.rdata"))
#' asv <- are_c[[1]]
#' zinbCalibrated = are_c[[2]][are_c[[2]]$model == "ZINB", "asv"]
#' asv <- are_c[[1]][, c("tissue", "plot", "row", "genotype", "biomass", "sd", zinbCalibrated)]
#'
#' threshMods <- thresh(asvTab = asv, phenoCols = "biomass", model = "hinge", cores = 10, calibratePheno = "genotype")
#' dim(threshMods)
#' head(threshMods)
#'
#' @export

thresh <- function(asvTab, phenoCols, asvCols = NULL, model = "hinge", cores = getOption("mc.cores", 1), calibratePheno = NULL) {
  if (is.null(asvCols)) {
    asvCols <- colnames(asvTab)[grepl("ASV", colnames(asvTab))]
  }
  threshOut <- do.call(rbind, lapply(phenoCols, function(phenotype) {
    print(phenotype)
    if (!is.null(calibratePheno)) {
      formString <- paste0(phenotype, "~", paste0(calibratePheno, collapse = "+"))
      asvTab[[phenotype]] <- residuals(lm(as.formula(formString), data = asvTab, na.action = na.exclude))
    } else {
      asvTab[[phenotype]] <- asvTab[[phenotype]]
    }
    thresh_df <- do.call(rbind, parallel::mclapply(asvCols, function(asv_col) {
      if (model == "hinge" | model == "M01") {
        model = "hinge"
        f1 <- as.formula(paste0(phenotype, "~1"))
        f2 <- as.formula(paste0("~", asv_col))
      } else if (model == "upperhinge" | model == "M10") {
        model = "upperhinge"
        f1 <- as.formula(paste0(phenotype, "~", asv_col))
        f2 <- as.formula(paste0("~1"))
      } else if (model == "segmented" | model == "M11") {
        model = "segmented"
        f1 <- as.formula(paste0(phenotype, "~1"))
        f2 <- as.formula(paste0("~", asv_col))
      }

      sub <- asvTab[, c(phenotype, asv_col)]
      tryCatch(
        {
          fit <- chngpt::chngptm(
            formula.1 = f1, formula.2 = f2, data = sub, type = model,
            family = "gaussian", est.method = "fastgrid", var.type = "bootstrap", save.boot = T
          )
          out <- data.frame(coef(summary(fit)))
          out$Source <- rownames(out)
          rownames(out) <- NULL
          out$changePoint <- as.numeric(fit$chngpt)
          out$asv <- asv_col
          out$phenotype <- phenotype
          out$model <- model
          if (!is.null(calibratePheno)) {
            out$calibratePheno = formString
          }
          out
        },
        warning = function(war) {},
        error = function(err) {}
      )
    }, mc.cores = cores))
    thresh_df
  }))
  return(threshOut)
}
