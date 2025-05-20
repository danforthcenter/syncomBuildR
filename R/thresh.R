#' Function to run changepoint models of individual microbes vs phenotypes or
#' communities/clusters vs phenotypes.
#'
#'
#' @param asvTab An asv table or analogous dataframe with a row per observation and columns for traits.
#' @param phenoCols A character vector of column names for phenotypes to be used in changepoint models.
#' @param asvCols A character vector of columns representing microbes (predictor variables).
#' Defaults NULL where all column names containing the string "ASV" will be used.
#' @param model Type of changepoint model in chngpt::chngptm labeling convention. Currently hinge,
#' upperhinge, and segmented are supported. See Figure 2.1 of the chngpt
#' [vignette](https://cran.r-project.org/web/packages/chngpt/vignettes/chngpt-vignette.pdf).
#' @param cores Number of cores to run in parallel.
#' @param calibratePheno An optional vector of column names to calibrate the phenotypes by.
#' This should generally correspond to those used in `cal` if the ASV table has been calibrated or
#' just represent confounders that you wish to remove effects from in the changepoint regression.
#' @param p.adjust.method A method to adjust for multiple testing, defaults to "none". Available options
#' are shown with \code{stats::p.adjust.methods}.
#' \code{lm} is used to get residuals of the phenotype after these effects are regressed out.
#'
#' @keywords changepoint, threshold, regression, phenotype
#' @import chngpt
#' @return A \code{thresh} object.
#'
#' @examples
#'
#' asv$biomass_z <- rnorm(nrow(asv))
#' tm <- thresh(asv, "biomass_z")
#' tm
#'
#' @export

thresh <- function(asvTab, phenoCols, asvCols = NULL, model = "hinge",
                   cores = getOption("mc.cores", 1), calibratePheno = NULL,
                   p.adjust.method = "none") {
  if (is.null(asvCols)) {
    asvCols <- colnames(asvTab)[grepl("ASV", colnames(asvTab))]
  }
  if (!is.null(calibratePheno)) {
    for (phenotype in phenoCols) {
      formString <- paste0(phenotype, "~", paste0(calibratePheno, collapse = "+"))
      asvTab[[phenotype]] <- residuals(lm(as.formula(formString),
                                          data = asvTab, na.action = na.exclude))
    }
  }
  threshOut <- parallel::mclapply(asvCols, function(asv_col) {
    thresh_pheno <- lapply(phenoCols, function(phenotype) {
      if (model == "hinge" | model == "M01") {
        model <- "hinge"
        f1 <- as.formula(paste0(phenotype, "~1"))
        f2 <- as.formula(paste0("~", asv_col))
      } else if (model == "upperhinge" | model == "M10") {
        model <- "upperhinge"
        f1 <- as.formula(paste0(phenotype, "~", asv_col))
        f2 <- as.formula(paste0("~1"))
      } else if (model == "segmented" | model == "M11") {
        model <- "segmented"
        f1 <- as.formula(paste0(phenotype, "~1"))
        f2 <- as.formula(paste0("~", asv_col))
      }
      sub <- asvTab[, c(phenotype, asv_col)]
      tryCatch(
        {
          fit <- chngpt::chngptm(
            formula.1 = f1, formula.2 = f2, data = sub, type = model,
            family = "gaussian", est.method = "fastgrid", var.type = "bootstrap", save.boot = TRUE
          )
          return(fit)
        },
        error = function(err) {}
      )
    })
    unpacked_pheno <- .unpack_chngptm_proto_thresh(thresh_pheno, names = phenoCols)
    unpacked_pheno$model <- thresh_pheno
    return(unpacked_pheno)
  }, mc.cores = cores)
  #* `Parse output into thresh object`
  names(threshOut) <- asvCols
  thresh <- Reduce(.merge_proto_thresh, threshOut)
  #* p-value adjustment
  pvals <- unlist(lapply(thresh$slope, function(x) x$pval))
  adj_pvals <- p.adjust(pvals, method = p.adjust.method)
  for (i in seq_along(adj_pvals)) {
    thresh$slope[[i]]$padj <- adj_pvals[i]
  }
  #* assign other thresh slots
  thresh[["predictor"]] <- asvCols
  thresh[["data"]] <- asvTab[, c(phenoCols, asvCols)] # also stored in model$best.fit$data
  thresh[["type"]] <- "chngptm"
  thresh[["unit"]] <- "individual"
  thresh[["control"]] <- list("call" = match.call(),
                              "p.adjust.method" = p.adjust.method,
                              "subsettable" = c(
                                "intercept", "changepoint", "slope",
                                "phenotype", "model", "predictor"
                              ))
  thresh <- as.thresh(thresh)
  return(thresh)
}

#' Function to unpack list returned by main loop of chngptm model fitting.
#' 
#' Should take a list of models and return a list of components of those models.
#' 
#' @param proto_thresh a list of chngptm models
#' @param names a vector of names for the models (the phenotypes they were fit to, `phenoCols`)
#'
#' @return a list of components from the models, simplified for downstream use and for standardization
#' between eventual backends.
#' @keywords internal
#' @noRd

.unpack_chngptm_proto_thresh <- function(proto_thresh, names) {
  intercepts <- lapply(proto_thresh, function(pt) {
    coef(summary(pt))[1, 1]
  })
  changepoints <- lapply(proto_thresh, function(pt) {
    as.numeric(pt$chngpt)
  })
  slopes <- lapply(proto_thresh, function(pt) {
    coefs <- coef(summary(pt))
    list("est" = coefs[2, 1],
         "pval" = coefs[2, 5]
         )
  })
  #* data, predictor, type, etc will be set outside of this function,
  #* this is just to turn the list from being indexed on phenotypes to being
  #* indexed on ASV
  out <- list(
    intercept = intercepts,
    changepoint = changepoints,
    slope = slopes,
    phenotype = names
  )
  return(out)
}

#' Function to merge lists by name
#' 
#' This is used to combine outputs from many ASVs (predictors) over potentially many phenotypes (outcomes)
#' The output of this function should be a list of lists where the lists have length N_asvs
#' 
#' @param x element in the list (note, in practice this is used with Reduce)
#' @param y element in the list
#' @return A list of lists where each element in the list has 1 object per predictor (ASV).
#' @keywords internal
#' @noRd

.merge_proto_thresh <- function(x, y) {
  if (is.null(names(x)) | is.null(names(y))) {
    return(c(x, y)) # if we're at the bottom level of the hierarchy, return the values
  }
  keys <- union(names(x), names(y)) # if we are not at the bottom level of the hierarchy then grab
  # the shared names and apply the x/y objects again with only the shared names. Recurse until at
  # bottom of hierarchy.
  out <- sapply(keys, function(k) {
    .merge_proto_thresh(x[[k]], y[[k]])
    }, simplify = F)
  return(out)
}


