#' Function to run changepoint models of individual microbes vs phenotypes or
#' communities/clusters vs phenotypes.
#'
#'
#' @param x An asv table or analogous dataframe with a row per observation and columns for traits.
#' @param phenoCols A character vector of column names for phenotypes to be used in changepoint models.
#' @param predCols A character vector of columns representing microbes (predictor variables).
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

thresh <- function(x, phenoCols, predCols = NULL, model = "hinge",
                   cores = getOption("mc.cores", 1), calibratePheno = NULL,
                   p.adjust.method = "none") {
  UseMethod("thresh")
}

#' @rdname thresh
#' @examples
#' taxa <- c(
#'   "Bacteria", "Proteobacteria", "Betaproteobacteria", "Burkholderiales",
#'   "Burkholderiaceae", "Paraburkholderia", NA
#' )
#' taxa <- matrix(rep(taxa, 10), nrow = 10, byrow = TRUE)
#' colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
#' rownames(taxa) <- paste0("ASV", 1:10)
#' # taxonomy data if used should have ASV names explicitly as a column
#' taxa_df <- as.data.frame(taxa)
#' taxa_df$asv <- rownames(taxa_df)
#' sp_dist <- asvDist(asv, method = "spearman", clr_transform = TRUE, edgeFilter = 0.5)
#' net_data <- asvNet(sp_dist, taxa_df, edge = "spearman_similarity")
#' net_data <- netClust(net = net_data, "components")
#' asv$biomass_z <- rnorm(nrow(asv))
#' 
#' 
#' 
#' # ARGS
#' x <- net_data
#' phenoCols = NULL
#' predCols = NULL
#' model = "hinge"
#' cores = getOption("mc.cores", 1)
#' calibratePheno = NULL
#' p.adjust.method = "none"
#' 
#' @export
thresh.sbcnet <- function(x, phenoCols, predCols = NULL, model = "hinge",
                          cores = getOption("mc.cores", 1), calibratePheno = NULL,
                          p.adjust.method = "none") {
  #* here predCols needs to do something to find a default if it isn't given
  #* I am thinking that the default is the thresh on all columns of x$nodes that
  #* contains the string "cluster"?
  #* 
  #* What all needs to happen:
  #* Need to find the cluster columns that I want to group then use
  #* Need to find the members of each cluster, group those in the ASV data, and sum them per sample
  #* Need to run thresh model on that summed data
  #* Need to return thresh object.
  nodes <- x[["nodes"]]
  if (is.null(predCols)) {
    predCols <- colnames(nodes)[grepl("cluster", colnames(nodes), ignore.case = TRUE)]
  }
  if (is.null(cluster)) {
    cluster <- unique(nodes[[clusterCol]])
  }
  if (is.null(asvCols)) {
    asvCols <- "ASV"
  }
  if (length(asvCols) == 1) {
    asvCols <- colnames(asvTab)[grepl(asvCols, colnames(asvTab))]
  }
  if (is.null(clusterCol)) {
    clusterCol <- colnames(nodes)[grepl("cluster", colnames(nodes))][1]
  }
  #* `take nodes in a given cluster and aggregate a count in the asv table`
  clust_ag <- do.call(cbind, lapply(cluster, function(clust) {
    asvs_in_cluster <- nodes[nodes[[clusterCol]] == clust, "asv"]
    setNames(
      data.frame(
        rowSums(
          as.data.frame(asvTab[, c(asvs_in_cluster)])
        )
      ), c(paste0("cluster_", clust))
    )
  }))
  clusterColumns <- colnames(clust_ag)
  clust_ag <- cbind(asvTab[, -which(colnames(asvTab) %in% asvCols)], clust_ag)
  #* `calibrate phenotype by calibratePheno`
  if (!is.null(calibratePheno)) {
    for (phenotype in phenoCols) {
      formString <- paste0(phenotype, "~", paste0(calibratePheno, collapse = "+"))
      clust_ag[[phenotype]] <- residuals(lm(as.formula(formString),
                                            data = clust_ag, na.action = na.exclude))
    }
  }
  netThreshOut <- do.call(rbind, parallel::mclapply(clusterColumns, function(col) {
    thresh_df <- do.call(rbind, lapply(phenoCols, function(phenotype) {
      if (model == "hinge" | model == "M01") {
        model <- "hinge"
        f1 <- as.formula(paste0(phenotype, "~1"))
        f2 <- as.formula(paste0("~", col))
      } else if (model == "upperhinge" | model == "M10") {
        model <- "upperhinge"
        f1 <- as.formula(paste0(phenotype, "~", col))
        f2 <- as.formula(paste0("~1"))
      } else if (model == "segmented" | model == "M11") {
        model <- "segmented"
        f1 <- as.formula(paste0(phenotype, "~1"))
        f2 <- as.formula(paste0("~", col))
      }
      sub <- clust_ag[, c(phenotype, col)]
      tryCatch(
        {
          fit <- chngpt::chngptm(
            formula.1 = f1, formula.2 = f2, data = sub, type = model,
            family = "gaussian", est.method = "fastgrid", var.type = "bootstrap", save.boot = TRUE
          )
          out <- data.frame(coef(summary(fit)))
          out$Source <- rownames(out)
          rownames(out) <- NULL
          out$changePoint <- as.numeric(fit$chngpt)
          out$cluster <- col
          out$phenotype <- phenotype
          out$model <- model
          out$clusterType <- clusterCol
          out$model_id <- paste(clusterCol, col, phenotype, model, sep = "_")
          if (!is.null(calibratePheno)) {
            out$calibratePheno <- formString
          }
          out
        },
        warning = function(war) {},
        error = function(err) {}
      )
    }))
    return(thresh_df)
  }, mc.cores = cores))
  return(netThreshOut)
}

#' @rdname thresh
#' @export
thresh.data.frame <- function(x, phenoCols, predCols = NULL, model = "hinge",
                              cores = getOption("mc.cores", 1), calibratePheno = NULL,
                              p.adjust.method = "none") {
  if (is.null(predCols)) {
    predCols <- colnames(x)[grepl("ASV", colnames(x))]
  }
  if (!is.null(calibratePheno)) {
    for (phenotype in phenoCols) {
      formString <- paste0(phenotype, "~", paste0(calibratePheno, collapse = "+"))
      x[[phenotype]] <- residuals(lm(as.formula(formString),
                                          data = x, na.action = na.exclude))
    }
  }
  threshOut <- parallel::mclapply(predCols, function(pred_col) {
    thresh_pheno <- lapply(phenoCols, function(phenotype) {
      if (model == "hinge" | model == "M01") {
        model <- "hinge"
        f1 <- as.formula(paste0(phenotype, "~1"))
        f2 <- as.formula(paste0("~", pred_col))
      } else if (model == "upperhinge" | model == "M10") {
        model <- "upperhinge"
        f1 <- as.formula(paste0(phenotype, "~", pred_col))
        f2 <- as.formula(paste0("~1"))
      } else if (model == "segmented" | model == "M11") {
        model <- "segmented"
        f1 <- as.formula(paste0(phenotype, "~1"))
        f2 <- as.formula(paste0("~", pred_col))
      }
      sub <- x[, c(phenotype, pred_col)]
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
  names(threshOut) <- predCols
  thresh <- Reduce(.merge_proto_thresh, threshOut)
  #* p-value adjustment
  pvals <- unlist(lapply(thresh$slope, function(x) x$pval))
  adj_pvals <- p.adjust(pvals, method = p.adjust.method)
  for (i in seq_along(adj_pvals)) {
    thresh$slope[[i]]$padj <- adj_pvals[i]
  }
  #* assign other thresh slots
  thresh[["predictor"]] <- predCols
  thresh[["data"]] <- x[, c(phenoCols, predCols)] # also stored in model$best.fit$data
  thresh[["type"]] <- "chngptm"
  thresh[["unit"]] <- "individual"
  thresh[["control"]] <- list("call" = match.call(),
                              "p.adjust.method" = p.adjust.method,
                              "subsettable" = c(
                                "intercept", "changepoint", "slope",
                                "phenotype", "model", "predictor"
                              ),
                              "calibration" = calibratePheno
  )
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


