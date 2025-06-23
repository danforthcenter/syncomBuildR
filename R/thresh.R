#' Function to run changepoint models of individual microbes vs phenotypes or
#' communities/clusters vs phenotypes.
#'
#'
#' @param x An asv table or analogous dataframe with a row per observation and columns for traits.
#' Alternatively an \code{scbnet} object.
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
#' \code{lm} is used to get residuals of the phenotype after these effects are regressed out.
#' @param p.adjust.method A method to adjust for multiple testing, defaults to "none". Available options
#' are shown with \code{stats::p.adjust.methods}.
#' @param keep_models Should full model objects be kept? This can make thresh objects much larger.
#' Defaults to FALSE.
#' @param ... Additional arguments passed to methods.
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
                   p.adjust.method = "none", keep_models = FALSE,  ...) {
  UseMethod("thresh")
}

#' @rdname thresh
#' @param asvTab An ASV table used in making the \code{scbnet} object.
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
#' net_data <- netClust(net = net_data, pullNode(net_data,
#'   node = "ASV10",
#'   edge = "spearman_similarity", edgeFilter = 0.7, nodeCol = "asv"
#' ))
#' net_data <- netClust(net = net_data, pullNode(net_data,
#'   node = "ASV1",
#'   edge = "spearman_similarity", edgeFilter = 0.7, nodeCol = "asv"
#' ))
#' asv$biomass_z <- rnorm(nrow(asv))
#' tn <- thresh(net_data, phenoCols = "biomass_z", asvTab = asv)
#'
#'
#' @method thresh scbnet
#' @export
thresh.scbnet <- function(x, phenoCols, predCols = NULL, model = "hinge",
                          cores = getOption("mc.cores", 1), calibratePheno = NULL,
                          p.adjust.method = "none", keep_models = FALSE, asvTab = NULL, ...) {
  nodes <- x[["nodes"]]
  if (is.null(predCols)) {
    predCols <- colnames(nodes)[grepl("cluster", colnames(nodes), ignore.case = TRUE)]
  }
  #* Need to find the members of each cluster, group those in the ASV data, and sum them per sample
  #* `Make aggregated count table`
  # for every pred_col in predCols (cluster columns)
  clust_ag_list <- lapply(predCols, function(pred_col) {
    # get all the levels of it, and for each cluster (level)
    col_df <- do.call(cbind, lapply(unique(nodes[[pred_col]]), function(clust) {
      # find the ASVs in that cluster
      asvs_in_cluster <- nodes[nodes[[pred_col]] == clust, "asv"]
      # make a column of those ASVs summed across samples, name it, and return it to be cbound
      setNames(
        data.frame(
          rowSums(
            as.data.frame(asvTab[, c(asvs_in_cluster)])
          )
        ), c(paste0(pred_col, "_", clust))
      )
    }))
    return(col_df)
  })
  # bind all into data.frame
  clust_ag <- do.call(cbind, clust_ag_list)
  # store names for use later
  clusterColumns <- lapply(clust_ag_list, colnames)
  names(clusterColumns) <- predCols
  # add metadata/phenotype data from asvTab
  clust_ag <- cbind(asvTab[, !grepl("ASV", colnames(asvTab))], clust_ag)
  #* `calibrate phenotype by calibratePheno`
  if (!is.null(calibratePheno)) {
    for (phenotype in phenoCols) {
      formString <- paste0(phenotype, "~", paste0(calibratePheno, collapse = "+"))
      clust_ag[[phenotype]] <- residuals(lm(as.formula(formString),
        data = clust_ag, na.action = na.exclude
      ))
    }
  }
  #* Need to run thresh model on that summed data
  proto_thresh <- parallel::mclapply(names(clusterColumns), function(pred_col) {
    thresh_pred_col <- lapply(clusterColumns[[pred_col]], function(clust) {
      thresh_pheno <- lapply(phenoCols, function(phenotype) {
        if (model == "hinge" | model == "M01") {
          model <- "hinge"
          f1 <- as.formula(paste0(phenotype, "~1"))
          f2 <- as.formula(paste0("~", clust))
        } else if (model == "upperhinge" | model == "M10") {
          model <- "upperhinge"
          f1 <- as.formula(paste0(phenotype, "~", clust))
          f2 <- as.formula(paste0("~1"))
        } else if (model == "segmented" | model == "M11") {
          model <- "segmented"
          f1 <- as.formula(paste0(phenotype, "~1"))
          f2 <- as.formula(paste0("~", clust))
        }
        sub <- clust_ag[, c(phenotype, clust)]
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
      # unpack models of each phenotype per a cluster within a clustering scheme
      names_thresh_pheno <- phenoCols[which(!unlist(lapply(thresh_pheno, is.null)))]
      thresh_pheno <- thresh_pheno[which(!unlist(lapply(thresh_pheno, is.null)))]
      unpacked_pheno <- .unpack_chngptm_proto_thresh(thresh_pheno, names = names_thresh_pheno,
                                                     pred = pred_col)
      unpacked_pheno$predictor <- pred_col
      if (keep_models) {
        unpacked_pheno$model <- thresh_pheno
      }
      return(unpacked_pheno)
    })
    thresh_pred_col
    names(thresh_pred_col) <- clusterColumns[[pred_col]]
    reduced_thresh_pc <- Reduce(.merge_proto_thresh, thresh_pred_col)
    return(reduced_thresh_pc)
  }, mc.cores = min(length(clusterColumns), cores))
  full_proto_thresh <- Reduce(.merge_proto_thresh, proto_thresh)
  #* p-value adjustment
  full_proto_thresh$pval <- p.adjust(full_proto_thresh$pval, method = p.adjust.method)
  #* add other slots
  full_proto_thresh[["data"]] <- clust_ag[, c(phenoCols, unname(unlist(clusterColumns)))]
  full_proto_thresh[["type"]] <- "chngptm"
  full_proto_thresh[["unit"]] <- "cluster"
  full_proto_thresh[["control"]] <- list(
    "call" = match.call(),
    "p.adjust.method" = p.adjust.method,
    "subsettable" = c(
      "intercept", "changepoint", "slope", "pval",
      "phenotype", "model", "predictor"
    ),
    "calibration" = calibratePheno
  )
  thresh <- as.thresh(full_proto_thresh)
  return(thresh)
}

#' @method thresh data.frame
#' @export
thresh.data.frame <- function(x, phenoCols, predCols = NULL, model = "hinge",
                              cores = getOption("mc.cores", 1), calibratePheno = NULL,
                              p.adjust.method = "none", keep_models = FALSE, ...) {
  if (is.null(predCols)) {
    predCols <- colnames(x)[grepl("ASV", colnames(x))]
  }
  if (!is.null(calibratePheno)) {
    for (phenotype in phenoCols) {
      formString <- paste0(phenotype, "~", paste0(calibratePheno, collapse = "+"))
      x[[phenotype]] <- residuals(lm(as.formula(formString),
        data = x, na.action = na.exclude
      ))
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
    names_thresh_pheno <- phenoCols[which(!unlist(lapply(thresh_pheno, is.null)))]
    thresh_pheno <- thresh_pheno[which(!unlist(lapply(thresh_pheno, is.null)))]
    unpacked_pheno <- .unpack_chngptm_proto_thresh(thresh_pheno, names = names_thresh_pheno,
                                                   pred = pred_col)
    if (keep_models) {
      unpacked_pheno$model <- thresh_pheno
    }
    return(unpacked_pheno)
  }, mc.cores = cores)
  #* `Parse output into thresh object`
  names(threshOut) <- predCols
  thresh <- Reduce(.merge_proto_thresh, threshOut)
  #* p-value adjustment
  thresh$pval <- p.adjust(thresh$pval, method = p.adjust.method)
  #* assign other thresh slots
  thresh[["data"]] <- x[, c(phenoCols, predCols)] # also stored in model$best.fit$data
  thresh[["type"]] <- "chngptm"
  thresh[["unit"]] <- "individual"
  thresh[["control"]] <- list(
    "call" = match.call(),
    "p.adjust.method" = p.adjust.method,
    "subsettable" = c(
      "intercept", "changepoint", "slope", "pval",
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

.unpack_chngptm_proto_thresh <- function(proto_thresh, names, pred) {
  intercepts <- lapply(proto_thresh, function(pt) {
    coef(summary(pt))[1, 1]
  })
  changepoints <- lapply(proto_thresh, function(pt) {
    as.numeric(pt$chngpt)
  })
  slopes <- lapply(proto_thresh, function(pt) {
    coefs <- coef(summary(pt))
    list(
      "est" = coefs[2, 1],
      "pval" = coefs[2, 5]
    )
  })
  out <- list(
    intercept = intercepts,
    changepoint = changepoints,
    slope = unlist(lapply(slopes, function(s) {s$est})),
    pval = unlist(lapply(slopes, function(s) {s$pval})),
    phenotype = names,
    predictor = rep(pred, length(names))
  )
  return(out)
}

#' Function to merge lists by name
#'
#' This is used to combine outputs from many ASVs (predictors) over potentially many
#' phenotypes (outcomes). The output of this function should be a list of lists where the lists
#' have length N_asvs
#'
#' @param x element in the list (note, in practice this is used with Reduce)
#' @param y element in the list
#' @return A list of lists where each element in the list has 1 object per predictor (ASV).
#' @keywords internal
#' @noRd

.merge_proto_thresh <- function(x, y) {
  if (is.null(names(x)) || is.null(names(y))) {
    return(c(x, y)) # if we're at the bottom level of the hierarchy, return the values
  }
  keys <- union(names(x), names(y)) # if we are not at the bottom level of the hierarchy then grab
  # the shared names and apply the x/y objects again with only the shared names. Recurse until at
  # bottom of hierarchy.
  out <- sapply(keys, function(k) {
    .merge_proto_thresh(x[[k]], y[[k]])
  }, simplify = FALSE)
  return(out)
}
