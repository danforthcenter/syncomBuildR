#' Function to run bayesian changepoint models of individual microbes vs phenotypes or
#' communities/clusters vs phenotypes.
#'
#'
#' @param x An asv table or analogous dataframe with a row per observation and columns for traits.
#' Alternatively an \code{scbnet} object.
#' @param phenoCols A character vector of column names for phenotypes to be used in changepoint models.
#' @param predCols A character vector of columns representing microbes (predictor variables).
#' Defaults NULL where all column names containing the string "ASV" will be used.
#' @param model Changepoint model formula as used in mcp::mcp.
#' See [mcp documentation](https://lindeloev.github.io/mcp/articles/formulas.html) for examples.
#' If you are using multiple phenoCols or multiple predCols specify Y and X as outcome/predictor.
#' @param cores Number of cores to run in parallel.
#' @param calibratePheno An optional vector of column names to calibrate the phenotypes by.
#' This should generally correspond to those used in `cal` if the ASV table has been calibrated or
#' just represent confounders that you wish to remove effects from in the changepoint regression.
#' \code{lm} is used to get residuals of the phenotype after these effects are regressed out.
#' @param keep_models Should full model objects be kept? This can make thresh objects much larger.
#' Defaults to FALSE.
#' @param hypothesis A hypothesis to test, defaults to "X > 0" which uses the same logic as
#' the model formula to test phenoCols.
#' @param ... Additional arguments passed to methods.
#'
#' @keywords changepoint, threshold, regression, phenotype
#' @import chngpt
#' @return A \code{thresh} object.
#'
#' @examples
#'
#' asv$biomass_z <- rnorm(nrow(asv))
#' tm <- mcp_thresh(asv, "biomass_z",
#'   model = list(
#'     Y ~ 1, # intercept,
#'     ~ 0 + X # joined slope after first changepoint
#'   )
#' )
#' tm
#'
#' @export

mcp_thresh <- function(
    x, phenoCols, predCols = NULL,
    model = list(Y ~ 1, ~ 0 + X),
    cores = getOption("mc.cores", 1), calibratePheno = NULL,
    keep_models = FALSE, hypothesis = "X > 0", ...) {
  UseMethod("mcp_thresh")
}

#' @rdname mcp_thresh
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
#' tn <- mcp_thresh(net_data, phenoCols = "biomass_z", asvTab = asv)
#' tn
#'
#' @method mcp_thresh scbnet
#' @importFrom mcp mcp
#' @export
mcp_thresh.scbnet <- function(
    x, phenoCols, predCols = NULL, model = list(Y ~ 1, ~ 0 + X),
    cores = getOption("mc.cores", 1), calibratePheno = NULL,
    keep_models = FALSE,
    hypothesis = "X > 0",
    asvTab = NULL, ...) {
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
      clust_ag[[phenotype]] <- residuals(
        lm(
          as.formula(formString),
          data = clust_ag, na.action = na.exclude
        )
      )
    }
  }
  proto_thresh <- parallel::mclapply(names(clusterColumns), function(pred_col) {
    thresh_pred_col <- lapply(clusterColumns[[pred_col]], function(clust) {
      thresh_pheno <- lapply(phenoCols, function(phenotype) {
        sub <- clust_ag[, c(phenotype, clust)]
        tryCatch(
          { # mcp
            model_char <- lapply(model, deparse1) # convert to character prettily
            model_char_iter <- lapply(model_char, \(s) {
              s <- gsub("\\bX\\b", phenotype, s)
              return(s)
            })
            model_char_iter <- lapply(model_char_iter, \(s) {
              s <- gsub("\\bY\\b", clust, s)
              return(s)
            })
            model_iter <- lapply(model_char_iter, as.formula)
            fit <- mcp::mcp(
              model = model_iter,
              data = sub,
              cores = cores,
              chains = max(c(2, cores)),
              iter = 2000, adapt = 1000
            )
            fit <- fit[names(fit) != "mcmc_loglik"] # drop loglik, not using it and it's huge.
            class(fit) <- "mcpfit"
            return(fit)
          },
          error = function(err) {}
        )
      })
      # unpack models of each phenotype per a cluster within a clustering scheme
      names_thresh_pheno <- phenoCols[which(!unlist(lapply(thresh_pheno, is.null)))]
      thresh_pheno <- thresh_pheno[which(!unlist(lapply(thresh_pheno, is.null)))]
      unpacked_pheno <- .unpack_chngptm_mcp_thresh(thresh_pheno,
        names = names_thresh_pheno,
        pred = pred_col,
        hyp = hypothesis
      )
      unpacked_pheno$predictor <- pred_col
      if (keep_models) {
        unpacked_pheno$model <- thresh_pheno
      } else {
        unpacked_pheno$model <- lapply(thresh_pheno, function(tp) {
          tp$mcmc_post <- utils::head(tp$mcmc_post, n = 25)
          return(tp)
        })
      }
      return(unpacked_pheno)
    })
    names(thresh_pred_col) <- clusterColumns[[pred_col]]
    reduced_thresh_pc <- Reduce(.merge_proto_thresh, thresh_pred_col)
    return(reduced_thresh_pc)
  }, mc.cores = min(length(clusterColumns), cores))
  full_proto_thresh <- Reduce(.merge_proto_thresh, proto_thresh)
  #* add other slots
  full_proto_thresh[["data"]] <- clust_ag[, c(phenoCols, unname(unlist(clusterColumns)))]
  full_proto_thresh[["type"]] <- "mcp"
  full_proto_thresh[["unit"]] <- "cluster"
  full_proto_thresh[["control"]] <- list(
    "call" = match.call(),
    "subsettable" = c(
      "post_means", "slope", "post_probs", "hypotheses",
      "phenotype", "model", "predictor"
    ),
    "calibration" = calibratePheno
  )
  thresh <- as.thresh(full_proto_thresh)
  return(thresh)
}

#' @method mcp_thresh data.frame
#' @importFrom mcp mcp
#' @export
mcp_thresh.data.frame <- function(x, phenoCols, predCols = NULL,
                                  model = list(Y ~ 1, ~ 0 + X),
                                  cores = getOption("mc.cores", 1),
                                  calibratePheno = NULL,
                                  keep_models = FALSE,
                                  hypothesis = "X > 0",
                                  ...) {
  if (is.null(predCols)) {
    predCols <- colnames(x)[grepl("ASV", colnames(x))]
  }
  if (!is.null(calibratePheno)) {
    for (phenotype in phenoCols) {
      formString <- paste0(phenotype, "~", paste0(calibratePheno, collapse = "+"))
      x[[phenotype]] <- residuals(
        lm(
          as.formula(formString),
          data = x, na.action = na.exclude
        )
      )
    }
  }
  threshOut <- parallel::mclapply(predCols, function(pred_col) {
    thresh_pheno <- lapply(phenoCols, function(phenotype) {
      sub <- x[, c(phenotype, pred_col)]
      tryCatch(
        {
          model_char <- lapply(model, deparse1)
          model_char_iter <- lapply(model_char, \(s) {
            s <- gsub("\\bX\\b", phenotype, s)
            return(s)
          })
          model_char_iter <- lapply(model_char_iter, \(s) {
            s <- gsub("\\bY\\b", pred_col, s)
            return(s)
          })
          model_iter <- lapply(model_char_iter, as.formula)
          fit <- mcp::mcp(
            model = model_iter,
            data = sub,
            cores = cores,
            chains = max(c(2, cores)),
            iter = 2000, adapt = 1000
          )
          fit <- fit[names(fit) != "mcmc_loglik"] # drop loglik, not using it and it's huge.
          class(fit) <- "mcpfit"
          return(fit)
        },
        error = function(err) {}
      )
    })
    names_thresh_pheno <- phenoCols[which(!unlist(lapply(thresh_pheno, is.null)))]
    thresh_pheno <- thresh_pheno[which(!unlist(lapply(thresh_pheno, is.null)))]
    unpacked_pheno <- .unpack_chngptm_mcp_thresh(thresh_pheno,
      names = names_thresh_pheno,
      pred = pred_col,
      hyp = hypothesis
    )
    if (keep_models) {
      unpacked_pheno$model <- thresh_pheno
    } else {
      unpacked_pheno$model <- lapply(thresh_pheno, function(tp) {
        tp$mcmc_post <- utils::head(tp$mcmc_post, n = 25)
        return(tp)
      })
    }
    return(unpacked_pheno)
  }, mc.cores = cores)
  #* `Parse output into thresh object`
  names(threshOut) <- predCols
  thresh <- Reduce(.merge_proto_thresh, threshOut)
  #* assign other thresh slots
  thresh[["data"]] <- x[, c(phenoCols, predCols)]
  thresh[["type"]] <- "mcp"
  thresh[["unit"]] <- "individual"
  thresh[["control"]] <- list(
    "call" = match.call(),
    "subsettable" = c(
      "post_means", "slope", "post_probs", "hypotheses",
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
#' @param pred predictor variables as a vector
#' @param hyp Hypothesis as a character string
#'
#' @return a list of components from the models, simplified for downstream use and for standardization
#' between eventual backends.
#' @keywords internal
#' @noRd

.unpack_chngptm_mcp_thresh <- function(proto_thresh, names, pred, hyp) {
  means <- lapply(proto_thresh, function(pt) {
    summary(pt$mcmc_post)[[1]][, 1]
  })
  hypotheses <- lapply(seq_along(proto_thresh), function(i) {
    pt <- proto_thresh[[i]]
    var_names <- names(summary(pt$mcmc_post)[[1]][, 1])
    outcome <- names[[i]]
    outcome_var <- var_names[which(grepl(outcome, var_names))]
    hyp1 <- gsub("\\bX\\b", outcome_var, hyp)
    return(
      mcp::hypothesis(pt, hyp1)
    )
  })
  post_probs <- unlist(lapply(hypotheses, function(h) {
    h$p
  }))
  hypotheses_string <- unlist(lapply(hypotheses, function(h) {
    h$hypothesis
  }))
  out <- list(
    post_means = means,
    post_probs = post_probs,
    hypotheses = hypotheses_string,
    phenotype = names,
    predictor = rep(pred, length(names))
  )
  return(out)
}
