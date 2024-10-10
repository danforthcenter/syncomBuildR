#' Function to aggregate clustered networks overall phenotypic impacts.
#'
#'
#' @param net Object returned from \code{\link{asvNet}} with clusters added using \code{\link{netClust}}
#' @param asvTab An asv table with phenotypes joined.
#' @param asvCols A character vector of columns representing microbes (nodes).
#' Defaults NULL where all column names containing the string "ASV" will be used.
#' @param clusterCol The column name to use for clusters. If NULL then the first column name containing
#' the string "cluster" is used.
#' @param cluster A vector of clusters to be regressed against. By default this is NULL and all
#' clusters will be used.
#' @param phenoCols A vector of phenotype columns to use in changepoint regression.
#' These should be column names in the ASV table
#' @param model Type of changepoint model in chngpt::chngptm labeling convention. Currently hinge,
#' upperhinge, and segmented are supported. See Figure 2.1 of the chngpt
#' [vignette](https://cran.r-project.org/web/packages/chngpt/vignettes/chngpt-vignette.pdf)
#' @param cores Number of cores to run in parallel, defaults to 1 if "mc.cores" option is not set.
#' @param calibratePheno An optional vector of column names to calibrate the phenotypes by.
#' This should generally correspond to those used in `cal` if the ASV table has been calibrated or
#' just represent confounders that you wish to remove effects from in the changepoint regression.
#' \code{lm} is used to get residuals of the phenotype after these effects are regressed out.
#'
#' @keywords network, changepoint
#'
#' @import parallel
#' @import chngpt
#' @return A data frame of model results, same as \code{thresh}.
#' @examples
#'
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
#' netThresh_output <- netThresh(net_data,
#'   asvTab = asv, asvCols = NULL,
#'   clusterCol = "component_cluster", cluster = NULL, phenoCols = "biomass_z",
#'   model = "hinge", calibratePheno = "tissue"
#' )
#'
#' @export

netThresh <- function(net, asvTab, asvCols = NULL, clusterCol = NULL, cluster = NULL, phenoCols = NULL,
                      model = "hinge", cores = getOption("mc.cores", 1), calibratePheno = NULL) {
  #* `calculated values`
  nodes <- net[["nodes"]]
  if (is.null(cluster)) {
    cluster <- unique(net[["nodes"]][[clusterCol]])
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
