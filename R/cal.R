#' Function to calibrate asv counts for some set of confounder variables
#'
#' @param asvTab ASV table similar to that returned by \code{\link{qc}}, often with additional metadata joined.
#' @param asvCols Numeric index of columns representing ASVs in asvTab. Defaults to all columns containing "ASV" in their column name. Note that missing values in these columns will be treated as 0s.
#' @param cal Vector of columns to use to calibrate asv counts. These can be categorical or continuous and are passed to MASS::glm.nb or pscl::zeroinfl, depending on the presense of 0s in the data.
#' @param ZI_cutoff A proportion of 0s that is considered too many. Higher values will fit models to more of the data, but those models may be worse fits for calibration. This defaults to 0.9 which is rather high.
#' @param cores Number of cores to use. Passed to parallel::mclapply. Defaults to 1 if the "mc.cores" option is not set.
#' @param verbose Logical, should fit metrics be returned? If TRUE then the output is a list with two elements (new ASV table and the metrics data), otherwise only the new ASV table is returned.
#' @param model Type of changepoint model in chngpt::chngptm labeling convention. See Figure 2.1 of the \link{[chngpt vignette](https://cran.r-project.org/web/packages/chngpt/vignettes/chngpt-vignette.pdf)}
#' @keywords calibrate, ZINB, NB
#' @import parallel
#' @import pscl
#' @import MASS
#' @return An ASV table as a wide dataframe with effects from some variables modeled out.
#'
#' @examples
#'
#' print(load("/home/jsumner/Desktop/stargate/SINC/sincUtils/syncomBuilder/qc_output.rdata"))
#' asv[1:10, 1:10]
#' print(load("/home/jsumner/Desktop/stargate/SINC/sincUtils/syncomBuilder/biomass.rdata"))
#' agBiomass <- aggregate(biomass ~ plot + row + genotype, biomass, mean)
#' agBiomass$sd <- aggregate(biomass ~ plot + row + genotype, biomass, sd)$biomass
#' asv_joined <- plyr::join(asv, agBiomass, by = "plot")
#' are_c <- cal(asv_joined[asv_joined$tissue == "ARE", ], cal = "genotype", cores = 10)
#' table(are_c[[2]]$model)
#' dim(are_c[[1]])
#'
#' @export


cal <- function(asvTab = NULL, asvCols = which(grepl("ASV", colnames(asvTab))),
                cal = NULL, ZI_cutoff = 0.9, cores = getOption("mc.cores", 1), verbose = TRUE) {
  if (is.null(cal)) {
    stop("Argument 'cal' is NULL but requires a value")
  }
  x <- parallel::mclapply(asvCols, function(i) {
    sub <- asvTab[, -asvCols]
    microbes <- asvTab[[i]]
    microbesNoNA <- ifelse(is.na(microbes), 0, microbes)
    sub$microbe <- microbesNoNA
    sub <- sub[complete.cases(sub[, cal]), ]
    if (mean(sub$microbe == 0) <= ZI_cutoff) {
      form <- as.formula(paste0("microbe ~ ", paste(cal, collapse = "+")))
      tryCatch(
        {
          if (all(sub$microbe > 0)) {
            mod1 <- MASS::glm.nb(data = sub, form)
            which_mod <- "NB"
          } else {
            mod1 <- pscl::zeroinfl(data = sub, form, dist = "negbin")
            which_mod <- "ZINB"
          }
          out <- sub
          out$ASVresid <- residuals(mod1)
          if (which_mod == "NB") {
            out$ASVnew <- out$ASVresid * exp(mod1$coefficients[1])
          } else {
            out$ASVnew <- out$ASVresid * exp(mod1$coefficients$count[1])
          }
          mod2 <- lm(data = out, microbe ~ ASVnew)
          out$ASVnew <- out$ASVnew * mod2$coefficients[2] + mod2$coefficients[1]
          out$ASVnew <- out$ASVnew + abs(min(out$ASVnew))
          out$nonZero = mean(sub$microbe == 0)
          out$model <- which_mod
          out$asvName <- colnames(asvTab)[i]
          out
        },
        warning = function(war) {},
        error = function(err) {}
      )
    } else {
      out <- sub
      out$ASVresid <- NA
      out$ASVnew <- sub$microbe
      out$nonZero = mean(sub$microbe == 0)
      out$model <- "none"
      out$asvName <- colnames(asvTab)[i]
      out
    }
  }, mc.cores = cores)
  names(x) <- colnames(asvTab)[asvCols]
  modelErrors <- which(unlist(lapply(x, is.null)))
  if (length(modelErrors) > 0) {
    modelErrorsCols <- colnames(asvTab)[asvCols][modelErrors]
    fitCols <- colnames(asvTab)[asvCols][-modelErrors]
    dfs <- x[-modelErrors]
  } else {
    modelErrorsCols <- character(0)
    fitCols <- colnames(asvTab)[asvCols]
    dfs <- x
  }

  first <- dfs[[1]][, c(colnames(asvTab[, -asvCols]), 'ASVnew')]
  colnames(first)[ncol(first)] <- dfs[[1]][1, "asvName"]
  newAsvTab <- cbind(first, do.call(cbind, lapply(2:length(dfs), function(d) {
    setNames(data.frame(dfs[[d]]$ASVnew), dfs[[d]][1, "asvName"])
  })))
  if (verbose) {
    metrics <- rbind(
      do.call(rbind, lapply(dfs, function(d) {
        data.frame(asv = d[1, "asvName"], model = d[1, "model"], nonZeroPct = d[1, 'nonZero'])
      })),
      data.frame(
        asv = modelErrorsCols, model = rep("error", length(modelErrorsCols)),
        nonZeroPct = unlist(lapply(modelErrorsCols, function(e) mean(asvTab[[e]] == 0)))
      )
    )
    outList <- list("asvTab" = newAsvTab, "details" = metrics)
  } else {
    outList <- newAsvTab
  }
  return(outList)
}
