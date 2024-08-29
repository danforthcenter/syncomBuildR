#' Function to plot results of a changepoint model with data returned from \code{thresh}
#'
#'
#' @param thresh Output from \code{\link{thresh}}
#' @param asv The asv table used in making \code{thresh} object.
#' @param asvCols A vector of ASV column names. Defaults to NULL in which case all columns containing
#' "ASV" are used and a list of ggplots is returned.
#' @param phenotype A vector of phenotype names in \code{thresh}. Defaults to NULL where all phenotypes
#' are used and a list of plots is returned per ASV.
#' @param unit The unit or scale of the changepoint models. This defaults to "asv" for use with thresh
#' and "cluster" should be used with netThresh output.
#' @param net The asvNet object if netThresh output is being plotted.
#'
#' @keywords changepoint, threshold, regression, phenotype, ggplot
#'
#' @import ggplot2
#' @import viridis
#'
#' @return A ggplot or list of ggplots showing changepoint models against some set of phenotypes.
#'
#'
#' @examples
#'
#' # a<-qc(); b<-cal(a); c<-thresh(b)
#' print(load("/home/jsumner/Desktop/stargate/SINC/sincUtils/syncomBuilder/threshOutput.rdata"))
#' head(threshMods)
#'
#' threshPlot(threshMods, asv, "ASV9")
#'
#' @export

threshPlot <- function(thresh, asv, asvCols = NULL, phenotype = NULL, unit = "asv", net = NULL) {
  if (is.null(asvCols)) {
    asvCols <- unique(thresh[[unit]])
  }
  if (!is.null(net)) {
    #* calculate a summarized ASV table by clusters
    nodes <- net[["nodes"]]
    clusterCol <- thresh$clusterType[1]
    clusters <- unique(nodes[[clusterCol]])
    clust_ag <- do.call(cbind, lapply(clusters, function(clust) {
      asvs_in_cluster <- nodes[nodes[[clusterCol]] == clust, "asv"]
      setNames(data.frame(rowSums(asv[, c(asvs_in_cluster)])), c(paste0("cluster_", clust)))
    }))
    asv <- cbind(asv[, -which(grepl("ASV", colnames(asv)))], clust_ag)
  }
  if (is.null(phenotype)) {
    phenotype <- unique(thresh$phenotype)
  }
  if ("calibratePheno" %in% colnames(thresh)) {
    cal <- unique(unlist(lapply(unique(thresh$calibratePheno), function(s) {
      trimws(strsplit(as.character(as.formula(s))[3], "[+|:|*]")[[1]])
    })))
  } else {
    cal <- NULL
  }
  outList <- lapply(asvCols, function(microbe) {
    thresh_sub <- thresh[thresh[[unit]] == microbe, ]
    phenoPlots <- lapply(phenotype, function(pheno) {
      asv_sub <- asv[, colnames(asv) %in% c(pheno, microbe, cal)]
      if (!is.null(cal)) {
        asv_sub[[pheno]] <- residuals(lm(
          as.formula(
            paste0(
              pheno, "~",
              gsub(".*~", "", thresh_sub[1, "calibratePheno"])
            )
          ),
          data = asv_sub, na.action = na.exclude
        ))
      }
      interceptData <- thresh_sub[thresh_sub$Source == "(Intercept)" & thresh_sub$phenotype == pheno, ]
      chngptData <- thresh_sub[thresh_sub$Source != "(Intercept)" & thresh_sub$phenotype == pheno, ]
      if (nrow(chngptData) < 1) {
        return(NULL)
      }
      postCptCol <- if (chngptData$p.value < 0.05) {
        viridis::plasma(1, begin = 0.7)
      } else {
        "black"
      }

      p <- ggplot2::ggplot(chngptData) +
        ggplot2::geom_point(
          data = asv_sub[asv_sub[[microbe]] <= chngptData$changePoint, ],
          ggplot2::aes(.data[[microbe]], .data[[pheno]]),
          color = "gray40", size = 2, alpha = 0.5
        ) +
        ggplot2::geom_point(
          data = asv_sub[asv_sub[[microbe]] > chngptData$changePoint, ],
          ggplot2::aes(.data[[microbe]], .data[[pheno]]),
          color = postCptCol, size = 2, alpha = 0.85
        ) +
        geom_vline(aes(xintercept = .data[["changePoint"]]), linetype = 5, color = "black") +
        geom_segment(
          data = chngptData,
          aes(
            x = min(asv_sub[[microbe]]),
            xend = .data[["changePoint"]],
            y = mean(interceptData$est),
            yend = mean(interceptData$est)
          ),
          color = "black"
        ) +
        geom_segment(
          data = chngptData, x = chngptData$changePoint,
          xend = max(asv_sub[[microbe]]),
          y = interceptData$est,
          yend = interceptData$est + chngptData$est * (max(asv_sub[[microbe]]) - chngptData$changePoint),
          color = postCptCol
        ) +
        ylab(paste0(ifelse(is.null(cal), "", "Calibrated "), pheno)) +
        xlab(microbe) +
        ggtitle(paste0(unit, sub(unit, "", chngptData[1, unit], ignore.case = TRUE))) +
        labs(subtitle = paste0("P-value: ", round(chngptData$p.value, 3))) +
        theme_light() +
        theme(
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 9)
        )
      return(p)
    })
    return(phenoPlots)
  })

  if (length(asvCols) == 1) {
    outList <- outList[[1]]
  }
  if (length(phenotype) == 1 && length(asvCols) == 1) {
    outList <- outList[[1]]
  }
  return(outList)
}
