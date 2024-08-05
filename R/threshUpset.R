#' Function to make upset plots of thresh output
#'
#' @param thresh output from the thresh function
#' @param cores Number of cores to use in summarizing data. This defaults to 1 if the "mc.cores" option
#' is not set.
#' @param method The p-adjustment method to be used. Can be any suffix from ?p.adjust.methods
#' @keywords changepoint, threshold, regression, phenotype
#' @import ComplexUpset
#' @import patchwork
#' @importFrom parallel mclapply
#' @importFrom ggplot2 labs
#' @importFrom data.table as.data.table dcast
#' @return A list of ggplots
#' 
#' @examples 
#' 
#' # a<-qc(); b<-cal(a); c<-thresh(b)
#' print(load("/home/jsumner/Desktop/stargate/SINC/sincUtils/syncomBuilder/cal_output.rdata"))
#' asv<-are_c[[1]]
#' zinbCalibrated = are_c[[2]][are_c[[2]]$model=="ZINB", "asv"]
#' asv<-are_c[[1]][, c("tissue","plot","row","genotype", "biomass","sd", zinbCalibrated)]
#' 
#' threshMods<-thresh(asvTab=asv, phenoCols="biomass", model="hinge", cores=10, calibratePheno = "genotype")
#' threshUpset(threshMods, 1)
#' 
#' @export

threshUpset <- function(thresh, cores = getOption("mc.cores", 1), method = "fdr") {
  df <- do.call(rbind, lapply(split(thresh,
                                    interaction(thresh$phenotype)), function(d){
                                      d$p.adj.bonferroni <- p.adjust(d$p.value., method = "bonferroni")
                                      d$p.adj.holm <- p.adjust(d$p.value., method = "holm")
                                      d$p.adj.hochberg <- p.adjust(d$p.value., method = "hochberg")
                                      d$p.adj.fdr <- p.adjust(d$p.value., method = "fdr")
                                      d$p.adj.by <- p.adjust(d$p.value., method = "BY")
                                      d$p.adj.none <- p.adjust(d$p.value., method = "none")
                                      return(d)
                                    }))

  adj_results <- as.data.frame(do.call(cbind, lapply(df[, grepl("p.adj", colnames(df))], function(col){
    as.numeric(col<0.05)
  })))
  p0 <- ComplexUpset::upset(adj_results, intersect = colnames(adj_results))
  p0[[2]] <- p0[[2]] +
    ggplot2::labs(title = paste0("P adjustment options in ", tiss))
  p.adj.col <- paste0("p.adj.", method)
  sig <- do.call(rbind, parallel::mclapply(unique(df$asv), function(asv){
    sub <- df[df$asv == asv, ]
    if(any(sub[[p.adj.col]]<0.05)){
      return(sub)
    } else{NULL}
  }, mc.cores=cores))

  sub_upsetData<-sig[, c("phenotype", p.adj.col, "asv")]
  sub_upsetData[[p.adj.col]] <- ifelse(sub_upsetData[[p.adj.col]] < 0.05, 1, 0)
  dt_l <- data.table::as.data.table(sub_upsetData)

  dw <- as.data.frame(data.table::dcast(dt_l, asv ~ phenotype, value.var = p.adj.col))

  p01 <- ComplexUpset::upset(dw, intersect = colnames(sub_upsetData)[-1])
  p01[[2]] <- p01[[2]] + ggplot2::labs(title = paste0("ASV ~ Phenotype Correlations in ", tiss),
                             subtitle=paste0(method, " P Values Intersections"))

  return(list("p_adj_plot" = p0, "plot" = p01))
}
