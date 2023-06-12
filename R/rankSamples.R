#' Function to pick samples for limited dilution/other bench science following changepoint regression of an ASV table or network.
#' 
#' @param asvTab An asv table as returned by \code{\link{qc}} or \code{\link{cal}}.
#' @param thresh Output from \code{\link{thresh}} or \code{\link{netThresh}}.
#' @param network Optionally a network object as returned from \code{\link{asvNet}} and clustered with \code{\link{netClust}}. Defaults to NULL in which case the asv table is used and thresh is assumed to be changepoint regressions at the asv level.
#' @param id The name of the column containing sample names.
#' @param groups Optionally specify which ASVs/Clusters should be used. If NULL, the default, then significant models (p value below p.cutoff) are used to find groups.
#' @param phenotypes Optionally specify which phenotypes should be used. If NULL, the default, then all in thresh are used.
#' @param p.cutoff P value cutoff for statistical significance if groups is NULL. Defaults to 0.05.
#' 
#' @return A data frame of samples ranked by the impact of their nodes (ASVs) or clusters
#' 
#' @examples 
#' 
#' # a<-qc(); b<-cal(a); c<-thresh(b); d<-asvDist(a) ; e<-net(d, thresh = c) ;  f<-netClust(e) ; g <- rankSamples(f)
#' 
#' print(load("/home/jsumner/Desktop/stargate/SINC/sincUtils/syncomBuilder/net_output_clustered.rdata"))
#' print(load("/home/jsumner/Desktop/stargate/SINC/sincUtils/syncomBuilder/cal_output.rdata"))
#' print(load("/home/jsumner/Desktop/stargate/SINC/sincUtils/syncomBuilder/threshOutput.rdata"))
#' 
#' asv<-are_c[[1]]
#' zinbCalibrated = are_c[[2]][are_c[[2]]$model=="ZINB", "asv"]
#' asv<-are_c[[1]][, c("tissue","plot","row","genotype", "biomass","sd", zinbCalibrated)]
#' 
#' netThresh_output<-netThresh(net_data, asvTab=asv, asvCols=NULL, 
#'    clusterCol="kmeans_cluster", cluster=NULL, phenoCols="biomass",
#'    model="hinge", calibratePheno="genotype")
#' 
#' asv$sample = paste0(asv$plot, asv$tissue)
#' 
#' network_ranking<-rankSamples(asvTab=asv, thresh=netThresh_output, network=net_data, id="sample", groups=NULL, phenotypes=NULL, p.cutoff=0.85)
#' head(network_ranking)
#' asvTable_ranking<-rankSamples(asvTab=asv, thresh=threshMods, network=NULL, id="sample", groups=NULL, phenotypes=NULL, p.cutoff=0.05)
#' head(asvTable_ranking)
#' 
#' @export

rankSamples<-function(asvTab, thresh, network=NULL, id, groups=NULL, phenotypes=NULL, p.cutoff=0.05){
  if(!is.null(network)){network_mode = T} else {network_mode = F}
  if(network_mode){thresh_group = "cluster"}else{thresh_group = "asv"}
  
  samples = asvTab[[id]]
  starting_ncol = ncol(asvTab)+1
  if(is.null(phenotypes)){
    phenotypes = unique(thresh$phenotype)
  }
  
  if("calibratePheno" %in% colnames(thresh)){
    cal<-unlist(lapply(unique(thresh$calibratePheno), function(s) {as.character(as.formula(s))[3]}))
  } else{cal=NULL}
  
  if(!is.null(cal)){
    asvTab[, phenotypes]<-lapply(phenotypes, function(pheno){
      residuals(lm(as.formula( thresh[thresh$phenotype==pheno,"calibratePheno"][1] ), data=asvTab))
    })
  }
  if(is.null(groups)){
    groups <- thresh[thresh$phenotype %in% phenotypes & thresh$Source != "(Intercept)" & thresh$p.value < p.cutoff, thresh_group]
  }
  if(network_mode){
    clusterCol = thresh$clusterType[1]
    nodes<-network$nodes
    asvTab[, paste0(groups, "_sum")]<-lapply(unique(sub("cluster_", "", groups)), function(clust){
      asvs_in_clust<-nodes[nodes[[clusterCol]]==clust, "asv"]
      rowSums(as.matrix(asvTab[, asvs_in_clust]))
    })
    asvTab[, paste0(groups, "_rank")]<-lapply(unique(paste0(groups, "_sum")), function(clust){
      -1*(rank(asvTab[[clust]])-nrow(asvTab)-1)
    })
  }
  asvTab[, paste0(phenotypes, "_rank")]<-lapply(phenotypes, function(pheno){
    -1*(rank(asvTab[[pheno]])-nrow(asvTab)-1)
  })
  
  asvTab$mean_pheno_rank <- rowMeans(as.matrix(asvTab[, paste0(phenotypes, "_rank")]))
  if(network_mode){
    out<-cbind(setNames(data.frame(asvTab[[id]]), id), asvTab[, starting_ncol:ncol(asvTab)])
    out$mean_rank<-rowMeans(out[, c(paste0(groups, "_rank"), "mean_pheno_rank")])
  } else{
    out<-cbind(setNames(data.frame(asvTab[[id]], asvTab[,groups]), c(id, groups)), asvTab[, starting_ncol:ncol(asvTab)])
    out$mean_rank<-rowMeans(out[, c(groups, "mean_pheno_rank")])
  }
  out<-out[sort(out$mean_rank, decreasing=F, index.return=T)$ix, ]
  out$overall_rank<-1:nrow(out)
  return(out)
}
