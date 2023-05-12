#' Function to make distance/dissimilarity data from an ASV table.
#' 
#' I should inform this function using emd from pcvr.
#' That pipeline is probably really similar to what I want to have happen here.
#' In that case this should be ready to take either a dist object or a dataframe or a matrix similar to dist object
#' That way I can support all of vegdist.
#' 
#' So I need two functions here, asvDist and net
#' 
#' 
#' 
#' 
#' @param data 
#' @keywords changepoint, threshold, regression, phenotype
#' @import chngpt
#' @return A dataframe summarizing changepoint models for individual ASVs vs phenotypes.
#' 
#' @examples 
#' 
#' # a<-qc(); b<-cal(a); c<-thresh(b); d<-asvDist(a)
#' print(load("/home/jsumner/Desktop/stargate/SINC/sincUtils/syncomBuilder/threshOutput.rdata"))
#' 
#' 
#' 
#' @export

asvDist<-function(){}

#* [args]
asvTab = asv
asvCols = colnames(asvTab)[grepl("ASV", colnames(asvTab))]
method = "euclidean" # distance or correlation type or function?
parallel=getOption("mc.cores",1)
transform = "clr" # center log transform or bust? Are there other options I'd want?
plot=T
#des = c("tissue", "plot", "genotype")

#* `calculated values`

if(parallel > 1){innerLapply <- function(...){parallel::mclapply(..., mc.cores=parallel)}}else{innerLapply<-lapply}

vegan_distances<-c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", 
                   "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis",
                   "chisq", "chord", "aitchison", "robust.aitchison")
correlation_coefficients<-c("spearman", "pearson")

#* `make sim matrix`
method="spearman"
method="euclidean"

if(method %in% vegan_distances){
  mat<-as.matrix(t(asvTab[,asvCols]))
  mat[rowSums(mat)>0, colSums(mat)>0] #* `check that matrix is not too sparse, colsum 0 would kill?`
  dist_mat<-as.matrix(vegan::vegdist(mat, method=method))
  sim_mat<-1/dist_mat
  # if(!is.null(des)){
  #   colnames(sim_mat)<-rownames(sim_mat)<-as.character(interaction(asvTab[,des]))
  # }
} else if (method %in% correlation_coefficients){
  mat<-as.matrix(asvTab[,asvCols])
  mat[rowSums(mat)>0, colSums(mat)>0]
  M<-Hmisc::rcorr(mat, type=method) # review structure 
  rhos<-M[["r"]]
  sim_mat <- (1 - rhos) / 2 # turn correlation into a distance
} else if(is.function(method)){
  sim_mat<-method(asvTab[,asvCols], ...)
}

sim_mat[1:10,1:10]
#* control for NAs and INFs in any of the matrices
sim_mat[is.na(sim_mat) | is.infinite(sim_mat)]<-NA










rcors <- function(df, method="spearman", transform="clr", colPattern="ASV", list=NULL) { # Correlation based
  if(is.numeric(df) && length(df)==1){
    i=df
    df<-list[[i]]
    group<-names(list)[i]
  } else{
    group<-"grp"
  }
  df<-df[,grepl(colPattern, colnames(df))]
  mat<-as.matrix(df)
  mat<-mat[, colSums(mat)>0]
  if(transform=="clr"){df<-compositions::clr(df)}
  if(method %in% c("spearman", "pearson")){
    mat<-as.matrix(df)
    mat<-mat[, colSums(mat)>0]
    M <- Hmisc::rcorr(mat,type=method)
    rhos<-M[["r"]] # get correlation coeffients
    dd <- stats::as.dist((1 - rhos) / 2)  # rescale so Rho = 1 is 0 for hclust-ing
  } else {stop("method must be one of spearman or pearson to make a network by correlations")}
  hc <- stats::hclust(dd, method = "complete")
  K=floor(length(hc$order)/15)
  x<-rhos[hc$order, hc$order]
  clusts<-cutree(hc,k=K)
  cList<-lapply(unique(clusts), function(i){
    nm<-names(clusts[clusts==i])
    x[nm,nm]
  })
  cList<-cList[order(unlist(lapply(cList,mean)),decreasing=T)]
  colOrder<-unlist(lapply(cList,colnames))
  out <- do.call(rbind,lapply(1:length(M), function(m) { 
    x<-as.data.frame(M[[m]])
    x$rowname<-rownames(x)
    x$trait<-names(M)[[m]]
    l<-tidyr::pivot_longer(x, cols = contains("ASV"))
    as.data.frame(l)
  }))
  out<-as.data.frame(tidyr::pivot_wider(out, names_from=trait, values_from=value))
  colnames(out)<-c("asv1", "asv2", names(M))
  head(out)
  out$p_sig<-ifelse(out$P<=0.05 & !is.na(out$P), out$P, NA)
  out$rho_strong<-ifelse(abs(out$r)>0.5 & !is.na(out$r) & !is.na(out$p_sig), out$r, NA)
  out$group<-group
  if(any(!unique(out$asv2)%in%colOrder) | any(!unique(out$asv1)%in%colOrder)){
    addAsvs<-unique( unique(out$asv1)[which(!unique(out$asv1)%in%colOrder)],
                     unique(out$asv1)[which(!unique(out$asv1)%in%colOrder)] )
    colOrder<-c(colOrder,addAsvs)
  }
  out$asv1 <- factor(out$asv1,levels=colOrder)
  out$asv2 <- factor(out$asv2,levels=colOrder)
  return(out)
}




