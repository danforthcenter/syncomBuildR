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
#' sp_dist<-asvDist(asv, method="spearman", clr_transform=T, edgeFilter=0.5)
#' unfiltered_sp_dist<-asvDist(asv, method="spearman", clr_transform=T, edgeFilter=NULL)
#' dim(sp_dist)
#' euc<-asvDist(asv, method="euclidean", clr_transform=F, edgeFilter=NULL)
#' bray<-asvDist(asv, method="bray", clr_transform=F, edgeFilter=NULL)
#' 
#' @export

asvDist<-function(asvTab, asvCols = NULL, method="spearman",
                  parallel=getOption("mc.cores", 1), clr_transform = T, 
                  edgeFilter = NULL, plot=F, outputname = NULL){
    #* `calculated values`
    
    if(is.null(asvCols)){asvCols=colnames(asvTab)[grepl("ASV", colnames(asvTab))]}
    if(parallel > 1){innerLapply <- function(...){parallel::mclapply(..., mc.cores=parallel)}}else{innerLapply<-lapply}
    vegan_distances<-c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", 
                       "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis",
                       "chisq", "chord", "aitchison", "robust.aitchison")
    correlation_coefficients<-c("spearman", "pearson")
    
    #* *pull asvs as matrix* [same for all options]
    
    mat<-as.matrix(t(asvTab[,asvCols]))
    old_dim<-dim(mat)
    mat<-mat[rowSums(mat)>0, colSums(mat)>0]
    new_dim<-dim(mat)
    if(any(old_dim!=new_dim)){warning("Rows and Columns that summed to 0 have been removed.")}
    
    #* *apply transformations* [CLR makes most sense, but there could be others in the future]
    
    if(clr_transform){mat<-compositions::clr(mat)}
    
    #* *make long data* [from similarity matrices]
    
    if(method %in% vegan_distances){
      dist_mat<-as.matrix(vegan::vegdist(mat, method=method))
      sim_mat<-as.data.frame(1/dist_mat)
      sim_mat$c1<-rownames(sim_mat)
      ldf<-data.table::melt(data.table::as.data.table(sim_mat), id.vars = c("c1"), variable.name="c2", value.name=method)
      cor_method=F
      } else if (method %in% correlation_coefficients){
      M<-Hmisc::rcorr(mat, type=method)
      M[[method]] <- (1 - M[["r"]]) / 2 # turn correlation into a distance
      ldf <- do.call(rbind,lapply(1:length(M), function(m) { 
        x<-as.data.frame(M[[m]])
        x$rowname<-rownames(x)
        x$trait<-names(M)[[m]]
        data.table::melt(data.table::as.data.table(x), id.vars = c("rowname", "trait"))
      }))
      ldf<-as.data.frame(data.table::dcast(ldf, ... ~ trait))
      colnames(ldf)<-c("c1", "c2", names(M))
      cor_method=T
    } else if(is.function(method)){
      ldf<-method(asvTab[,asvCols], ...)
      cor_method=F
      method=newMethodNameFrom_ldf
    }
    
    #* *plotting*
    if(plot){
      # consider grouping things here maybe with hclust.
    p<-ggplot2::ggplot(ldf, ggplot2::aes(x=c1,y=c2,fill=ldf[[method]]))+
      ggplot2::geom_tile(color=NA,linewidth=0)+
      viridis::scale_fill_viridis(na.value = "grey100")+
      ggplot2::labs(fill=method)+
      ggplot2::theme(axis.title=ggplot2::element_blank(),
                     axis.text=ggplot2::element_blank(),
                     axis.line=ggplot2::element_blank())+
      ggplot2::theme_minimal()
    print(p)
    }
    
    #* *filter long data* [based on strength or P value or both]
    #* First remove NA/Inf values since those are bad for us
    
    ldf[[method]][is.na(ldf[[method]]) | is.infinite(ldf[[method]])]<-NA
    ldf<-ldf[!is.na(ldf[[method]]), ]
    if(!is.null(edgeFilter)){
      if(is.character(edgeFilter)){
        cutoff<-quantile(ldf[[method]], probs = as.numeric(edgeFilter))
        ldf<-ldf[ldf[[method]] >= as.numeric(cutoff), ]
      } else if(is.numeric(edgeFilter)){
        ldf<-ldf[ldf[[method]] >= edgeFilter, ]
      } else{stop("edgeFilter must be character or numeric, see ?asvDist for details.")}
    }
    #* *return data*
    return(ldf)
}











