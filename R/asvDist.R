#' Function to make distance/dissimilarity data from an ASV table.
#' 
#' 
#' @param asvTab ASV table to use. Should have rows as samples and columns representing ASVs.
#' @param asvCols Either column names, positions, or a single character string to use as regex to find all asv columns.
#' @param method A distance method compatible with \code{\link{vegan::vegdist}}, a correlation coefficient (spearman or pearson), or a custom function that returns a matrix of pairwise distances.
#' @param parallel How many cores to use in parallel. The parallelized component only comes into play with correlation coefficient methods and is not particularly heavy, so more cores will have a minimal speed improvement.
#' @param clr_transform Should ASV counts be transformed with a centered log ratio? This defaults to TRUE to alleviate the compositional nature of the data.
#' @param edgeFilter A filter for the edges returned. This defaults to NULL where no filtering is done. Optionally this can be numeric, in which case distances/correlations below that number are removed, or it can be a number as a string in which case that quantile is used (ie, "0.5" filters for above the median.).
#' @param plot Logical, should the data be plotted? Defaults to FALSE.
#' @param returnASVs Should ASVs or samples be considered the experimental unit? Defaults to TRUE in which case distances/correlations between ASVs across samples are returned. If FALSE then samples are compared by their ASV composition. Note that different methods are appropriate depending on whether samples or ASVs are being considered "nodes".
#' @keywords changepoint, threshold, regression, phenotype
#' @import chngpt
#' @import vegan
#' @import parallel
#' @import data.table
#' @import ggplot2
#' @return A dataframe showing pairwise correlations between individual ASVs/samples.
#' 
#' @examples 
#' 
#' # a<-qc(); b<-cal(a); c<-thresh(b); d<-asvDist(a)
#' print(load("/home/jsumner/Desktop/stargate/SINC/sincUtils/syncomBuilder/threshOutput.rdata"))
#' sp_dist<-asvDist(asv, method="spearman", clr_transform=TRUE, edgeFilter=0.5)
#' unfiltered_sp_dist<-asvDist(asv, method="spearman", clr_transform=TRUE, edgeFilter=NULL)
#' dim(sp_dist)
#' euc<-asvDist(asv, method="euclidean", clr_transform=FALSE, edgeFilter=NULL)
#' bray<-asvDist(asv, method="bray", clr_transform=FALSE, edgeFilter=NULL)
#' save(sp_dist, euc, bray, file="/home/jsumner/Desktop/stargate/SINC/sincUtils/syncomBuilder/asvDist_output.rdata")
#' 
#' 
#' @export

asvDist<-function(asvTab, asvCols = NULL, method="spearman",
                  parallel=getOption("mc.cores", 1), clr_transform = TRUE, 
                  edgeFilter = NULL, plot=FALSE, returnASVs =TRUE){
    #* `calculated values`
    
    if(is.null(asvCols)){asvCols=colnames(asvTab)[grepl("ASV", colnames(asvTab))]}
    #if(parallel > 1){innerLapply <- function(...){parallel::mclapply(..., mc.cores=parallel)}}else{innerLapply<-lapply}
    vegan_distances<-c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", 
                       "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis",
                       "chisq", "chord", "aitchison", "robust.aitchison")
    correlation_coefficients<-c("spearman", "pearson")
    
    #* *pull asvs as matrix* [same for all options]
    if(returnASVs){
    mat<-as.matrix(t(asvTab[,asvCols]))
    } else{mat<-as.matrix(asvTab[,asvCols])
    }
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
      M<-Hmisc::rcorr(t(mat), type=method)
      M[[method]] <- (1 - M[["r"]]) / 2 # turn correlation into a distance
      ldf <- do.call(rbind,parallel::mclapply(1:length(M), function(m) { 
        x<-as.data.frame(M[[m]])
        x$rowname<-rownames(x)
        x$trait<-names(M)[[m]]
        data.table::melt(data.table::as.data.table(x), id.vars = c("rowname", "trait"))
      }, mc.cores=parallel))
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











