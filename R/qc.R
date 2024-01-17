#' Ease of use quality control function for DADA2 output.
#' 
#' @param file A file path to an Rdata object containing DADA2 output in such as that from DADA2 examples (function?)
#' @param asvTab ASV table from DADA2. If a file path is given then this defaults to seqtab.print if left NULL for consistency with common DADA2 examples and previous SINC work.
#' @param taxa A taxonomy file for the ASVs in asvTab. This defaults to taxa_rdp if left NULL.
#' @param asvAbnd Lower bound on sum per microbe across all samples, ASVs with fewer total counts than this number will be removed.
#' @param sampleAbnd Lower bound on sum of all microbes per sample, samples with fewer total counts than this number will be removed. Note this filtering happens before `asvAbnd` filtering.
#' @param normalize Method for normalization, currently supported options are NULL (where no normalization is done),
#' "rescale" (the default) where samples are rescaled to have the same sample abundance using the exponent of the mean log sample sum,
#'  and "pqn" where probabilistic quotient normalization is used.
#' @param rmTx A list of taxa filters. Each element in the list should specify a level of the taxonomy, an affector, and a comma separated string of values to filter that taxonomy for. See details.
#' @param separate An optional dataframe or vector of the same length as the number of samples in asvTab. If provided then this separates the asvTab by these groupings and abundance filtering will be done separately in each subset of the data. Subsets will be bound together and missing values will be filled with NA should subsets contain different selections of ASVs.
#' @param metadata An optional dataframe of metadata to add to the asv table. If left NULL while `separate` is a dataframe then this will include that data as metadata.
#' @keywords DADA2
#' @importFrom plyr rbind.fill
#' @return An ASV table as a wide dataframe
#' @details 
#' 
#' The `rmTx` argument specifies a list of filters, with one element per each level of taxonomy that should be filtered. Each element in the list should take the form of "Taxalevel in example, example2".
#' Which will remove ASVs where their "Taxalevel" is "example" or "example2". Alternatively the affector can be "contains" in which case regex is used.
#' In that case the string "Class contains plant, other" would search for matches to "plant|other" in the Class column of the taxonomy. More general regex patterns are supported. 
#' 
#'
#' @examples 
#' 
#' print(load("~/Desktop/stargate/SINC/sincUtils/syncomBuilder/field_2021_small_ex.rdata"))
#' separate<-data.frame(tissue = stringr::str_extract(samps, "[:alpha:]+"))
#' metadata = cbind(separate, data.frame(sample = samps,
#'   plot = stringr::str_extract(samps, "[0-9]+")))
#' 
#' file=NULL; asvTab = seqtab.print; taxa=NULL; asvAbnd = 100; sampleAbnd=1000; normalize="pqn"
#' rmTx=list("Order in Chloroplast", "Phylum in Plantae"); separate=separate
#' metadata=cbind(separate, data.frame(plot = stringr::str_extract(samps, "[0-9]+")))
#' 
#' asv<-qc(asvTab = seqtab.print, taxa=taxa_rdp, asvAbnd = 100, sampleAbnd=1000, rescale=T, 
#' rmTx=list("Order in Chloroplast", "Phylum in Plantae"), separate=NULL)
#' asv<-qc(asvTab = seqtab.print, taxa=taxa_rdp, asvAbnd = 100, sampleAbnd=1000, rescale=T,
#' rmTx=list("Order in Chloroplast", "Phylum in Plantae"), separate=separate, metadata=metadata)
#' save(asv, file="../qc_output.rdata")
#' 
#' @export

qc<-function(file=NULL, asvTab = NULL, taxa=NULL, asvAbnd = 100, sampleAbnd=1000, normalize="rescale",
             rmTx=list("Order in Chloroplast", "Phylum in Plantae"), separate=NULL, metadata=NULL){
  if(!is.null(file)){load(file)}
  if(is.null(taxa)){taxa<-taxa_rdp}
  if(is.null(asvTab)){asvTab<-seqtab.print}
  if(!is.data.frame(asvTab)){asvTab<-as.data.frame(asvTab)}
  if(!is.data.frame(taxa)){taxa<-as.data.frame(taxa)}
  split<-if(!is.null(separate)){T}else{F}
  if(is.null(metadata) & is.data.frame(separate)){metadata = separate}
  
  if(!is.null(rmTx)){
    indsL<-lapply(rmTx, function(stri){
      taxaLevel = strsplit(stri," ")[[1]][1]
      affector<-strsplit(stri," ")[[1]][2]
      values<-trimws(gsub(",", "", strsplit(stri," ")[[1]][-c(1:2)]))
      
      if(affector %in% c("in", "is", "=")){
        indL<-lapply(values, function(value){
          index<-taxa[[taxaLevel]]!=value # flag F if the row in taxa should be removed
          index[is.na(index)]<-T # if taxaLevel[i] is NA then keep it
          index
        })
        ind<-apply(matrix(unlist(indL), ncol=length(indL) ), MARGIN=1, FUN=all) # if all are T then keep the row
      } else if(affector=="contains"){
        valReg<-paste0(values, collapse="|")
        ind<-!grepl(valReg, taxa[[taxaLevel]])
      }
      ind
    })
    keep_index<-apply(matrix(unlist(indsL), ncol=length(indsL) ), MARGIN=1, FUN=all) 
    taxa<-taxa[keep_index,]
    asvTab<-asvTab[,keep_index]
  }
  
  if(!is.null(sampleAbnd)){
    originalRowNames<-rownames(asvTab)
    if(!split){
    asvTab<-asvTab[rowSums(asvTab) >= sampleAbnd, ]
    }else{
      asvTab<-do.call(rbind, 
              lapply(split(asvTab, interaction(separate, drop=T) ), function(d){
                d[rowSums(d) >= sampleAbnd, ] }))
    rownames(asvTab)<-sub( paste0(levels(interaction(separate, drop=T)),".", collapse="|"),"",rownames(asvTab))
    sepColNames<-colnames(separate)
    separate<-as.data.frame(separate[which(rownames(asvTab) %in% originalRowNames),])
    colnames(separate)<-sepColNames
    }
    if(!is.null(metadata)){
      metadata<-metadata[which(rownames(asvTab) %in% originalRowNames),]
    }
    
  }
  if(!is.null(asvAbnd)){
    if(!split){asvTab<-asvTab[, colSums(asvTab)>=asvAbnd]
    }else{
      asvTab<-do.call(plyr::rbind.fill, 
                      lapply(split(asvTab, interaction(separate, drop=T)), function(d){
                        d[,colSums(d)>=asvAbnd] }))
    }
  }
  if(!is.null(normalize)){
    if(normalize == "rescale"){
      scaleFactor<-exp(mean(log(rowSums(asvTab, na.rm=T))))
      asvTab<-as.data.frame(t(apply(asvTab ,MARGIN=1, FUN = function(i) round(scaleFactor*i/sum(i, na.rm=T)) )))
    } else if (normalize == "pqn"){
      median_vec <- apply(asvTab / rowMeans(asvTab, na.rm = TRUE), 2, median, na.rm = TRUE)+1
      normalized_asvTab <- t(apply(asvTab, 1, function(x) x/median_vec))
      rownames(normalized_asvTab) <- rownames(asvTab)
      asvTab <- normalized_asvTab
    } else{
      warning(paste0("Available options for normalization are NULL, 'pqn', and 'rescale', ", normalize, " is not implemented"))
    }
  }
  
  if(!is.null(metadata)){
    asvTab<-cbind(metadata, asvTab)
  }
  return(asvTab)
}
