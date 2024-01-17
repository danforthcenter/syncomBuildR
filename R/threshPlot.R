#' Function to plot results of a changepoint model with data returned from \code{thresh}
#' 
#' 
#' @param thresh Output from \code{\link{thresh}}
#' @param asv The asv table used in making \code{thresh} object.
#' @param asvCols A vector of ASV column names. Defaults to NULL in which case all columns containing "ASV" are used and a list of ggplots is returned.
#' @param phenotype A vector of phenotype names in \code{thresh}. Defaults to NULL where all phenotypes are used and a list of plots is returned per ASV.
#' 
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
#' 
#' @export

threshPlot<-function(thresh, asv, asvCols=NULL, phenotype=NULL){
  if(is.null(asvCols)){ asvCols<-colnames(asv)[grepl("ASV", colnames(asv))] }
  if(is.null(phenotype)){phenotype<-unique(thresh$phenotype)}
  if("calibratePheno" %in% colnames(thresh)){
    cal<-unlist(lapply(unique(thresh$calibratePheno), function(s) {as.character(as.formula(s))[3]}))
  } else{cal=NULL}
  outList<-lapply(asvCols, function(microbe){
    thresh_sub<-thresh[thresh$asv == microbe, ]
    asv_sub<-asv[, colnames(asv) %in% c(phenotype, microbe, cal)]
    phenoPlots<-lapply(unique(thresh_sub$phenotype), function(pheno){
      #print(c(pheno, microbe))
      if(!is.null(cal)){
        asv_sub[[pheno]]<-residuals(lm(as.formula(thresh_sub[1,"calibratePheno"]), data=asv_sub))
      }
      interceptData<-thresh_sub[thresh_sub$Source =="(Intercept)" & thresh_sub$phenotype==pheno, ]
      chngptData<-thresh_sub[thresh_sub$Source != "(Intercept)" & thresh_sub$phenotype==pheno, ]
      postCptCol<-if(chngptData$p.value < 0.05){viridis::plasma(1, begin=0.7)}else{"black"}
      
      p<-ggplot2::ggplot(chngptData)+
        ggplot2::geom_point(data=asv_sub[asv_sub[[microbe]] < chngptData$changePoint,],
                            ggplot2::aes(.data[[microbe]], .data[[pheno]]),
                   color="gray40",size=2, alpha=0.5)+
        ggplot2::geom_point(data=asv_sub[asv_sub[[microbe]] > chngptData$changePoint,],
                            ggplot2::aes(.data[[microbe]], .data[[pheno]]),
                   color=postCptCol,size=2, alpha=0.85)+
        geom_vline(aes(xintercept = changePoint),linetype=5, color = "black")+
        geom_segment(data=chngptData, 
                     aes(x=min(asv_sub[[microbe]]),
                         xend=changePoint,
                         y=mean(interceptData$est),
                         yend=mean(interceptData$est)),
                     color="black")+
        geom_segment(data=chngptData, x=chngptData$changePoint,
                     xend=max(asv_sub[[microbe]]),
                     y=interceptData$est,
                     yend=interceptData$est + chngptData$est*(max(asv_sub[[microbe]])-chngptData$changePoint ),
                     color= postCptCol )+
        
       ylab(paste0(ifelse(is.null(cal), "","Calibrated "), pheno))+
        xlab(microbe)+
        ggtitle(paste0("ASV ",sub("ASV","",chngptData[1,"asv"])))+
        labs(subtitle = paste0("P-value: ", round(chngptData$p.value, 3)))+
        theme_light()+
        theme(axis.text = element_text(size = 8),
             axis.title= element_text(size = 9))
      return(p)
    })
  })
  
  if(length(asvCols)==1){
    outList<-outList[[1]]
  }
  if(length(phenotype)==1 & length(asvCols)==1){
    outList<-outList[[1]]
  }
  return(outList)
}
