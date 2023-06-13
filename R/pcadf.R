#' Function to run a PCA, plot and return the data with PC coordinates
#' 
#' @param df Dataframe to ordinate
#' @param cols columns to reduce dimensions of. Can be specified with names or positions. Defaults to all column names containing "ASV".
#' @param color column name used to color points in the pca plot.
#' @param returnData Logical, should data be returned? Defaults to TRUE where data and a ggplot are returned.
#' @param ncp Optional, number of principal components to return attached to dataframe if data is returned. Defaults to all.
#' @keywords pca
#' 
#' @import ggplot2
#' @import FactoMinerR
#' 
#' @examples
#' print(load("/home/jsumner/Desktop/stargate/SINC/sincUtils/syncomBuilder/cal_output.rdata"))
#' asv<-are_c[[1]]
#' zinbCalibrated = are_c[[2]][are_c[[2]]$model=="ZINB", "asv"]
#' asv<-are_c[[1]][, c("tissue","plot","row","genotype", "biomass","sd", zinbCalibrated)]
#' x<-pcadf(df=asv, cols=NULL, color=c("tissue", "genotype"), returnData=T, ncp=NULL)
#' x$plot
#' 
#' @export

pcadf<-function(df=NULL, cols=NULL, color=NULL, returnData=T, ncp=NULL){
  # df=asv; cols=NULL; color=c("tissue", "genotype"); returnData=T; ncp=NULL
  if(is.null(cols)){cols=which(grepl("ASV",colnames(df)))
  }else if(is.character(cols) & length(cols)==1){
    cols = which(grepl(cols,colnames(df)))
  } else if(!is.numeric(cols)){
    cols<-which(colnames(df) %in% cols) }
  if(is.null(ncp)){ncp=min(dim(df[,cols]))-1}
  pca<-FactoMineR::PCA(df[,cols], ncp=ncp, graph=F)
  pc1Var<-round(pca$eig[1,2], 3)
  pc2Var<-round(pca$eig[2,2], 3)
  coords<-as.data.frame(pca$ind)
  coords<-coords[,grepl("coord", colnames(coords))]
  colnames(coords)<-gsub("coord.Dim.", "pc", colnames(coords))
  pca.df<-cbind(df[, -cols], coords)
  if(length(color)>1){
    pca.df[[paste0(color, collapse=".")]] = interaction(pca.df[,c(color)])
    color=paste0(color, collapse=".")
  }
  if(is.null(color)){pca.df$dummyVariableForColor = 1; color="dummyVariableForColor"}
  p<-ggplot2::ggplot(pca.df, ggplot2::aes(x=pc1, y=pc2, color = .data[[color]]))+
    ggplot2::geom_point(alpha=0.85)+
    ggplot2::labs(x=paste0("PC 1 (",pc1Var,"%)"),y=paste0("PC 2 (",pc2Var,"%)"))+
    ggplot2::guides(color=ggplot2::guide_legend(override.aes=list(alpha=1)))+
    ggplot2::theme_minimal()+
    ggplot2::theme(legend.position="bottom",
          plot.title = ggplot2::element_text(size=14),
          axis.line.y.left = ggplot2::element_line(),
          axis.line.x.bottom = ggplot2::element_line())
  if(color=="dummyVariableForColor"){
   p<-p+ggplot2::theme(legend.position="none")
  }
  if(returnData){return(list("data"=pca.df, "plot"=p))}else{return(p)}
}


