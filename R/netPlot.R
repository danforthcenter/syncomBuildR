#' Function to plot networks generated from \code{asvNet} with various emphases.
#' 
#' 
#' @param net Object returned from \code{\link{asvNet}}
#' @param fill Optional column name to fill points by. Accepts single column names, vectors of p-value columns, or "thresh" which will match all p-value columns if the network was fit with threshold model data from \link{\code{thresh}}.
#' @param shape Optional column name to use for node shapes. Accepts the same options as fill. 
#' @param size Size for points. Passed to ggplot2::geom_point.
#' @param edgeWeight Optional weighting for edges. Must be present in the "edges" of net. Default of NULL will show equal size edges between all connected nodes.
#' @param edgeFilter Optional value to filter edges for. If non-NULL then only edges with edgeWeight greater than this value are kept. This can be a character vector or a numeric. Character vectors are interpreted as quantiles ("0.5" corresponds to the top 50% are kept).
#' @param thresh_below Significant cutoff if p-value columns are used for fill or shape. Defaults to 0.05.
#' @keywords network, changepoint
#' @import igraph
#' @import data.table
#' @return A named list with three elements:
#' "Nodes" is a dataframe of nodes and their metadata
#' "Edges" is a dataframe of edges connecting nodes.
#' "graph" is the igraph object used to generate the dataframes.
#' 
#' 
#' @examples 
#' 
#' # a<-qc(); b<-cal(a); c<-thresh(b); d<-asvDist(a) ; e<-net(d, thresh = c)
#' print(load("/home/jsumner/Desktop/stargate/SINC/sincUtils/syncomBuilder/net_output.rdata"))
#' 
#' net = net_data
#' fill=NULL
#' shape=NULL
#' size=3
#' edgeWeight="spearman"
#' edgeFilter=NULL
#' thresh_below = 0.05
#' net.plot(net, fill, shape, size, edgeWeight, edgeFilter, thresh_below)
#' net.plot(net, fill="thresh", shape, size, edgeWeight, edgeFilter, thresh_below)
#' 
#' 
#' @export
#' 

net.plot<-function(net, fill=NULL, shape=NULL, size = 3, edgeWeight=NULL, edgeFilter = NULL, thresh_below=0.05){
  
  nodes<-net[["nodes"]]
  edges<-net[["edges"]]
  multi_thresh_fill=F
  single_thresh_fill=F
  multi_thresh_shape=F
  single_thresh_shape=F
  #* make fill work
  if(is.null(fill)){fill="NOFILL"
    edges$NOFILL="a"
    nodes$NOFILL="a"
  } else if(fill=="thresh"){
    fill = colnames(nodes)[grepl("[hinge|upperhinge|segmented]_p", colnames(nodes))]
    multi_thresh_fill=T
  } else if (length(fill)>1){
    multi_thresh_fill = T
  } else if(is.numeric(nodes[[fill]]) & grepl("[hinge|upperhinge|segmented]_p$", fill)){
    single_thresh_fill=T
  }
  if(multi_thresh_fill | single_thresh_fill){
    nodes[[paste0(fill, "_bin")]]<-unlist(lapply(fill, function(col){ as.numeric(nodes[[col]]<=thresh_below) }))
    nodes[["significantThresholdModels"]]<-rowSums(as.data.frame(nodes[[paste0(fill, "_bin")]]))
    fill="significantThresholdModels"
  }
  #* make shape work
  if(is.null(shape)){
    shape = "NOSHAPE"
    nodes$NOSHAPE="a"
  } else if(shape=="thresh"){
    shape = colnames(nodes)[grepl("[hinge|upperhinge|segmented]_p", colnames(nodes))]
    multi_thresh_shape=T
  } else if (length(shape)>1){
    multi_thresh_shape = T
  } else if(is.numeric(nodes[[shape]]) & grepl("[hinge|upperhinge|segmented]_p$", shape)){
    single_thresh_shape=T
  }
  if(multi_thresh_shape | single_thresh_shape){
    nodes[[paste0(shape, "_bin")]]<-unlist(lapply(shape, function(col){ as.numeric(nodes[[col]]<=thresh_below) }))
    nodes[["significantThresholdModels"]]<-factor(rowSums(as.data.frame(nodes[[paste0(shape, "_bin")]])))
    shape="significantThresholdModels"
  }
  if(is.null(edgeWeight)){
    edgeWeight="NOEDGEWEIGHT"
    edges$NOEDGEWEIGHT=1
  }
  if(!is.null(edgeFilter)){
    if(is.character(edgeFilter)){
      cutoff<-quantile(edges[[edgeWeight]], probs = as.numeric(edgeFilter))
      edges<-edges[edges[[edgeWeight]] >= as.numeric(cutoff), ]
    } else if(is.numeric(edgeFilter)){
      edges<-edges[edges[[edgeWeight]] >= edgeFilter, ]
    } else{stop("edgeFilter must be character or numeric, see ?net.plot for details.")}
  }
  p<-ggplot2::ggplot(nodes)+
    ggplot2::geom_segment(data=edges, ggplot2::aes(x=from.x, xend = to.x, y=from.y,
                                                   yend = to.y, linewidth=.data[[edgeWeight]]),colour="black",alpha=0.1) +
    ggplot2::geom_point(data=nodes, size=size, ggplot2::aes(x=V1,y=V2,
                                                            fill = .data[[fill]], color=.data[[fill]],
                                                            shape=.data[[shape]]), alpha=1, show.legend=T)+
    ggplot2::scale_linewidth(range=c(0.1,1.5))+
    #* note that scaling shape should work, but there is a documented ggplot2 bug where this messes up the legend, so 
    #* until that is fixed I will not specify fillable shapes.
    #scale_shape_discrete(breaks=c(21:(21-1+length(unique(nodes[[shape]]))) ), guide="legend", position="bottom")+
    ggplot2::guides(linewidth="none", shape=ggplot2::guide_legend(nrow=1), fill="none")+
    ggplot2::theme_void()+
    ggplot2::theme(legend.position="bottom")
  if(fill=="NOFILL"){ p<-p+ggplot2::guides(color="none")+ggplot2::scale_color_manual(values="gray80") }
  if(shape=="NOSHAPE"){ p<-p+ggplot2::guides(shape="none") }
  if(fill=="significantThresholdModels"){p<-p+ggplot2::scale_color_continuous(breaks=seq(0,max(nodes$significantThresholdModels, na.rm=T),1)) }
  return(p)
}
