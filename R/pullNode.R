#' Function to highlight subsets of networks generated from \code{asvNet} or \link{netcomi2scb}.
#' 
#' 
#' @param net Network list returned from \link{asvNet} or \link{netcomi2scb}.
#' @param node A node name to extract connections to.
#' @param edge Optional weighting for edges. Must be present in the "edges" of net. Default of NULL will show equal size edges between all connected nodes.
#' @param edgeFilter Optional value to filter edges for. If non-NULL then only edges with edgeWeight greater than this value are kept.
#' This can be a character vector or a numeric.
#' Character vectors are interpreted as quantiles ("0.5" corresponds to the top 50 percent are kept).
#' @param plot Logical, should data be plotted using \link{net.plot}.
#' 
#' @return A list with a subset network and optionally a ggplot.
#' 
#' 
#' @examples 
#' 
#' # a<-qc(); b<-cal(a); c<-thresh(b); d<-asvDist(a) ; e<-net(d, thresh = c)
#' print(load("/home/jsumner/Desktop/stargate/SINC/sincUtils/syncomBuilder/net_output.rdata"))
#' 
#' pullNode(net_data, node = "ASV10", edge = "spearman", plot=TRUE, nodeCol="asv")
#' 
#' @export

pullNode <- function(net, node, edge = NULL, edgeFilter=NULL, plot=TRUE, nodeCol = "name"){
  #* grab network components
  nodes <- net$nodes
  edges <- net$edges
  #* filter nodes for the desired node
  edges_sub <- edges[edges$from == node | edges$to ==node, ]
  nodes_sub <- nodes[nodes[[nodeCol]] == node | nodes[[nodeCol]] %in% edges_sub$to, ]
  #* filter and sort edges
  if(!is.null(edge)){
    edges_sub <- edges_sub[order(edges_sub[[edge]], decreasing = TRUE), ]
    if(!is.null(edgeFilter)){
      if(is.character(edgeFilter)){
        cutoff<-quantile(edges_sub[[edge]], probs = as.numeric(edgeFilter))
        edges_sub<-edges_sub[edges_sub[[edge]] >= as.numeric(cutoff), ]
      } else if(is.numeric(edgeFilter)){
        edges_sub<-edges_sub[edges_sub[[edge]] >= edgeFilter, ]
      } else{stop("edgeFilter must be character or numeric.")}
    }
  }
  #* plotting
  if(plot){
    if("netNumber" %in% colnames(nodes_sub) && length(unique(nodes_sub$netNumber))>1 ){
      facet <- "netNumber"
    } else {facet <- NULL}
    nodes_sub$fill <- ifelse(nodes_sub[[nodeCol]]==node, node, "other")
    p <- net.plot(list("nodes"=nodes_sub, "edges"=edges_sub), fill = "fill", size=3, edgeWeight=edge,
                  edgeFilter = edgeFilter, facet = facet)
    out <- list("net" = list("nodes"=nodes_sub, "edges"=edges_sub), "plot"=p)
  } else {
    out <- list("nodes"=nodes_sub, "edges"=edges_sub)
  }
  #* return
  return(out)
  
}

