#' Function to cluster networks generated from \code{asvNet}.
#' 
#' 
#' @param net Object returned from \code{\link{asvNet}}
#' @param method Method to use for clustering. Current supported options are components, dbscan, and kmeans. See details.
#' @param ... Additional arguments passed to function picked by method.
#' @keywords network, changepoint
#' @import igraph
#' @import dbscan
#' @return A named list (same as net) with three elements:
#' "Nodes" is a dataframe of nodes and their metadata, now with an added cluster.
#' "Edges" is a dataframe of edges connecting nodes.
#' "graph" is the igraph object used to generate the dataframes.
#' 
#' @details Each method will use a different function to cluster data according to the layout in the graph. Note that layouts in the graph are determined by \code{igraph::layout.auto}.
#' "component" uses \code{igraph::components} to cluster data and requires no additional arguments.
#' "dbscan" uses \code{dbscan::dbscan} to cluster data. This requires at least that the eps argument is set. See \code{?dbscan::dbscan}.
#' "kmeans" uses \code{stats::kmeans} to cluster data. This requires at least that the centers argument is set.
#' 
#' @examples 
#' 
#' # a<-qc(); b<-cal(a); c<-thresh(b); d<-asvDist(a) ; e<-net(d, thresh = c)
#' print(load("/home/jsumner/Desktop/stargate/SINC/sincUtils/syncomBuilder/net_output.rdata"))
#' table(netClust(net=net_data, "component")[["nodes"]]$component_cluster)
#' table(netClust(net=net_data, "dbscan", eps=3)[["nodes"]]$dbscan_cluster)
#' table(netClust(net=net_data, "kmeans", centers=3)[["nodes"]]$kmeans_cluster)
#' 
#' net_data<-netClust(net=net_data, "component")
#' net_data<-netClust(net=net_data, "dbscan", eps=3)
#' net_data<-netClust(net=net_data, "kmeans", centers=3)
#' 
#' 
#' @export
#' 

netClust<-function(net, method="components", ...){
  method=match.arg(method, c("components", "dbscan", "kmeans"))
  if(method == "components"){
    net[["nodes"]]$component_cluster<-igraph::components(net[["graph"]], ...)$membership 
  } else if(method=="dbscan"){
    net[["nodes"]]$dbscan_cluster<-dbscan::dbscan( net[["nodes"]][,c("V1", "V2")], ...)$cluster 
  } else if(method == "kmeans"){
    net[["nodes"]]$kmeans_cluster<-kmeans(net[["nodes"]][,c("V1", "V2")], ...)$cluster
  }
  return(net)
}

