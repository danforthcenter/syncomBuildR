#' Function to cluster networks generated from \code{asvNet}.
#'
#'
#' @param net Object returned from \link{asvNet}
#' @param method Method to use for clustering. This can be a method from "components", "dbscan",
#' and "kmeans" or the output from pullNode, in which case that node and it's connections are labelled
#' as a cluster.
#' @param ... Additional arguments passed to function picked by method.
#' @keywords network, changepoint
#' @importFrom igraph components
#' @importFrom dbscan dbscan
#' @return A named list (same as net) with three elements, same as \link{asvNet}:
#' \itemize{
#'    \item{"Nodes" is a dataframe of nodes and their metadata}
#'    \item{"Edges" is a dataframe of edges connecting nodes.}
#'    \item{"graph" is the igraph object used to generate the dataframes.}
#' }
#'
#' @details Each method will use a different function to cluster data according to the layout
#'          in the graph. Note that layouts in the graph are determined by \code{igraph::layout_nicely}.
#' \itemize{
#'    \item{"component" uses \code{igraph::components} to cluster data and
#'          requires no additional arguments.}
#'    \item{"dbscan" uses \code{dbscan::dbscan} to cluster data.
#'          This requires at least that the eps argument is set. See \code{?dbscan::dbscan}.}
#'    \item{"kmeans" uses \code{stats::kmeans} to cluster data.
#'          This requires at least that the centers argument is set.}
#' }
#'
#'
#' @examples
#'
#' # a<-qc(); b<-cal(a); c<-thresh(b); d<-asvDist(a) ; e<-net(d, thresh = c)
#' print(load("~/scripts/SINC/sincUtils/syncomBuilder/net_output.rdata"))
#' table(netClust(net = net_data, "component")[["nodes"]]$component_cluster)
#' table(netClust(net = net_data, "dbscan", eps = 3)[["nodes"]]$dbscan_cluster)
#' table(netClust(net = net_data, "kmeans", centers = 3)[["nodes"]]$kmeans_cluster)
#'
#' net_data <- netClust(net = net_data, "component")
#' net_data <- netClust(net = net_data, "dbscan", eps = 3)
#' net_data <- netClust(net = net_data, "kmeans", centers = 3)
#'
#' @export
#'

netClust <- function(net, method = "components", ...) {
  if (is.character(method)) {
    if (method == "components") {
      net[["nodes"]]$component_cluster <- as.character(
        igraph::components(net[["graph"]], ...)$membership
      )
    } else if (method == "dbscan") {
      net[["nodes"]]$dbscan_cluster <- as.character(
        dbscan::dbscan(net[["nodes"]][, c("V1", "V2")], ...)$cluster
      )
    } else if (method == "kmeans") {
      net[["nodes"]]$kmeans_cluster <- as.character(
        kmeans(net[["nodes"]][, c("V1", "V2")], ...)$cluster
      )
    }
  } else { # pull node option
    if (any(unlist(lapply(method, function(l) any(ggplot2::is.ggplot(l)))))) {
      method <- method$net
    }
    pulled_nodes <- method$nodes
    net_nodes <- net[["nodes"]]
    fills <- unique(pulled_nodes$fill)
    fills <- fills[!grepl("other", fills)]
    net_nodes[[paste0(
      "pullNode_cluster_",
      paste(fills, collapse = "_")
    )]] <- ifelse(net_nodes$asv %in% pulled_nodes$asv,
      "In", "Out"
    )
    net[["nodes"]] <- net_nodes
  }
  return(net)
}
