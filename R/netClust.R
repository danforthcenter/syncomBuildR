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
#' @return An \code{scbnet} object with cluster data added to the nodes data.frame.
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
#' taxa <- c(
#'   "Bacteria", "Proteobacteria", "Betaproteobacteria", "Burkholderiales",
#'   "Burkholderiaceae", "Paraburkholderia", NA
#' )
#' taxa <- matrix(rep(taxa, 10), nrow = 10, byrow = TRUE)
#' colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
#' rownames(taxa) <- paste0("ASV", 1:10)
#' # taxonomy data if used should have ASV names explicitly as a column
#' taxa_df <- as.data.frame(taxa)
#' taxa_df$asv <- rownames(taxa_df)
#' sp_dist <- asvDist(asv, method = "spearman", clr_transform = TRUE, edgeFilter = 0.5)
#' net_data <- asvNet(sp_dist, taxa_df, edge = "spearman_similarity")
#' net_data <- netClust(net = net_data, "components")
#' net_data <- netClust(net = net_data,  method = pullNode(net_data, node = c("ASV10", "ASV9"),
#'                      edge = "spearman_similarity", plot = TRUE, nodeCol = "asv"))
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
