#' Function to highlight subsets of networks generated from \code{asvNet} or \link{netcomi2scb}.
#'
#'
#' @param net Network list returned from \link{asvNet} or \link{netcomi2scb}.
#' @param node A node name (or vector of names) to extract connections to.
#' @param edge Optional weighting for edges. Must be present in the "edges" of net. Default of NULL
#' will show equal size edges between all connected nodes.
#' @param edgeFilter Optional value to filter edges for. If non-NULL then only edges with edgeWeight
#' greater than this value are kept.
#' This can be a character vector or a numeric.
#' Character vectors are interpreted as quantiles ("0.5" corresponds to the top 50 percent are kept).
#' @param plot deprecated, use \code{\link{plot.scbnet}}.
#' @param nodeCol Column of node information to find the node in. Defaults to "asv".
#' @param keepNames Logical, should other nodes be filled in the plot as their node name (TRUE) or all
#' as "other" (FALSE)? Defaults to FALSE.
#' @param order The order of connection to include, defaults to 1. If Node A has an Edge to Node B,
#' and Node B has an edge to Node C, but Node A does not have an edge to Node C then A to C is a 2nd
#' order connection, and Node C would only be included if order >= 2.
#' @return An \code{scbnet} object of the subset network
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
#'
#' sub_net <- pullNode(net_data, node = "ASV10",
#' edge = "spearman_similarity", edgeFilter = 0.7, nodeCol = "asv")
#'
#' @export

pullNode <- function(net, node, edge = NULL, edgeFilter = NULL,
                     plot = FALSE, nodeCol = "asv", keepNames = FALSE, order = 1) {
  #* grab network components
  nodes <- net$nodes
  edges <- net$edges
  #* filter and sort edges by some metric/value
  if (!is.null(edge)) {
    edges <- edges[order(edges[[edge]], decreasing = TRUE), ]
    if (!is.null(edgeFilter)) {
      if (is.character(edgeFilter)) {
        cutoff <- quantile(edges[[edge]], probs = as.numeric(edgeFilter))
        edges <- edges[edges[[edge]] >= as.numeric(cutoff), ]
      } else if (is.numeric(edgeFilter)) {
        edges <- edges[edges[[edge]] >= edgeFilter, ]
      } else {
        stop("edgeFilter must be character or numeric.")
      }
    }
  }
  #* order logic
  total_from_names <- node
  for (i in seq_len(order)) {
    edges_sub <- edges[edges$from %in% total_from_names, ]
    if (i > 1) {
      total_from_names <- unique(c(total_from_names, edges_sub$from, edges_sub$to))
    }
  }
  #* filter edges
  edges_sub <- edges[edges$from %in% total_from_names | edges$to %in% total_from_names, ]
  #* filter nodes for the desired node
  nodes_sub <- rbind(
    nodes[nodes[[nodeCol]] %in% node, ],
    nodes[nodes[[nodeCol]] %in% c(edges_sub$to, edges_sub$from), ]
  )
  nodes_sub <- nodes_sub[!duplicated(nodes_sub), ]
  if (keepNames) {
    nodes_sub$fill <- nodes_sub[[nodeCol]]
  } else {
    nodes_sub$fill <- ifelse(nodes_sub[[nodeCol]] %in% node,
                             paste0(node, collapse = ", "), "other"
    )
  }
  #* subset igraph object
  g <- net$graph
  gs <- igraph::induced_subgraph(g, vids = which(names(igraph::V(g)) %in% nodes_sub[[nodeCol]]))
  gs <- igraph::delete_edges(gs, seq_len(length(igraph::E(gs)))) # or graph - E(graph) works
  gs <- igraph::add_edges(gs, edges = c(rbind(edges_sub$from, edges_sub$to)))
  #* plotting
  if (plot) {
    warning("Plot argument is deprecated, use generic plot.scbnet instead.")
    if ("netNumber" %in% colnames(nodes_sub) && length(unique(nodes_sub$netNumber)) > 1) {
      facet <- "netNumber"
    } else {
      facet <- NULL
    }
    p <- plot.scbnet(list("nodes" = nodes_sub, "edges" = edges_sub),
      fill = "fill", size = 3, edgeWeight = edge,
      edgeFilter = NULL, facet = facet
    )
    print(p)
  }
  out <- as.scbnet(list("nodes" = nodes_sub, "edges" = edges_sub, "graph" = gs))
  return(out)
}
