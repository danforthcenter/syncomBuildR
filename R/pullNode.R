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
#' @param plot Logical, should data be plotted using \link{net.plot}.
#' @param keepNames Logical, should other nodes be filled in the plot as their node name (TRUE) or all
#' as "other" (FALSE)?
#' Defaults to FALSE.
#'
#' @return A list with a subset network and optionally a ggplot.
#'
#'
#' @examples
#'
#' # a<-qc(); b<-cal(a); c<-thresh(b); d<-asvDist(a) ; e<-net(d, thresh = c)
#' print(load("/home/jsumner/Desktop/stargate/SINC/sincUtils/syncomBuilder/net_output.rdata"))
#'
#' pullNode(net_data, node = "ASV10", edge = "spearman", plot = TRUE, nodeCol = "asv")
#'
#' @export

pullNode <- function(net, node, edge = NULL, edgeFilter = NULL,
                     plot = TRUE, nodeCol = "asv", keepNames = FALSE) {
  #* grab network components
  nodes <- net$nodes
  edges <- net$edges
  #* filter edges for node connection
  edges_sub <- edges[edges$from %in% node | edges$to %in% node, ]
  #* filter and sort edges by some metric/value
  if (!is.null(edge)) {
    edges_sub <- edges_sub[order(edges_sub[[edge]], decreasing = TRUE), ]
    if (!is.null(edgeFilter)) {
      if (is.character(edgeFilter)) {
        cutoff <- quantile(edges_sub[[edge]], probs = as.numeric(edgeFilter))
        edges_sub <- edges_sub[edges_sub[[edge]] >= as.numeric(cutoff), ]
      } else if (is.numeric(edgeFilter)) {
        edges_sub <- edges_sub[edges_sub[[edge]] >= edgeFilter, ]
      } else {
        stop("edgeFilter must be character or numeric.")
      }
    }
  }
  #* filter nodes for the desired node
  nodes_sub <- rbind(
    nodes[nodes[[nodeCol]] %in% node, ],
    nodes[nodes[[nodeCol]] %in% c(edges_sub$to, edges_sub$from), ]
  )
  nodes_sub <- nodes_sub[!duplicated(nodes_sub), ]
  #* plotting
  if (plot) {
    if ("netNumber" %in% colnames(nodes_sub) && length(unique(nodes_sub$netNumber)) > 1) {
      facet <- "netNumber"
    } else {
      facet <- NULL
    }
    if (keepNames) {
      nodes_sub$fill <- nodes_sub[[nodeCol]]
    } else {
      nodes_sub$fill <- ifelse(nodes_sub[[nodeCol]] %in% node,
        paste0(node, collapse = ", "), "other"
      )
    }
    p <- net.plot(list("nodes" = nodes_sub, "edges" = edges_sub),
      fill = "fill", size = 3, edgeWeight = edge,
      edgeFilter = NULL, facet = facet
    )
    out <- list("net" = list("nodes" = nodes_sub, "edges" = edges_sub), "plot" = p)
  } else {
    out <- list("nodes" = nodes_sub, "edges" = edges_sub)
  }
  #* return
  return(out)
}
