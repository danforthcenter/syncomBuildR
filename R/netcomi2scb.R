#' Function to translate NetCoMi objects into dataframes usable by other syncombuildr functions.
#'
#' The main purpose of this function is to bring NetCoMi's metrics and testing into a format
#' that allows for more customized visualization such as with \link{net.plot}
#'
#'
#' @param microNetObj A microNetObj object as returned by \code{NetCoMi::netConstruct}
#' @param microNetProps A microNetProps object as returned by \code{NetCoMi::netAnalyze}
#' @param microNetComp A microNetComp object as returned by \code{NetCoMi::netCompare}
#' @param cutoff Optional value to filter edges for. If non-NULL then only edges with edgeWeight
#' greater than this value are kept.
#' This can be a character vector or a numeric.
#' Character vectors are interpreted as quantiles ("0.5" corresponds to the top 50% are kept).
#' Defaults to "0.9" .
#' @param combine Logical, should multiple networks have X/Y coordinates for nodes standardized to
#' make comparable plots later on? Defaults to TRUE.
#' @keywords network, changepoint
#' @importFrom igraph graph_from_data_frame layout_nicely as_data_frame
#' @return If the microNetObj contains one network then a named list with three elements is returned:
#' \itemize{
#'    \item{"Nodes" is a dataframe of nodes and their metadata}
#'    \item{"Edges" is a dataframe of edges connecting nodes.}
#'    \item{"graph" is the igraph object used to generate the dataframes.}
#' }
#' If the microNetObj contains two networks then a list of such lists is returned.
#'
#' @examples
#' \dontrun{
#' if ("NetCoMi" %in% installed.packages()) {
#' microNetObj <- NetCoMi::netConstruct(as.matrix(asv[, grepl("ASV", colnames(asv))]),
#'   measure = "spearman", sparsMethod = "t-test", alpha = 0.7
#' )
#' net <- netcomi2scb(microNetObj)
#' }
#' }
#' @export

netcomi2scb <- function(microNetObj, microNetProps = NULL, microNetComp = NULL, cutoff = "0.9",
                        combine = TRUE) {
  if (microNetObj$twoNets) {
    nNets <- 2
  } else {
    nNets <- 1
  }
  res <- lapply(1:nNets, function(w) {
    mat <- microNetObj[[paste0("adjaMat", w)]]

    long <- do.call(rbind, lapply(seq_len(nrow(mat)), function(i) {
      do.call(rbind, lapply(seq_len(ncol(mat)), function(j) {
        data.frame(value = mat[i, j], from = colnames(mat)[j], to = rownames(mat)[i])
      }))
    }))
    if (!is.null(cutoff)) {
      if (is.character(cutoff)) {
        cutoff_value <- as.numeric(quantile(as.numeric(mat[mat < 1]), as.numeric(cutoff)))
      } else {
        cutoff_value <- cutoff
      }
      long <- long[long$value > cutoff_value, ]
    }
    g <- igraph::graph_from_data_frame(long, "undirected")

    nd <- as.data.frame(igraph::layout_nicely(g)) # layout_as_star , layout_nicely ?
    nd$name <- igraph::as_data_frame(g, "vertices")$name
    eg <- microNetObj[[paste0("edgelist", w)]]
    colnames(eg)[1:2] <- c("from", "to")
    eg$from.x <- nd$V1[match(eg$from, nd[["name"]])]
    eg$from.y <- nd$V2[match(eg$from, nd[["name"]])]
    eg$to.x <- nd$V1[match(eg$to, nd[["name"]])]
    eg$to.y <- nd$V2[match(eg$to, nd[["name"]])]

    if (!is.null(microNetProps)) {
      centralities <- as.data.frame(t(do.call(rbind, microNetProps$centralities)))
      centralities$name <- rownames(centralities)
      nd <- merge(nd, centralities, by = "name")
    }
    if (!is.null(microNetComp)) {
      pvals <- as.data.frame(t(do.call(rbind, microNetComp[["pvalDiffCentrAdjust"]])))
      pvals$name <- rownames(pvals)
      negLog10pvals <- as.data.frame(-log10(t(do.call(rbind, microNetComp[["pvalDiffCentrAdjust"]]))))
      colnames(negLog10pvals) <- paste0("negLog10", colnames(negLog10pvals))
      negLog10pvals$name <- rownames(negLog10pvals)
      nd <- merge(nd, pvals, by = "name")
      nd <- merge(nd, negLog10pvals, by = "name")
    }
    return(list("nodes" = nd, "edges" = eg, "graph" = g))
  })
  if (nNets == 1) {
    res <- res[[1]]
  } else if (combine) {
    res <- .combine_nets(nets = res)
  }
  return(res)
}



#' Helper function to combine multiple netcomi networks into one igraph/dataframe object.
#'
#' The main purpose of this function is to bring NetCoMi's metrics and testing into a format
#' that allows for more customized visualization such as with \link{net.plot}
#'
#' @param nets List in the style of netcomi2scb output.
#' @param metrics Character vector of metrics to combine.
#' @keywords internal
#' @noRd

.combine_nets <- function(nets, metrics = c("degree", "between", "close", "eigen")) {
  nets <- lapply(seq_len(length(nets)), function(i) {
    net <- nets[[i]]
    net <- lapply(net, function(n) {
      if (is.data.frame(n)) {
        n$netNumber <- i
        for (metric in metrics) {
          n[[metric]] <- n[[paste0(metric, i)]]
        }
      }
      return(n)
    })
  })
  for (i in 2:length(nets)) {
    iterNet <- nets[[i]]
    originNet <- nets[[1]]
    iterNet$edges$from.x <- originNet$nodes$V1[match(iterNet$edges$from, originNet$nodes[["name"]])]
    iterNet$edges$from.y <- originNet$nodes$V2[match(iterNet$edges$from, originNet$nodes[["name"]])]
    iterNet$edges$to.x <- originNet$nodes$V1[match(iterNet$edges$to, originNet$nodes[["name"]])]
    iterNet$edges$to.y <- originNet$nodes$V2[match(iterNet$edges$to, originNet$nodes[["name"]])]
    iterNet$nodes$V1 <- originNet$nodes$V1[match(originNet$nodes[["name"]], iterNet$nodes[["name"]])]
    iterNet$nodes$V2 <- originNet$nodes$V2[match(originNet$nodes[["name"]], iterNet$nodes[["name"]])]
    nets[[i]]$edges <- iterNet$edges
    nets[[i]]$nodes <- iterNet$nodes
  }
  allNodes <- do.call(rbind, lapply(nets, function(net) {
    net$nodes
  }))
  allEdges <- do.call(rbind, lapply(nets, function(net) {
    net$edges
  }))

  nets[[length(nets) + 1]] <- list("nodes" = allNodes, "edges" = allEdges, "graph" = NULL)

  return(nets)
}
