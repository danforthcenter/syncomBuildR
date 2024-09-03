#' Function to plot networks generated from \code{asvNet} with various emphases.
#'
#'
#' @param net Object returned from \link{asvNet}.
#' @param fill Optional column name to fill points by. Accepts single column names, vectors of p-value
#' columns, or "thresh" which will match all p-value columns if the network was fit with threshold model
#' data from \link{thresh}.
#' @param shape Optional column name to use for node shapes. Accepts the same options as fill.
#' @param size Size for points. Passed to ggplot2::geom_point.
#' @param edgeWeight Optional weighting for edges. Must be present in the "edges" of net. Default of
#' NULL will show equal size edges between all connected nodes.
#' @param edgeFilter Optional value to filter edges for. If non-NULL then only edges with edgeWeight
#' greater than this value are kept.
#' This can be a character vector or a numeric.
#' Character vectors are interpreted as quantiles ("0.5" corresponds to the top 50 percent are kept).
#' @param thresh_below Significant cutoff if p-value columns are used for fill or shape.
#' Defaults to 0.05.
#' @param facet Optionally a variable to facet the plot on. This is meant to be used to separate
#' multiple networks created from \link{netcomi2scb}
#' in which case "netNumber" should be used.
#' @param method Method for visualization, defaults to "ggplot" but also accepts "plotly".
#' @import ggplot2
#' @return A plot of the kind determined by the method argument.
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
#'
#' sp_dist <- asvDist(asv, method = "spearman", clr_transform = TRUE, edgeFilter = 0.5)
#' net_data <- asvNet(sp_dist, taxa_df, edge = "spearman")
#'
#' net.plot(net_data, size = 3, edgeWeight = "spearman", thresh_below = 0.05)
#'
#' @export
#'

net.plot <- function(net, fill = NULL, shape = NULL, size = 3, edgeWeight = NULL,
                     edgeFilter = NULL, thresh_below = 0.05, facet = NULL, method = "ggplot") {
  nodes <- net[["nodes"]]
  edges <- net[["edges"]]
  multi_thresh_fill <- FALSE
  single_thresh_fill <- FALSE
  multi_thresh_shape <- FALSE
  single_thresh_shape <- FALSE
  #* make fill work
  fill_helper_out <- .netPlotFillHelper(
    multi_thresh_fill, single_thresh_fill, fill,
    edges, nodes
  )
  nodes <- fill_helper_out$nodes
  edges <- fill_helper_out$edges
  multi_thresh_fill <- fill_helper_out$mtf
  single_thresh_fill <- fill_helper_out$stf
  fill <- fill_helper_out$fill

  if (multi_thresh_fill || single_thresh_fill) {
    nodes[[paste0(fill, "_bin")]] <- unlist(lapply(fill, function(col) {
      as.numeric(nodes[[col]] <= thresh_below)
    }))
    nodes[["significantThresholdModels"]] <- rowSums(as.data.frame(nodes[[paste0(fill, "_bin")]]))
    fill <- "significantThresholdModels"
  }
  #* make shape work
  shape_helper_out <- .netPlotShapeHelper(
    multi_thresh_shape, single_thresh_shape, shape,
    edges, nodes
  )
  nodes <- shape_helper_out$nodes
  edges <- shape_helper_out$edges
  multi_thresh_shape <- shape_helper_out$mts
  single_thresh_shape <- shape_helper_out$sts
  shape <- shape_helper_out$shape

  if (multi_thresh_shape || single_thresh_shape) {
    nodes[[paste0(shape, "_bin")]] <- unlist(lapply(shape, function(col) {
      as.numeric(nodes[[col]] <= thresh_below)
    }))
    nodes[["significantThresholdModels"]] <- factor(
      rowSums(as.data.frame(nodes[[paste0(shape, "_bin")]]))
    )
    shape <- "significantThresholdModels"
  }
  if (is.null(edgeWeight)) {
    edgeWeight <- "NOEDGEWEIGHT"
    edges$NOEDGEWEIGHT <- 1
  }
  if (!is.null(edgeFilter)) {
    if (is.character(edgeFilter)) {
      cutoff <- quantile(edges[[edgeWeight]], probs = as.numeric(edgeFilter))
      edges <- edges[edges[[edgeWeight]] >= as.numeric(cutoff), ]
    } else if (is.numeric(edgeFilter)) {
      edges <- edges[edges[[edgeWeight]] >= edgeFilter, ]
    } else {
      stop("edgeFilter must be character or numeric, see ?net.plot for details.")
    }
  }
  if (method == "plotly") {
    if (is.null(facet)) {
      p <- .plotlyNetPlot(nodes, edges, edgeWeight, fill, shape, facet)
    } else {
      p <- .facetedPlotlyNetPlot(nodes, edges, edgeWeight, fill, shape, facet)
    }
  } else {
    p <- .ggNetPlot(nodes, edges, edgeWeight, fill, shape, facet, size)
  }
  return(p)
}

#' @keywords internal
#' @noRd

.netPlotFillHelper <- function(multi_thresh_fill, single_thresh_fill, fill,
                               edges, nodes) {
  if (is.null(fill)) {
    fill <- "NOFILL"
    edges$NOFILL <- "a"
    nodes$NOFILL <- "a"
  } else if (fill == "thresh") {
    fill <- colnames(nodes)[grepl("[hinge|upperhinge|segmented]_p", colnames(nodes))]
    multi_thresh_fill <- TRUE
  } else if (length(fill) > 1) {
    multi_thresh_fill <- TRUE
  } else if (is.numeric(nodes[[fill]]) && grepl("[hinge|upperhinge|segmented]_p$", fill)) {
    single_thresh_fill <- TRUE
  }
  return(list(
    "mtf" = multi_thresh_fill,
    "stf" = single_thresh_fill,
    "edges" = edges,
    "nodes" = nodes,
    "fill" = fill
  ))
}

#' @keywords internal
#' @noRd

.netPlotShapeHelper <- function(multi_thresh_shape, single_thresh_shape, shape,
                                edges, nodes) {
  if (is.null(shape)) {
    shape <- "NOSHAPE"
    nodes$NOSHAPE <- "a"
  } else if (shape == "thresh") {
    shape <- colnames(nodes)[grepl("[hinge|upperhinge|segmented]_p", colnames(nodes))]
    multi_thresh_shape <- TRUE
  } else if (length(shape) > 1) {
    multi_thresh_shape <- TRUE
  } else if (is.numeric(nodes[[shape]]) && grepl("[hinge|upperhinge|segmented]_p$", shape)) {
    single_thresh_shape <- TRUE
  }
  return(list(
    "mts" = multi_thresh_shape,
    "sts" = single_thresh_shape,
    "edges" = edges,
    "nodes" = nodes,
    "shape" = shape
  ))
}

#' @keywords internal
#' @noRd

.ggNetPlot <- function(nodes, edges, edgeWeight, fill, shape, facet, size) {
  p <- ggplot2::ggplot(nodes) +
    ggplot2::geom_segment(
      data = edges, ggplot2::aes(
        x = .data[["from.x"]], xend = .data[["to.x"]], y = .data[["from.y"]],
        yend = .data[["to.y"]], linewidth = .data[[edgeWeight]]
      ),
      colour = "black", alpha = 0.1
    ) +
    ggplot2::geom_point(
      data = nodes, size = size, ggplot2::aes(
        x = .data[["V1"]], y = .data[["V2"]],
        fill = .data[[fill]], color = .data[[fill]],
        shape = .data[[shape]]
      ),
      alpha = 1, show.legend = TRUE
    ) +
    ggplot2::scale_linewidth(range = c(0.1, 1.5)) +
    ggplot2::guides(linewidth = "none", shape = ggplot2::guide_legend(nrow = 1), fill = "none") +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "bottom")
  if (fill == "NOFILL") {
    p <- p + ggplot2::guides(color = "none") + ggplot2::scale_color_manual(values = "gray80")
  }
  if (shape == "NOSHAPE") {
    p <- p + ggplot2::guides(shape = "none")
  }
  if (fill == "significantThresholdModels") {
    p <- p + ggplot2::scale_color_continuous(
      breaks = seq(0, max(nodes$significantThresholdModels, na.rm = TRUE), 1)
    )
  }
  if (!is.null(facet)) {
    p <- p + ggplot2::facet_wrap(as.formula(paste0("~", facet))) +
      ggplot2::theme(
        panel.background = element_rect(color = "black"),
        strip.background = element_rect(fill = "gray40"),
        strip.text = element_text(color = "white", size = 12)
      )
  }
  return(p)
}

#' @keywords internal
#' @importFrom plotly plot_ly add_segments add_markers layout
#' @noRd

.plotlyNetPlot <- function(nodes, edges, edgeWeight = NULL, fill = NULL, shape = NULL, facet = NULL) {
  p <- plotly::plot_ly(data = nodes)
  p <- plotly::add_segments(
    p = p, data = edges,
    x = ~from.x,
    y = ~from.y,
    xend = ~to.x,
    yend = ~to.y,
    text = ~ paste0(
      from, " to ", to, "\n", edgeWeight, ": ",
      round(edges[[edgeWeight]], 3)
    ),
    hoverinfo = "text",
    line = list(width = 0.1),
    name = "Edges",
    inherit = FALSE
  )
  if (fill == "NOFILL") {
    p <- plotly::add_markers(
      p = p, data = nodes,
      x = ~V1,
      y = ~V2,
      text = ~ paste0(asv, "\n", "strength: ", round(nodes[["strength"]], 1)),
      hoverinfo = "text",
      type = "scatter",
      mode = "markers",
      name = "Nodes",
      marker = list(
        size = ~ log(strength, base = exp(1) / 2),
        opacity = 0.75,
        color = "black",
        line = list(color = "black")
      )
    )
  } else {
    p <- plotly::add_markers(
      p = p, data = nodes,
      x = ~V1,
      y = ~V2,
      text = ~ paste0(
        asv, "\n", "strength: ",
        round(nodes[["strength"]], 1), "\n",
        nodes[[fill]]
      ),
      hoverinfo = "text",
      type = "scatter",
      mode = "markers",
      name = "Nodes",
      marker = list(
        size = ~ log(strength, base = exp(1) / 2),
        opacity = 0.75,
        color = ~ nodes[[fill]],
        line = list(color = "black")
      )
    )
  }
  p <- plotly::layout(p,
    xaxis = list(visible = FALSE),
    yaxis = list(visible = FALSE)
  )
  return(p)
}

#' @keywords internal
#' @importFrom plotly subplot
#' @noRd

.facetedPlotlyNetPlot <- function(nodes, edges, edgeWeight = NULL, fill = NULL,
                                  shape = NULL, facet = NULL) {
  nodes <- nodes[!is.na(nodes[[facet]]), ]
  ps <- lapply(unique(nodes[[facet]]), function(fct) {
    sub_nodes <- nodes[nodes[[facet]] == fct, ]
    sub_edges <- edges[edges$to %in% nodes$asv || edges$from %in% nodes$asv, ]
    .plotlyNetPlot(sub_nodes, sub_edges,
      edgeWeight = edgeWeight, fill = fill,
      shape = shape, facet = facet
    )
  })
  p <- plotly::subplot(ps,
    nrows = ceiling(length(unique(nodes[[facet]])) / 3),
    shareX = TRUE, shareY = TRUE
  )
  return(p)
}
