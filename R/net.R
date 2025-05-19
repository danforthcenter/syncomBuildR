#' Function to perform network analysis on a distance/dissimilarity matrix of a microbial community.
#'
#'
#' @param df Long dataframe such as that returned from \link{asvDist}.
#' @param metadata Optional dataframe of metadata to attach to nodes. Metadata should belong to data at
#' the node level, so taxonomy for ASV nodes or genotype per sample nodes, but not genotype per ASV
#' since an ASV is present in many samples.
#' @param edge Column name of df to use for edge weighting. Typically this is the same as \code{method}
#' from \code{asvDist} used to make df.
#' @param thresh Optional output from \link{thresh}. If provided then P values of changepoints
#' will be added as node metadata.
#' @param metadata_join Column name of metadata that identifies nodes. Defaults to "asv".
#' @param thresh_join Column name of thresh that identifies nodes. Defaults to "asv" and there are no
#' cases currently where we expect it to be a good idea to change this.
#' @keywords network, changepoint
#' @importFrom igraph graph_from_data_frame graph_from_adjacency_matrix layout.auto get.data.frame
#' as_data_frame betweenness degree E strength harmonic_centrality eigen_centrality
#' authority_score page_rank
#' @import data.table
#' @return An \code{scbnet} object.
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
#' net_data <- asvNet(sp_dist, taxa_df, edge = "spearman_similarity")
#'
#' @export

asvNet <- function(df, metadata = NULL, edge = "spearman_similarity", thresh = NULL,
                   metadata_join = "asv", thresh_join = "asv") {
  if (is.data.frame(df)) {
    i <- unlist(lapply(seq_len(nrow(df)), function(i) {
      paste0(sort(as.character(df[i, 1:3])), collapse = ".")
    }))
    g <- igraph::graph_from_data_frame(df[!duplicated(i), ], directed = FALSE)
  } else if (is.matrix(df)) {
    g <- igraph::graph_from_adjacency_matrix(df, "undirected")
  }
  nd <- as.data.frame(igraph::layout.auto(g))
  eg <- igraph::as_data_frame(g, "edges")
  #* link metadata to nodes
  nd[[metadata_join]] <- igraph::as_data_frame(g, "vertices")$name
  if (!is.null(metadata)) {
    nd <- merge(nd, metadata, by = metadata_join)
  }
  #* link thresh to nodes
  if (!is.null(thresh)) {
    agThresh <- aggregate(p.value. ~ asv + phenotype + model, thresh[thresh$Source != "(Intercept)", ],
      FUN = identity
    )
    subThresh <- data.table::dcast(data.table::as.data.table(agThresh), asv ~ phenotype + model,
      value.var = "p.value."
    )
    colnames(subThresh)[2:ncol(subThresh)] <- paste0(colnames(subThresh)[2:ncol(subThresh)], "_p")
    nd <- merge(nd, subThresh, by = thresh_join)
  }
  #* Calculate network metrics
  nd$betweenness <- igraph::betweenness(g)
  nd$degree <- igraph::degree(g)
  if (!is.null(edge) && is.numeric(eg[[edge]])) {
    igraph::E(g)$weight <- eg[[edge]] + 0.1
    nd$strength <- igraph::strength(g)
  }
  nd$harmonic_centrality <- igraph::harmonic_centrality(g)
  nd$eigen_centrality <- igraph::eigen_centrality(g)[[1]]
  nd$authority_score <- igraph::authority_score(g)[[1]]
  nd$page_rank <- igraph::page_rank(g)[[1]]
  nd$k_coreness <- igraph::coreness(g)[[1]]
  #* add coordinates for plotting edges
  eg$from.x <- nd$V1[match(eg$from, nd[[metadata_join]])]
  eg$from.y <- nd$V2[match(eg$from, nd[[metadata_join]])]
  eg$to.x <- nd$V1[match(eg$to, nd[[metadata_join]])]
  eg$to.y <- nd$V2[match(eg$to, nd[[metadata_join]])]
  eg <- eg[!duplicated(eg), ]
  nd <- nd[!duplicated(nd), ]
  out <- as.scbnet(list("nodes" = nd, "edges" = eg, "graph" = g))
  return(out)
}
