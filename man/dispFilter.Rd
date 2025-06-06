% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/disparityFilter.R
\name{dispFilter}
\alias{dispFilter}
\title{Disparity Filtering for microbial networks}
\usage{
dispFilter(
  net,
  weights = NULL,
  alpha = 0.05,
  cores = getOption("mc.cores", 1),
  replot = TRUE
)
}
\arguments{
\item{net}{Object returned from \link{asvNet}.}

\item{weights}{Optional weights specification. This defaults to NULL in which case "weight" is used
as calculated by igraph. This can also be a numeric vector of edge weights
(ordered for \code{net$edges}) or a column name from \code{net$edges}.}

\item{alpha}{the significance level to filter edges for. Defaults to 0.05 for no serious reason.}

\item{cores}{Number of cores to run in parallel, defaults to 1 if "mc.cores" option is unset.}

\item{replot}{Logical, should nodes be rearranged to represent the network better visually?
Defaults to TRUE.}
}
\value{
A modified version of net with filtered edges (and nodes if any were now isolated).
}
\description{
Disparity Filtering for microbial networks
}
\examples{

taxa <- c(
  "Bacteria", "Proteobacteria", "Betaproteobacteria", "Burkholderiales",
  "Burkholderiaceae", "Paraburkholderia", NA
)
taxa <- matrix(rep(taxa, 10), nrow = 10, byrow = TRUE)
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(taxa) <- paste0("ASV", 1:10)
# taxonomy data if used should have ASV names explicitly as a column
taxa_df <- as.data.frame(taxa)
taxa_df$asv <- rownames(taxa_df)

sp_dist <- asvDist(asv, method = "spearman", clr_transform = TRUE, edgeFilter = 0.5)
net_data <- asvNet(sp_dist, taxa_df, edge = "spearman_similarity")
# we could use spearman_similarity, but not specifying will use the graph's default weights
net_data_disp_filt <- dispFilter(net_data, alpha = 0.4)

}
