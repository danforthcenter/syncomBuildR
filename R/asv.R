#' Example ASV count table data
#'
#' A very small toy asv table to use in examples.
#'
#' @format ## `asv`
#' A data frame with 32 rows and 13 columns:
#' \describe{
#'   \item{plot}{Some sample metadata such a plot in a field, etc}
#'   \item{tissue}{More sample metadata such as the tissue type}
#'   \item{sample}{A unique identifier for each sample (row)}
#'   \item{ASV1, ASV2, ..., ASV10}{Simulated ASV Counts}
#' }
#' @source Made from the example for \link{qc}
"asv"


#' Example ASV count table data
#'
#' A set of taxa information for use in examples.
#'
#' @format ## `taxa`
#' A matrix with 10 rows and 6 columns, corresponding to the asv data:
#' \describe{
#'   \item{Kingdom}{Simulated kingdom data}
#'   \item{Phylum}{Simulated phylum data}
#'   \item{Class}{Simulated class data}
#'   \item{Order}{Simulated order data}
#'   \item{Family}{Simulated family data}
#'   \item{Genus}{Simulated genus data}
#' }
#' @source Pulled from a SINC dataset
"taxa"
