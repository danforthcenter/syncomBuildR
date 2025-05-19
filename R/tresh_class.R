#' Class \code{thresh} for output from the \code{syncomBuildR::asvNet} function.
#'
#' Class for networks as used in syncomBuildR.
#'
#' @name thresh-class
#' @docType class
#'
#' @details
#' See \code{methods(class = "thresh")} for an overview of available methods.
#'
#' @slot nodes Data frame of node information
#' @slot edges Data frame of edge information
#' @slot graph igraph object
#' @slot control Additional (optional) information about the network
#'
#' @seealso
#'   \code{\link{asvNet}}
#'
NULL

as.thresh <- function(x) {
  class(x) <- "thresh"
  return(x)
}

#' Print an \code{thresh} object.
#'
#' @aliases print.thresh
#'
#' @param x An object of class \code{thresh}.
#' @param ... further arguments, passed to print.default.
#'
#' @seealso \code{\link{summary.thresh}}
#' @method print thresh
#' @export
print.thresh <- function(x, ...) {
  return(summary.thresh(x, ...))
}


#' Summarize an \code{thresh} object.
#'
#' @aliases summary.thresh
#'
#' @param object An object of class \code{thresh}.
#' @param ... further arguments, passed to summary.default
#'
#' @method summary thresh
#' @importFrom utils head tail
#' @export

summary.thresh <- function(object, ...) {
  #* 
  return(invisible(object))
}



