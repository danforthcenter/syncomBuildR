#' Class \code{scbnet} for output from the \code{syncomBuildR::asvNet} function.
#'
#' Class for networks as used in syncomBuildR.
#'
#' @name scbnet-class
#' @docType class
#'
#' @details
#' See \code{methods(class = "scbnet")} for an overview of available methods.
#'
#' @slot nodes Data frame of node information
#' @slot edges Data frame of edge information
#' @slot graph igraph object
#' @slot control Additional information about the network
#'
#' @seealso
#'   \code{\link{asvNet}}
#'
NULL

as.scbnet <- function(x) {
  class(x) <- "scbnet"
  return(x)
}

#' Print an \code{scbnet} object.
#'
#' @aliases print.scbnet
#'
#' @param x An object of class \code{scbnet}.
#' @param ... further arguments, passed to print.default.
#'
#' @seealso \code{\link{summary.scbnet}}
#' @method print scbnet
#' @export
print.scbnet <- function(x, ...) {
  return(summary.scbnet(x, ...))
}


#' Summarize an \code{scbnet} object.
#'
#' @aliases summary.scbnet
#'
#' @param object An object of class \code{scbnet}.
#' @param ... further arguments, passed to summary.default
#'
#' @method summary scbnet
#' @importFrom utils head tail
#' @export

summary.scbnet <- function(object, ...) {
  #* `Dimensions`
  cat(
    paste0(
      nrow(object$nodes), " Nodes connected by ", nrow(object$edges), " edges.\n"
    )
  )
  #* `Get information from nodes`
  node_names <- sort(as.numeric(gsub("[a-zA-Z]*", "", object$nodes$name)), index.return = TRUE)
  if (any(nchar(node_names$x) < 1) || any(duplicated(node_names$x))) {
    node_names <- sort(object$nodes$name)
  } else {
    node_names <- object$nodes$name[node_names$ix]
  }
  cat(
    paste0(
      "Nodes:\n", paste(head(node_names, 3), collapse = ", "), ", ..., ", tail(node_names, 1), "\n"
    )
  )
  #* `Get information from edges`
  cat("Edges:\n")
  summary(object$edges[ , -which(
    colnames(object$edges) %in% c("from", "to", "to.x", "to.y", "from.x", "from.y")
    )]
    )
  #* `Check for control information`
  #* Not sure if this will always exist or if it's something that I'll only make in response to certain
  #* things. Thinking that this holds stuff like where the graph came from and what has happened upstream?
  #*
  #* For now this is pending
  return(invisible(object))
}



