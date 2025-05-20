#' Class \code{thresh} for output from \code{syncomBuildR::thresh} and related functions.
#'
#' Class for changepoint models as used in syncomBuildR.
#'
#' @name thresh-class
#' @docType class
#'
#' @details
#' See \code{methods(class = "thresh")} for an overview of available methods.
#'
#' @slot intercept The intercept term (y value pre changepoint, per hinge model)
#' @slot changepoint The changepoint term (x value at which slope takes effect)
#' @slot slope The slope term
#' @slot phenotype The phenotype (outcome) that the model was fit to.
#' @slot model the model object
#' @slot predictor the ASV/Cluster (predictor) that the model was fit using.
#' @slot data The data used to fit the model
#' @slot type The backend used to fit the model as well as the type of model. Currently "chngptm" is
#' supported.
#' @slot unit The unit/scope. Current options are "individual" for single ASVs or "cluster" for groups
#' of ASVs.
#' @slot control Additional information such as the call.
#'
#' @seealso
#'   \code{\link{thresh}}
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
  #* `Main Points`
  has_adj <- is.null(object$control$p.adjust.method) || object$control$p.adjust.method == "none"
  N_significant <- sum(
    unlist(lapply(object$slope, function(s) {s$padj < 0.05}))
  )
  N_insignificant <- length(object$slope) - N_significant
  cat(
    paste0("Thresh fit with ", object$type, " to ", object$unit, "s\n",
           "\tIncludes ", length(unique(object$predictor)), " predictors over ",
           length(unique(object$phenotype)), " phenotypes using ", nrow(object$data), " observations\n",
           "\t", N_significant, " significant slopes",
           ifelse(has_adj,
                  "",
                  paste0(" (", object$control$p.adjust.method, " adjusted)")
           )
           )
  )
  cat("\n\n")
  if (length(object$intercept) > 1) {
    sum_mat <- cbind(
      sapply(object[c("intercept", "changepoint")],
             function(x) {
               summary(unlist(x))
             }),
      summary(unlist(lapply(object$slope, function(x) {x$est}))),
      summary(unlist(lapply(object$slope, function(x) {x$padj})))
    )
    colnames(sum_mat)[3:4] <- c("slope", "pval")
    print(sum_mat)
  } else {
    sum_mat <- x
  }
  return(invisible(sum_mat))
}

#' Bracket subsetting for thresh objects
#' @aliases `[.thresh`
#'
#' Thresh objects have variable lengths and it's common to want to subset them based on
#' phenotypes or predictors. Specialized bracket methods are implemented to support that.
#' @param x thresh object
#' @param i indexing, logical and numeric vectors are handled specially to subset the thresh object.
#' Other options are handled per a normal list
#' @examples
#' asv$biomass_z <- rnorm(nrow(asv))
#' tm <- thresh(asv, "biomass_z")
#' tm["predictor"] # pulls list
#' tm[1:2] # subsets lists to make a smaller thresh object
#' @method [ thresh
#' @export

"[.thresh" <- function(x, i) {
  if (!is.logical(i) && !is.numeric(i)) {
    return(NextMethod())
  }
  # apply subsetting to subsettable slots
  x[x$control$subsettable] <- lapply(x[x$control$subsettable], function(j) {
    j[i]
    })
  # other slots don't change
  return(x)
}
