#' Function to assess differential abundance of microbe counts or pseudo counts
#'
#' @param df ASV table with metadata and counts
#' @param col ASV column to model
#' @param predictors Column names of predictors, can be discrete or continuous.
#' @param zi_cutoff Proportion of zeros to consider the distribution to be zero inflated,
#' defaults to 0.1.
#' @param priors Prior distributions as \code{brmsprior} objects, defaults are shown in details.
#' If left NULL (the default) then wide student priors will be used
#' (\code{brms::set_prior('student_t(3,0,5)',  class="b")}).
#' @param intercept Logical, should an intercept be used? Defaults to TRUE.
#' @param backend Argument to control \code{brms::brm()} backend interface to Stan.
#' @param chains Argument to control number of chains used by \code{brms::brm()}.
#' @param cores Argument to control number of cores to use in parallel \code{brms::brm()}.
#' @param iter Argument to control number of iterations \code{brms::brm()}.
#' @param ... Further arguments passed to \code{brms::brm()}.
#'
#' @details
#'
#' The prior distributions for a ZINB model are shown below, the only difference with
#' an NB model is that there is no zi parameter.
#'
#'                   prior     class      coef group resp dpar nlpar lb ub       source
#'                   (flat)         b                                            default
#'                   (flat)         b levels                                (vectorized)
#'   student_t(3, 3.1, 2.5) Intercept                                            default
#'        gamma(0.01, 0.01)     shape                                  0         default
#'               beta(1, 1)        zi                                  0  1      default
#'
#' @return A brmsfit object
#' @examples
#'
#'
#' if ("brms" %in% installed.packages()) {
#'   print(load("~/Desktop/stargate/SINC/sincUtils/syncomBuilder/qc_output.rdata"))
#'   asv[1:10, 1:10]
#'
#'   df <- asv[, 1:20]
#'   sapply(df[, 3:20], function(c) mean(c == 0))
#'
#'   col <- "ASV9"
#'   predictors <- "tissue"
#'   zi_cutoff <- 0.1
#'
#'   df[[predictors]] <- relevel(factor(df[[predictors]]), "AS")
#'
#'   ex <- b_da(df, col, predictors, zi_cutoff, intercept = FALSE)
#' }
#'
b_da <- function(df, col, predictors, zi_cutoff = 0.1, priors = NULL,
                 intercept = TRUE, backend = "cmdstanr",
                 cores = 4, chains = 4, iter = 2000, ...) {
  if (any(missing(df), missing(col), missing(predictors))) {
    stop("df, col, and predictors must be specified")
  }
  if (mean(df[[col]] == 0) > zi_cutoff) {
    mod_type <- "zinb"
  } else {
    mod_type <- "nb"
  }

  if (is.null(priors)) {
    priors <- brms::set_prior("student_t(3,0,5)", class = "b")
  }

  if (!intercept) {
    predictors <- c(0, predictors)
  }

  if (mod_type == "zinb") {
    form <- brms::bf(as.formula(paste0(col, "~", paste(predictors, collapse = "+"))),
      family = "zero_inflated_negbinomial"
    )
    mod <- brms::brm(form, df,
      prior = priors, backend = backend, chains = chains, cores = cores,
      iter = iter
    )
  } else {
    form <- brms::bf(as.formula(paste0(col, "~", paste(predictors, collapse = "+"))),
      family = "negbinomial"
    )
    mod <- brms::brm(form, df,
      prior = priors, backend = backend, chains = chains, cores = cores,
      iter = iter
    )
  }
  return(mod)
}
