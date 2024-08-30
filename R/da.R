#' Function to assess differential abundance of microbe counts or pseudo counts
#'
#' @param df ASV table with metadata and counts
#' @param col ASV column to model
#' @param predictors Column names of predictors, can be discrete or continuous.
#' @param zi_cutoff Proportion of zeros to consider the distribution to be zero inflated,
#' defaults to 0.1.
#'
#'
#' @importFrom pscl zeroinfl
#' @importFrom MASS glm.nb
#'
#' @return A data frame of model results.
#'
#' @examples
#'
#' # Using uncalibrated data
#'
#' print(load("~/scripts/SINC/sincUtils/syncomBuilder/qc_output.rdata"))
#' asv[1:10, 1:10]
#'
#' df <- asv[, 1:20]
#' sapply(df[, 3:20], function(c) mean(c == 0))
#'
#' col <- "ASV9"
#' predictors <- "tissue"
#' zi_cutoff <- 0.1
#'
#' df[[predictors]] <- relevel(factor(df[[predictors]]), "AS")
#'
#' da(df, col, predictors, zi_cutoff)
#' @export

da <- function(df, col, predictors, zi_cutoff = 0.1) {
  if (any(missing(df), missing(col), missing(predictors))) {
    stop("df, col, and predictors must be specified")
  }
  if (mean(df[[col]] == 0) > zi_cutoff) {
    mod_type <- "zinb"
  } else {
    mod_type <- "nb"
  }

  if (mod_type == "zinb") {
    m <- pscl::zeroinfl(as.formula(paste0(col, "~", paste(predictors, collapse = "+"))),
      data = df,
      dist = "negbin", link = "logit"
    )
    dc <- as.data.frame(summary(m)$coefficients$count)
    dc$comp <- "count"
    dc$coef <- rownames(dc)
    dz <- as.data.frame(summary(m)$coefficients$zero)
    dz$comp <- "zero"
    dz$coef <- rownames(dz)
    m_df <- rbind(dc, dz)
    m_df$model <- "zinb"
    rownames(m_df) <- NULL
    m_df
  } else {
    m2 <- MASS::glm.nb(as.formula(paste0(col, "~", paste(predictors, collapse = "+"))), data = df)
    m_df <- as.data.frame(summary(m2)$coefficients)
    m_df$comp <- "count"
    m_df$model <- "zinb"
    m_df$coef <- rownames(m_df)
    rownames(m_df) <- NULL
    m_df
  }
  return(m_df)
}
