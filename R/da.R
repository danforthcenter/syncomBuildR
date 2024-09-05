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
#' @return A data frame of model results and Log2 fold changes.
#'
#' @examples
#'
#' da(asv, col = "ASV2", predictors = "tissue", zi_cutoff = 0.1)
#'
#' @export

da <- function(df, col, predictors, zi_cutoff = 0.1) {
  if (any(missing(df), missing(col), missing(predictors))) {
    stop("df, col, and predictors must be specified")
  }
  log2_values <- log2(df[[col]]+0.0001)
  vals <- split(log2_values, df[, c(predictors)])
  df[[col]] <- round(df[[col]]) # in case of calibrated pseudo counts
  if (mean(df[[col]] == 0) > zi_cutoff) {
    mod_type <- "zinb"
  } else {
    mod_type <- "nb"
  }
  if (mod_type == "zinb") {
    m <- pscl::zeroinfl(as.formula(paste0(col, "~", paste(predictors, collapse = "*"))),
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
    log2FC <- unlist(lapply(seq_along(vals), function(i) {
      mean(vals[[i]]) - mean(vals[[1]])
    }))
    log2_mu <- unlist(lapply(vals, mean))
    m_df$log2_mu <- c(log2_mu, NA, log2_mu)
    m_df$log2FC <- c(log2FC, NA, log2FC)
  } else {
    m2 <- MASS::glm.nb(as.formula(paste0(col, "~", paste(predictors, collapse = "*"))), data = df)
    m_df <- as.data.frame(summary(m2)$coefficients)
    m_df$comp <- "count"
    m_df$model <- "nb"
    m_df$coef <- rownames(m_df)
    rownames(m_df) <- NULL
    m_df$log2_mu <- unlist(lapply(vals, mean))
    m_df$log2FC <- unlist(lapply(seq_along(vals), function(i) {
      mean(vals[[i]]) - mean(vals[[1]])
    }))
  }
  m_df$col <- col
  return(m_df)
}
