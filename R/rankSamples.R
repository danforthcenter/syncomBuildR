#' Function to pick samples for limited dilution/other bench science following changepoint
#' regression of an ASV table or network.
#'
#' @param x ASV table as a data.frame or a \code{thresh} object. For a data.frame the cols argument
#' must be specified. For a \code{thresh} object all ASVs that were significantly related to a
#' phenotype are used (significance determined by the cutoff argument).
#' @param ... Additional arguments passed to methods.
#'
#' @return A data frame of samples ranked by the impact of their nodes (ASVs) or clusters
#'
#' @examples
#'
#' set.seed(123)
#' asv$biomass_z <- rnorm(nrow(asv))
#' x <- thresh(asv, "biomass_z")
#' rankSamples(x, cutoff = 0.5)
#'
#' @export

rankSamples <- function(x, ...) {
  UseMethod("rankSamples")
}

#' @rdname rankSamples
#' @param cutoff P value cutoff for a thresh model being significant
#' @param asvTab Optional asv table. If missing then the data slot of the thresh object is used.
#' @method rankSamples thresh
#' @export

rankSamples.thresh <- function(x, cutoff = 0.05, asvTab = NULL, ...) {
  if (is.null(asvTab)) {
    asvTab <- x$data
  }
  #* find significant pvalues from thresh
  sig_pvals <- which(unlist(lapply(x$slope, function(o) o$padj)) < cutoff)
  if (length(sig_pvals) < 1) {
    stop("No models had p-values below p.cutoff")
  }
  #* pull which ASVs those are from
  sig_asvs <- x$predictor[sig_pvals]
  #* see how many times each ASV appears (make weights)
  tab <- as.data.frame(table(sig_asvs))
  colnames(tab)[1] <- "asv"
  #* sum weighted counts per each row of ASV data
  scored <- do.call(rbind, lapply(seq_len(nrow(asvTab)), function(i) {
    .rank_samples_summing(asvTab[i, ], tab_df = tab)
  })
  )
  #* arrange by score
  ranked <- scored[order(scored$score, decreasing = TRUE), ]
  ranked$rank <- seq_len(nrow(ranked))
  return(ranked)
}

#' @rdname rankSamples
#' @param cols Columns to use in scoring. These should be some subset of ASV columns.
#' @method rankSamples data.frame
#' @export

rankSamples.data.frame <- function(x, cols = NULL, ...) {
  tab <- as.data.frame(table(cols))
  colnames(tab)[1] <- "asv"
  #* sum weighted counts per each row of ASV data
  scored <- do.call(rbind, lapply(seq_len(nrow(x)), function(i) {
    .rank_samples_summing(x[i, ], tab_df = tab)
  })
  )
  #* arrange by score
  ranked <- scored[order(scored$score, decreasing = TRUE), ]
  ranked$rank <- seq_len(nrow(ranked))
  return(ranked)
}

#' Helper function for scoring samples
#' @param r a row of an asv table
#' @param tab_df a data.frame made from a table of ASV labels with columns "asv" and "Freq"
#' @noRd
#' @keywords internal
.rank_samples_summing <- function(r, tab_df) {
  r2 <- r[, grepl("ASV", colnames(r))]
  rt <- as.data.frame(t(r2))
  colnames(rt) <- "count"
  rt$asv <- rownames(rt)
  j <- plyr::join(tab_df, rt, by = "asv")
  j$score <- j$Freq * j$count
  score <- sum(j$score)
  r$score <- score
  out <- r[, !grepl("ASV", colnames(r))]
  return(out)
}
