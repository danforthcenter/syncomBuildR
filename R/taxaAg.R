#' Make named vector of taxonomy sums
#'
#' @param x ASV table, asv columns thought start with "ASV" prefix.
#' @param tx Taxonomy table in dada2 format. Rownames should correspond to ASV names in the same format as x's column names.
#' @param taxaLevel A level of the taxonomy in `tx` to aggregate data to. Must match a colname of `tx`.
#'
#' @keywords taxonomy, dirichlet, conjugate
#' @details
#' This function can be used to make counts for use in a dirichlet updating function such as `diriUpdate`
#'
#' @export

taxaSums <- function(x = NULL, tx = NULL, taxaLevel = "Class") {
  asvList <- colnames(x)[grepl("ASV", colnames(x))] # take the asv names from both samples
  tx$asv <- rownames(tx) # put the asv names in the taxa data
  x_sum <- data.frame(asv = asvList, sum = colSums(x[, colnames(x) %in% asvList])) # make a long version of the data with a columns for sample sum and ASV number
  x_sum$taxa <- unlist(lapply(x_sum$asv, function(a) tx[tx$asv == a, taxaLevel])) # label taxonomy
  x_sum$taxa <- replace(x_sum$taxa, is.na(x_sum$taxa), "Not_Assigned") # replace taxa NAs
  ag <- aggregate(x_sum, sum ~ taxa, FUN = sum) # sum the counts by taxonomy
  out <- as.numeric(ag$sum) # transpose the asv data
  names(out) <- ag$taxa # add names for taxa
  return(out)
}

#' Aggregate an ASV table to some taxonomic level
#'
#' @param x ASV table, asv columns thought start with "ASV" prefix.
#' @param tx Taxonomy table in dada2 format. Rownames should correspond to ASV names in the same format as x's column names.
#' @param taxaLevel A level of the taxonomy in `tx` to aggregate data to. Must match a colname of `tx`.
#'
#' @return A list of two dataframes. "data" is similar to an ASV table but aggregated to the given taxaLevel.
#' "meta" is a summary of the aggregation containing the number of ASVs and the total number of counts in each family in long format.
#'
#' @keywords network dirichlet
#'
#' @export

taxaAg <- function(x = NULL, tx = NULL, taxaLevel = "Class") {
  asvList <- colnames(x)[grepl("ASV", colnames(x))] # take the asv names from both samples
  tx$asv <- rownames(tx) # put the asv names in the taxa data
  tx[[taxaLevel]] <- replace(tx[[taxaLevel]], is.na(tx[[taxaLevel]]), "Not_Assigned")
  agList <- lapply(unique(tx[[taxaLevel]]), function(group) {
    subtx <- tx[tx[[taxaLevel]] == group, "asv"]
    if (length(subtx) > 1) {
      ag <- stats::setNames(data.frame(rowSums(x[, subtx])), group)
    } else {
      ag <- stats::setNames(data.frame(x[, subtx]), group)
    }
    meta <- data.frame(group = group, taxaLevel = taxaLevel, nAsvs = length(subtx), sumAsvs = sum(ag[[1]]))
    return(list("aggregate" = ag, "meta" = meta))
  })
  ag <- do.call(cbind, lapply(agList, function(a) {
    a$aggregate
  }))
  meta <- do.call(rbind, lapply(agList, function(a) {
    a$meta
  }))
  out_df <- cbind(x[, !grepl("ASV", colnames(x))], ag)
  out <- list("data" = out_df, "meta" = meta)
  return(out)
}
