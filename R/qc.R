#' Ease of use quality control function for DADA2 output.
#'
#' @param file A file path to an Rdata object containing DADA2 output in such as that from DADA2
#' examples (function?)
#' @param asvTab ASV table from DADA2. If a file path is given then this defaults to seqtab.print if
#' left NULL for consistency with common DADA2 examples and previous SINC work.
#' @param taxa A taxonomy file for the ASVs in asvTab. This defaults to taxa_rdp if left NULL.
#' @param asvAbnd Lower bound on sum per microbe across all samples, ASVs with fewer total counts than
#' this number will be removed.
#' @param sampleAbnd Lower bound on sum of all microbes per sample, samples with fewer total counts
#' than this number will be removed. Note this filtering happens before `asvAbnd` filtering.
#' @param normalize Method for normalization, currently supported options are NULL
#' (where no normalization is done),
#' "rescale" (the default) where samples are rescaled to have the same sample abundance using the
#' exponent of the mean log sample sum,
#'  and "pqn" where probabilistic quotient normalization is used.
#' @param rmTx A list of taxa filters. Each element in the list should specify a level of the taxonomy,
#' an affector, and a comma separated string of values to filter that taxonomy for. See details.
#' @param separate An optional dataframe or vector of the same length as the number of samples in
#' asvTab. If provided then this separates the asvTab by these groupings and abundance filtering will
#' be done separately in each subset of the data. Subsets will be bound together and missing values
#' will be filled with NA should subsets contain different selections of ASVs.
#' @param metadata An optional dataframe of metadata to add to the asv table. If left NULL while
#' `separate` is a dataframe then this will include that data as metadata. This must include 1 row
#' per row of your asv table.
#' @param return_removed Logical, abundance filtered data be returned? Defaults to FALSE. If TRUE then
#' a list is returned.
#' @keywords DADA2
#' @importFrom plyr rbind.fill
#' @return An ASV table as a wide dataframe
#' @details
#'
#' The `rmTx` argument specifies a list of filters, with one element per each level of taxonomy that
#' should be filtered. Each element in the list should take the form of "Taxalevel in example, example2".
#' Which will remove ASVs where their "Taxalevel" is "example" or "example2". Alternatively the affector
#' can be "contains" in which case regex is used.
#' In that case the string "Class contains plant, other" would search for matches to "plant|other" in
#' the Class column of the taxonomy. More general regex patterns are supported.
#'
#'
#' @examples
#'
#' print(load("~/scripts/SINC/sincUtils/syncomBuilder/field_2021_small_ex.rdata"))
#' separate <- data.frame(tissue = stringr::str_extract(samps, "[:alpha:]+"))
#' metadata <- cbind(separate, data.frame(
#'   sample = samps,
#'   plot = stringr::str_extract(samps, "[0-9]+")
#' ))
#'
#'
#' asv <- qc(
#'   asvTab = seqtab.print, taxa = taxa_rdp, asvAbnd = 100, sampleAbnd = 1000, normalize = "rescale",
#'   rmTx = list("Order in Chloroplast", "Phylum in Plantae"), separate = NULL
#' )
#' asv <- qc(
#'   asvTab = seqtab.print, taxa = taxa_rdp, asvAbnd = 100, sampleAbnd = 1000, normalize = "rescale",
#'   rmTx = list("Order in Chloroplast", "Phylum in Plantae"), separate = separate, metadata = metadata,
#'   return_removed = TRUE
#' )
#' save(asv, file = "../qc_output.rdata")
#'
#' @export

qc <- function(file = NULL, asvTab = NULL, taxa = NULL, asvAbnd = 100, sampleAbnd = 1000,
               normalize = "rescale",
               rmTx = list("Order in Chloroplast", "Phylum in Plantae", "Kingdom in Eukaryota"),
               separate = NULL, metadata = NULL, return_removed = FALSE) {
  if (!is.null(file)) {
    seqtab.print <- NULL
    taxa_rdp <- NULL
    load(file)
    asvTab <- seqtab.print
    taxa <- taxa_rdp
  }
  asvTab <- as.data.frame(asvTab)
  taxa <- as.data.frame(taxa)
  split <- !is.null(separate)
  if (is.null(metadata) && is.data.frame(separate)) {
    metadata <- separate
  }

  if (!is.null(rmTx)) {
    indsL <- lapply(rmTx, function(stri) {
      taxaLevel <- strsplit(stri, " ")[[1]][1]
      affector <- strsplit(stri, " ")[[1]][2]
      values <- trimws(gsub(",", "", strsplit(stri, " ")[[1]][-c(1:2)]))

      if (affector %in% c("in", "is", "=")) {
        indL <- lapply(values, function(value) {
          index <- taxa[[taxaLevel]] != value # flag F if the row in taxa should be removed
          index[is.na(index)] <- TRUE # if taxaLevel[i] is NA then keep it
          index
        })
        ind <- apply(matrix(unlist(indL), ncol = length(indL)), MARGIN = 1, FUN = all)
      } else if (affector == "contains") {
        valReg <- paste0(values, collapse = "|")
        ind <- !grepl(valReg, taxa[[taxaLevel]])
      }
      ind
    })
    keep_index <- apply(matrix(unlist(indsL), ncol = length(indsL)), MARGIN = 1, FUN = all)
    taxa <- taxa[keep_index, ]
    asvTab <- asvTab[, keep_index]
  }

  if (!is.null(sampleAbnd)) {
    originalRowNames <- rownames(asvTab)
    if (!split) {
      asvTabSampFilt <- asvTab[rowSums(asvTab) >= sampleAbnd, ]
    } else {
      asvTabSampFilt <- do.call(
        rbind,
        lapply(split(asvTab, interaction(separate, drop = TRUE)), function(d) {
          d[rowSums(d) >= sampleAbnd, ]
        })
      )
      rownames(asvTabSampFilt) <- sub(
        paste0(
          levels(interaction(separate, drop = TRUE)),
          ".",
          collapse = "|"
        ),
        "",
        rownames(asvTabSampFilt)
      )
      sepColNames <- colnames(separate)
      separate <- as.data.frame(separate[which(rownames(asvTabSampFilt) %in% originalRowNames), ])
      colnames(separate) <- sepColNames
    }
    if (!is.null(metadata)) {
      metadata <- metadata[rowSums(asvTab) >= sampleAbnd, ]
    }
    if (return_removed) {
      removed <- list()
      removed$rows <- asvTab[rowSums(asvTab) < sampleAbnd, ]
    }
  } else {
    asvTabSampFilt <- asvTab
  }
  if (!is.null(asvAbnd)) {
    if (!split) {
      asvTabFilt <- asvTabSampFilt[, colSums(asvTabSampFilt) >= asvAbnd]
    } else {
      asvTabFilt <- do.call(
        plyr::rbind.fill,
        lapply(split(asvTabSampFilt, interaction(separate, drop = TRUE)), function(d) {
          d[, colSums(d) >= asvAbnd]
        })
      )
    }
    if (return_removed) {
      removed <- list()
      removed$cols <- asvTab[, -which(colnames(asvTabFilt) %in% colnames(asvTab))]
    }
  } else {
    asvTabFilt <- asvTabSampFilt
  }
  if (!is.null(normalize)) {
    if (normalize == "rescale") {
      scaleFactor <- exp(mean(log(rowSums(asvTabFilt, na.rm = TRUE))))
      asvTabFilt <- as.data.frame(
        t(apply(asvTabFilt, MARGIN = 1, FUN = function(i) {
          round(scaleFactor * i / sum(i, na.rm = TRUE))
        }))
      )
    } else if (normalize == "pqn") {
      median_vec <- apply(asvTabFilt / rowMeans(asvTab, na.rm = TRUE), 2, median, na.rm = TRUE) + 1
      normalizedAsvTab <- t(apply(asvTabFilt, 1, function(x) x / median_vec))
      rownames(normalizedAsvTab) <- rownames(asvTabFilt)
      asvTab <- normalizedAsvTab
    } else {
      warning(paste0(
        "Available options for normalization are NULL, 'pqn', and 'rescale', ",
        normalize, " is not implemented"
      ))
    }
  }

  if (!is.null(metadata)) {
    asvTab <- cbind(metadata, asvTabFilt)
  }
  if (return_removed) {
    out <- list("asv" = asvTab, "removed" = removed)
  } else {
    out <- asvTab
  }
  return(out)
}
