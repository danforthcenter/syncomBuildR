#' Ease of use quality control function for DADA2 output.
#'
#' @param file A file path to an Rdata object containing DADA2 output formatted per sinc standards
#' (ie, following the dada2 vignette conventions).
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
#' @examples
#'
#' set.seed(123)
#' plot_numbers <- rep(seq(10, 21, 1), each = 3)
#' tissues <- c("re", "rr", "s")
#' sample_numbers <- paste0("S", 1:36)
#' samples <- paste(plot_numbers, tissues, sample_numbers, sep = "_")
#' metadata <- data.frame(
#'   plot = plot_numbers,
#'   tissue = rep(tissues, times = 12),
#'   sample = samples
#' )
#' asvs <- matrix(floor(rlnorm(360, log(1), 1.5)), ncol = 10)
#' colnames(asvs) <- paste0("ASV", 1:10)
#' taxa <- c(
#'   "Bacteria", "Proteobacteria", "Betaproteobacteria", "Burkholderiales",
#'   "Burkholderiaceae", "Paraburkholderia", NA
#' )
#' taxa <- matrix(rep(taxa, 10), nrow = 10, byrow = TRUE)
#' colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
#' rownames(taxa) <- colnames(asvs)
#' asv <- qc(
#'   asvTab = asvs, taxa = taxa, asvAbnd = 10,
#'   sampleAbnd = 10, normalize = "rescale",
#'   rmTx = list("Order in Chloroplast", "Phylum in Plantae"), separate = NULL,
#'   metadata = metadata
#' )
#'
#' @export

qc <- function(file = NULL, asvTab = NULL, taxa = NULL, asvAbnd = 100, sampleAbnd = 1000,
               normalize = "rescale",
               rmTx = list("Order in Chloroplast", "Phylum in Plantae", "Kingdom in Eukaryota"),
               separate = NULL, metadata = NULL, return_removed = FALSE) {
  file_helper_output <- .qc_file_helper(file, asvTab, taxa)
  asvTab <- as.data.frame(file_helper_output$asvTab)
  taxa <- as.data.frame(file_helper_output$taxa)
  split <- !is.null(separate)
  if (is.null(metadata) && is.data.frame(separate)) {
    metadata <- separate
  }
  #* `Taxa filtering`
  taxa_filtering_output <- .qc_taxa_filter(asvTab, taxa, rmTx)
  asvTab <- taxa_filtering_output$asvTab
  taxa <- taxa_filtering_output$taxa
  #* `Sample Abundance filtering`
  sample_filtering_output <- .qc_sample_abundance_filter(asvTab, sampleAbnd, metadata, return_removed,
                                                         separate, split)
  asvTab <- sample_filtering_output$asvTab
  metadata <- sample_filtering_output$metadata
  removed <- sample_filtering_output$removed
  #* `ASV Abundance filtering`
  asv_filtering_output <- .qc_asv_abundance_filter(asvTab, removed, split, separate)
  asvTab <- asv_filtering_output$asvTab
  removed <- asv_filtering_output$removed
  #* `Normalization`
  asvTab <- .qc_normalization(asvTab, normalize)
  #* `Return ASV table`
  if (!is.null(metadata)) {
    asvTab <- cbind(metadata, asvTab)
  }
  if (return_removed) {
    out <- list("asv" = asvTab, "removed" = removed)
  } else {
    out <- asvTab
  }
  return(out)
}

#' @keywords internal
#' @noRd

.qc_file_helper <- function(file, asvTab, taxa) {
  if (!is.null(file)) {
    seqtab.print <- NULL
    taxa_rdp <- NULL
    load(file)
    asvTab <- seqtab.print
    taxa <- taxa_rdp
  }
  return(list("asvTab" = asvTab, "taxa" = taxa))
}

#' @keywords internal
#' @noRd

.qc_taxa_filter <- function(asvTab, taxa, rmTx) {
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
  return(list("asvTab" = asvTab, "taxa" = taxa))
}

#' @keywords internal
#' @noRd

.qc_sample_abundance_filter <- function(asvTab, sampleAbnd, metadata,
                                        return_removed, separate, split) {
  removed <- NULL
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
  return(list("asvTab" = asvTabSampFilt, "removed" = removed, "metadata" = metadata))
}

#' @keywords internal
#' @noRd

.qc_asv_abundance_filter <- function(asvTab, asvAbnd, removed, split, separate) {
  if (!is.null(asvAbnd)) {
    if (!split) {
      asvTabFilt <- asvTab[, colSums(asvTab) >= asvAbnd]
    } else {
      asvTabFilt <- do.call(
        plyr::rbind.fill,
        lapply(split(asvTab, interaction(separate, drop = TRUE)), function(d) {
          d[, colSums(d) >= asvAbnd]
        })
      )
    }
    if (!is.null(removed)) {
      removed$cols <- asvTab[, -which(colnames(asvTabFilt) %in% colnames(asvTab))]
    }
  } else {
    asvTabFilt <- asvTab
  }
  return(list("asvTab" = asvTabFilt, "removed" = removed))
}

#' @keywords internal
#' @noRd

.qc_normalization <- function(asvTab, normalize) {
  if (!is.null(normalize)) {
    if (normalize == "rescale") {
      scaleFactor <- exp(mean(log(rowSums(asvTab, na.rm = TRUE))))
      normalizedAsvTab <- as.data.frame(
        t(apply(asvTab, MARGIN = 1, FUN = function(i) {
          round(scaleFactor * i / sum(i, na.rm = TRUE))
        }))
      )
    } else if (normalize == "pqn") {
      median_vec <- apply(asvTab / rowMeans(asvTab, na.rm = TRUE), 2, median, na.rm = TRUE) + 1
      normalizedAsvTab <- t(apply(asvTab, 1, function(x) x / median_vec))
      rownames(normalizedAsvTab) <- rownames(asvTab)
    } else {
      stop(paste0(
        "Available options for normalization are NULL, 'pqn', and 'rescale', ",
        normalize, " is not implemented"
      ))
    }
  } else {
    normalizedAsvTab <- asvTab
  }
  return(normalizedAsvTab)
}
