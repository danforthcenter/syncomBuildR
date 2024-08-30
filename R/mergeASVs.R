#' Function combine ASV tables that have been aligned with blastASVs
#'
#' @param blast A dataframe as returned by \link{blastASVs}.
#' @param asv A list of two ASV tables.
#' @param taxa Optional taxonomy information in the form returned by dada2 for both experiments.
#' Defaults to NULL.
#' @param shared_meta Optional column names present in both ASV tables that you wish to include in the
#' new ASV table output.
#' @param check_taxa Logical, should taxa be checked for overlap and the results included in the link
#' part of the output? Note even if this is false you will be able to tell if the taxa are misaligned
#' in the taxa part of the output if taxa is non NULL,
#' @param min_pident Minimum pident match returned by blastASVs to merge. Defaults to 99.
#' @param cores Optionally a number of cores to run in parallel using parallel::mclapply. Defaults to 1
#' if the "mc.cores" option is unset
#' @keywords combine ASVs
#' @importFrom parallel mclapply
#' @importFrom utils tail
#' @import stats
#' @return A named list of dataframes including "asv" and "links", optionally with "taxa"
#'
#' @examples
#' \dontrun{
#' res <- read.csv("blast_2022_vs_2021.csv")
#' print(load("~/scripts/SINC/field_2021/microbiome/sinc_field_2021_asvTable.rdata"))
#' asv21 <- asv
#' taxa21 <- taxa
#' colnames(asv21)[c(1, 7)] <- c("genotype", "tissue")
#' print(load("~/scripts/SINC/field_2022/microbiome/2022_asvTable.rdata"))
#'
#' blast <- res
#' asv <- list("d1" = asv22, "d2" = asv21)
#' taxa <- list("d1" = taxa22, "d2" = taxa21)
#' shared_meta <- c("genotype", "tissue") # shared metadata to keep from both tables
#'
#' ex <- mergeASVs(blast, asv, taxa, shared_meta)
#' }
#' @export

mergeASVs <- function(blast, asv, taxa, shared_meta, check_taxa = TRUE,
                      min_pident = 99.5, cores = getOption("mc.cores", 1)) {
  #* label which experiment things are from in a dataframe
  new_asv_metadata <- data.frame("experiment" = rep(c("d1", "d2"), c(nrow(asv$d1), nrow(asv$d2))))
  #* assemble other metadata and bind it to the new_asv_metadata
  if (!is.null(shared_meta)) {
    shared_metadata <- stats::setNames(as.data.frame(do.call(cbind, lapply(shared_meta, function(col) {
      unlist(lapply(asv, function(d) {
        d[, col]
      }))
    }))), c(shared_meta))
    new_asv_metadata <- cbind(new_asv_metadata, shared_metadata)
  }
  #* split the blastASVs results by sequence in table 1
  dList <- split(blast, blast$seq1_id)
  #* Main looping function that will run on each ASV in parallel according to cores/"mc.cores" option
  #* Begin main loop
  linkage_data <- parallel::mclapply(seq_along(dList)[1:10], function(i) {
    d <- dList[[i]]
    #* in case there is no match return original counts from first ASV table padded with NAs
    asvColumn <- setNames(data.frame(
      c(asv[[1]][[d$seq1_id[1]]], rep(NA, nrow(asv[[2]])))
    ), c(paste0("ASV", i)))
    #* get link between ASVs, picking first match and filtering for pident in case you had
    #* different options in blastASVs
    link <- d[which(d$match_rank == min(d$match_rank, na.rm = TRUE)), ]
    link <- link[link$pident >= min_pident, ]
    #* if there is a match then do more stuff and make an ASV Column of counts with a new ASV number
    if (nrow(link) >= 1) {
      link$used <- FALSE
      link[1, "used"] <- TRUE
      ASV1 <- link[link$used, "seq1_id"]
      ASV2 <- link[link$used, "db_id"]
      asvColumn <- stats::setNames(
        data.frame(c(asv[[1]][[ASV1]], asv[[2]][[ASV2]])),
        c(paste0("ASV", i))
      )
      #* record new ASV number
      link$new_ASV_number <- paste0("ASV", i)
      #* if there was a link and we are checking taxa then check taxa here and add to the link data.
      if (check_taxa & !is.null(taxa)) {
        tx1 <- taxa[[1]][ASV1, ]
        tx2 <- taxa[[2]][ASV2, ]
        #* record how low their agreement gets and whether they disagree on any non NA fields
        taxaRes <- unlist(lapply(colnames(tx1), function(taxaLevel) {
          if (any(is.na(c(tx1[[taxaLevel]], tx2[[taxaLevel]])))) {
            NA
          } else if (tx1[[taxaLevel]] == tx2[[taxaLevel]]) {
            paste0("agree to ", taxaLevel)
          } else {
            paste0("disagree on ", taxaLevel)
          }
        }))
        link$agreement <- utils::tail(taxaRes[grepl("agree to", taxaRes)], n = 1)
        link$disagreement <- utils::tail(taxaRes[grepl("disagree on", taxaRes)], n = 1)
      }
      #* make list to be returned.
      outList <- list(link = link, asvColumn = asvColumn)
      #* if the taxa information is provided then make a new taxa table.
      if (!is.null(taxa)) {
        #* make crossed taxa table, POTENTIALLY with 2 rows per ASV if there is a problem.
        #* if that happens it will be flagged later and a warning is returned too.
        tx1 <- taxa[[1]][ASV1, ]
        tx2 <- taxa[[2]][ASV2, ]
        new_taxa_iter <- rbind(tx1, tx2)
        rownames(new_taxa_iter) <- NULL
        new_taxa_iter <- new_taxa_iter[!duplicated(new_taxa_iter), ]
        new_taxa_iter$ASVnumber <- paste0("ASV", i)
        outList$taxa <- new_taxa_iter
      }
    } else { # if there is no link
      #* make list to be returned in case there is no link
      asvColumn <- setNames(data.frame(
        c(asv[[1]][[d$seq1_id[1]]], rep(NA, nrow(asv[[2]])))
      ), c(paste0("ASV", i)))
      link <- data.frame(
        seq1_id = d$seq1_id[1], db_id = NA, pident = NA,
        length = NA, mismatch = NA, gapopen = NA, qstart = NA,
        qend = NA, sstart = NA, send = NA, evalue = NA, bitscore = NA,
        match_rank = NA, used = TRUE, new_ASV_number = paste0("ASV", i)
      )
      if (check_taxa & !is.null(taxa)) {
        link$agreement <- NA
        link$disagreement <- NA
      }
      outList <- list(link = link, asvColumn = asvColumn)
      if (!is.null(taxa)) {
        new_taxa_iter <- taxa[[1]][d$seq1_id[1], ]
        new_taxa_iter$ASVnumber <- paste0("ASV", i)
        outList$taxa <- new_taxa_iter
      }
    }
    return(outList)
  }, mc.cores = cores) # end main loop
  #* bind results from the main loop part of the function
  links <- do.call(rbind, lapply(linkage_data, function(lst) {
    lst$link
  }))
  new_asv_counts <- do.call(cbind, lapply(linkage_data, function(lst) {
    lst$asvColumn
  }))
  #* add metadata to new_asv_counts
  new_asv_table <- cbind(new_asv_metadata, new_asv_counts)
  #* Final output is a list, initializing it here
  final_output <- list(
    "asv" = new_asv_table,
    "links" = links
  )
  if (!is.null(taxa)) {
    new_taxa <- do.call(rbind, lapply(linkage_data, function(lst) {
      lst$taxa
    }))
    #* check for bad taxa matches
    if (any(as.numeric(table(new_taxa$ASVnumber)) > 1)) {
      warning(paste0(
        "Some non-unique taxa information, ",
        "check taxonomy results carefully for duplicated ASVnumbers"
      ))
    }
    #* add to final_output list
    final_output$taxa <- new_taxa
  }
  #* find columns from second table that are not included already
  #* this would happen whenever there are some matches below the pident cutoff
  n1 <- ncol(final_output$new_asv_counts)
  remaining_columns <- setdiff(colnames(asv[[2]]), links$db_id)
  if (length(remaining_columns) > 0) {
    #* need to make another set of ASV counts, grab taxa information, and make links.
    pt2_asv_counts <- do.call(cbind, lapply(seq_along(remaining_columns), function(i) {
      col <- remaining_columns[i]
      n2 <- n1 + i
      stats::setNames(
        data.frame(c(rep(NA, nrow(asv[[1]])), asv[[2]][[col]])),
        paste0("ASV", n2)
      )
    }))
    final_output$asv <- cbind(final_output$asv, pt2_asv_counts)
    pt2_link <- data.frame(
      seq1_id = NA, db_id = remaining_columns, pident = NA,
      length = NA, mismatch = NA, gapopen = NA, qstart = NA,
      qend = NA, sstart = NA, send = NA, evalue = NA, bitscore = NA,
      match_rank = NA, used = TRUE,
      new_ASV_number = paste0(
        "ASV",
        (n1 + 1):(n1 + length(remaining_columns))
      )
    )
    if (check_taxa && !is.null(taxa)) {
      pt2_link$agreement <- NA
      pt2_link$disagreement <- NA
    }
    final_output$links <- cbind(final_output$links, pt2_link)
    if (!is.null(taxa)) {
      pt2_taxa <- taxa[[2]][remaining_columns, ]
      pt2_taxa$ASVnumber <- remaining_columns
      final_output$taxa <- rbind(final_output$taxa, pt2_taxa)
    }
  }
  return(final_output)
}
