#' Function to demultiplex sequences either to single files or to multiple files for use with dada2.
#'
#' @param fwd_reads A ShortReadQ object or a file path to a fastq file that will be read in with
#' \code{ShortRead::readFastq}.
#' @param rev_reads A ShortReadQ object or a file path to a fastq file that will be read in with
#' \code{ShortRead::readFastq}.
#' @param barcodes A dataframe of barcodes and metadata for naming or a file path to such a flat file
#' readable by \code{utils::read.delim}.
#' @param fwd The column name of the barcodes dataframe representing forward barcodes.
#' Defaults to "FWD".
#' @param rev The column name of the barcodes dataframe representing reverse barcodes.
#' Defaults to "REV".
#' @param name One or more column names from the barcodes dataframe to be used in naming samples.
#' @param mode One of "dada" or "general". This controls how files are written out.
#' If mode is "dada" then a fastq file will be written out for
#' each row of the barcodes dataframe and the files will be named using the name argument.
#' If mode is "general" then each input fastq file will yield one
#' output fastq file with sample names in that file corresponding to the name argument.
#' Defaults to "dada".
#' @param cores Optionally the number of cores to run in parallel. This defaults to 1 if the "mc.cores"
#' option is not set.
#' @param writeOut Logical, should fastq/stats files be written out. Defaults to TRUE.
#' This can also be a path where data should be written to. Note that files are written as mode "w"
#' and this will throw an error if the file already exists.
#' @param stat Logical, should stats files be generated if writeOut is TRUE. Defaults to FALSE.
#'
#'
#' @keywords demultiplex, dada2
#'
#' @importFrom utils write.csv
#' @importFrom methods is
#' @importClassesFrom ShortRead ShortReadQ
#' @importMethodsFrom ShortRead append
#' @importFrom ShortRead readFastq writeFastq sread countFastq
#' @import parallel
#'
#' @return Writes out fq files and summary statistics depending on writeOut and stat arguments.
#' If both are TRUE then summary stats are returned as a dataframe.
#'
#' @export

demultiplex <- function(fwd_reads, rev_reads, barcodes, fwd = "FWD", rev = "REV",
                        name = c("Name", "Sample_name"), mode = c("dada", "general"),
                        cores = getOption("mc.cores", 1), writeOut = TRUE, stat = FALSE) {
  if (is.character(fwd_reads)) {
    fwd_reads <- ShortRead::readFastq(fwd_reads)
  }
  if (is.character(rev_reads)) {
    rev_reads <- ShortRead::readFastq(rev_reads)
  }
  if (!methods::is(fwd_reads, "ShortReadQ") || !methods::is(rev_reads, "ShortReadQ")) {
    stop("fwd and rev reads must be a ShortReadQ object or a file path to a fastq file.")
  }
  if (is.character(barcodes)) {
    barcodes <- utils::read.delim(barcodes)
  }
  if (!methods::is(barcodes, "data.frame")) {
    stop("barcodes must be a data.frame or a file path to a file readable by read.delim.")
  }
  if (!all(c(fwd, rev, name) %in% colnames(barcodes))) {
    stop("fwd, rev, and name must all be column names in the barcodes dataframe.")
  }
  mode <- mode[1]
  if (mode == "dada") {
    out <- .demultiplex_for_dada(fwd_reads, rev_reads, barcodes, fwd, rev, name, writeOut, stat, cores)
  } else if (mode == "general") {
    out <- .demultiplex_general(fwd_reads, rev_reads, barcodes, fwd, rev, name, writeOut, stat, cores)
  } else {
    stop("mode must be one of 'dada' or 'general'")
  }
  if (!is.logical(writeOut)) {
    writeOut <- TRUE
  }
  if (stat && writeOut) {
    out
  }
}

#' ***********************************************************************************************
#' *************** `DADA2 method` ****************************************
#' ***********************************************************************************************
#'
#' @description
#' Internal function for demultiplexion to the style we have used for dada2
#' @keywords internal
#' @noRd

.demultiplex_for_dada <- function(fwd_reads, rev_reads, barcodes, fwd = "FWD", rev = "REV",
                                  name = c("Name", "Sample_name"), writeOut = TRUE,
                                  stat = FALSE, cores = 1) {
  # format barcodes to uppercase
  barcodes[[fwd]] <- toupper(barcodes[[fwd]])
  barcodes[[rev]] <- toupper(barcodes[[rev]])
  fwd_reads_char <- as.character(ShortRead::sread(fwd_reads))
  rev_reads_char <- as.character(ShortRead::sread(rev_reads))

  test <- parallel::mclapply(seq_len(nrow(barcodes)), function(i) {
    # get a name for the sub sample
    file_name <- paste(barcodes[i, c(name)], collapse = "_")
    # get barcodes
    fwd_bar <- barcodes[i, fwd]
    rev_bar <- barcodes[i, rev]

    fwd_reads_w_fwd_bars <- which(grepl(paste0("^", fwd_bar), fwd_reads_char))
    rev_reads_w_rev_bars <- which(grepl(paste0("^", rev_bar), rev_reads_char))
    fwd_index_both <- intersect(fwd_reads_w_fwd_bars, rev_reads_w_rev_bars) # logical AND

    fwd_reads_w_rev_bars <- which(grepl(paste0("^", rev_bar), fwd_reads_char))
    rev_reads_w_fwd_bars <- which(grepl(paste0("^", fwd_bar), rev_reads_char))
    rev_index_both <- intersect(fwd_reads_w_rev_bars, rev_reads_w_fwd_bars) # logical AND
    # union of both groups of reads
    total_index <- union(fwd_index_both, rev_index_both) # logical OR
    fwd_clone <- fwd_reads[total_index]
    rev_clone <- rev_reads[total_index]
    # write fastq
    if (!is.logical(writeOut)) {
      write <- writeOut
      if (substr(write, nchar(write), nchar(write)) != "/") {
        write <- paste0(write, "/")
      }
      writeOut <- TRUE
    } else {
      write <- NULL
    }
    if (writeOut) {
      ShortRead::writeFastq(
        object = fwd_clone, file = paste0(write, file_name, "_f.fastq.gz"),
        mode = "w", compress = TRUE
      )
      ShortRead::writeFastq(
        object = rev_clone, file = paste0(write, file_name, "_r.fastq.gz"),
        mode = "w", compress = TRUE
      )
    }
    if (stat & writeOut) {
      stats <- data.frame(
        name = file_name,
        totalReads = length(fwd_reads),
        groupedReads_fwd = ShortRead::countFastq(paste0(write, file_name, "_f.fastq.gz"))[1, 1],
        groupedReads_rev = ShortRead::countFastq(paste0(write, file_name, "_r.fastq.gz"))[1, 1]
      )
      stats$ratio <- signif(stats$groupedReads_fwd / stats$totalReads, digits = 4)
      utils::write.csv(stats, file = paste0(write, file_name, ".csv"), row.names = FALSE)
      return(stats)
    }
  }, mc.cores = cores)
  if (!is.logical(writeOut)) {
    writeOut <- TRUE
  }
  if (stat && writeOut) {
    return(do.call(rbind, test))
  }
}

#' ***********************************************************************************************
#' *************** `Generalized method` ****************************************
#' ***********************************************************************************************
#'
#' @description
#' Internal function for demultiplexion to the style previously used in U/VSEARCH by mingsheng
#' @keywords internal
#' @noRd

.demultiplex_general <- function(fwd_reads, rev_reads, barcodes, fwd = "FWD", rev = "REV",
                                 name = c("Name", "Sample_name"), writeOut = TRUE,
                                 stat = FALSE, cores = 1) {
  # placeholder until I change how general works
  return(.demultiplex_for_dada(fwd_reads, rev_reads, barcodes, fwd, rev,
                               name, writeOut,
                               stat, cores))
  # format barcodes to uppercase
  barcodes[[fwd]] <- toupper(barcodes[[fwd]])
  barcodes[[rev]] <- toupper(barcodes[[rev]])
  # per barcode write out a fq and optionally a summary table
  lapply(seq_len(nrow(barcodes)), function(i) {
    # get a name for the sub sample
    name <- paste(barcodes[i, c(name)], collapse = "_")
    # get barcodes
    fwd_bar <- barcodes[i, fwd]
    rev_bar <- barcodes[i, rev]
    # filter for barcodes at beginning of sequence
    rowReads <- reads[grepl(paste0("^", fwd_bar), ShortRead::sread(reads))]
    columnReads <- reads[grepl(paste0("^", rev_bar), ShortRead::sread(reads))]
    # make reverse complements of fwd and reverse barcodes
    fwd_bar_reverse_complement <- as.character(
      Biostrings::reverseComplement(
        Biostrings::DNAString(fwd_bar)
      )
    )
    rev_bar_reverse_complement <- as.character(
      Biostrings::reverseComplement(
        Biostrings::DNAString(rev_bar)
      )
    )
    # search for reverse complements of opposite barcode at end of sequences
    rowReads <- rowReads[grepl(paste0(rev_bar_reverse_complement, "$"), ShortRead::sread(rowReads))]
    columnReads <- columnReads[grepl(
      paste0(fwd_bar_reverse_complement, "$"),
      ShortRead::sread(columnReads)
    )]
    # append row and column reads
    clone <- append(rowReads, columnReads)
    # write fastq
    if (!is.logical(writeOut)) {
      write <- writeOut
      if (substr(write, nchar(write), nchar(write)) != "/") {
        write <- paste0(write, "/")
      }
      writeOut <- TRUE
    } else {
      write <- NULL
    }
    if (writeOut) {
      if (i == 1) {
        md <- "w"
      } else {
        md <- "a"
      }
      ShortRead::writeFastq(
        object = clone, file = paste0(write, name, ".fastq.gz"),
        mode = md, compress = TRUE
      )
    }
  })
  if (stat && writeOut) {
    stats <- data.frame(
      name = name,
      totalReads = length(reads),
      groupedReads = ShortRead::countFastq(paste0(name, ".fastq.gz"))[1, 1]
    )
    stats$ratio <- signif(stats$groupedReads / stats$totalReads, digits = 4)
    utils::write.csv(stats, file = paste0(name, ".csv"), row.names = FALSE)
    return(stats)
  }
}
