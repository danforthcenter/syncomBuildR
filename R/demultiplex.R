#' 
#' Function to demultiplex sequences either to single files or to multiple files for use with dada2.
#' 
#' @param reads A ShortReadQ object or a file path to a fastq file that will be read in with \code{ShortRead::readFastq}.
#' @param barcodes A dataframe of barcodes and metadata for naming or a file path to such a flat file readable by \code{utils::read.delim}.
#' @param fwd The column name of the barcodes dataframe representing forward barcodes. Defaults to "FWD".
#' @param rev The column name of the barcodes dataframe representing reverse barcodes. Defaults to "REV".
#' @param name One or more column names from the barcodes dataframe to be used in naming samples.
#' @param mode One of "dada" or "general". This controls how files are written out. If mode is "dada" then a fastq file will be written out for 
#' each row of the barcodes dataframe and the files will be named using the name argument. If mode is "general" then each input fastq file will yield one
#' output fastq file with sample names in that file corresponding to the name argument. Defaults to "dada".
#' @param cores Optionally the number of cores to run in parallel. This defaults to 1 if the "mc.cores" option is not set.
#' @param writeOut Logical, should fastq/stats files be written out. Defaults to TRUE.
#' This can also be a path where data should be written to. Note that files are written as mode "w" and this will 
#' throw an error if the file already exists.
#' @param stat Logical, should stats files be generated if writeOut is TRUE. Defaults to FALSE.
#' 
#' 
#' @keywords demultiplex, dada2
#' 
#' @importFrom ShortRead readFastq sread writeFastq countFastq 
#' @importFrom Biostrings reverseComplement DNAString
#' @import parallel
#' 
#' @return Writes out fq files and summary statistics depending on writeOut and stat arguments.
#' If both are TRUE then summary stats are returned as a dataframe.
#' 
#' @examples 
#' 
#' ## Not Run:
#' if(FALSE){
#' 
#' library(ShortRead)
#' library(Biostrings)
#' library(parallel)
#' setwd("~/Desktop/stargate/SINC/sincUtils/syncomBuilder/nastya96WellEx")
#' 
#' barcodes <- read.delim("barcode_tab.tsv")
#' name = c("Name", "Well")
#' reads <- ShortRead::readFastq("Plate11_R1_001.fastq.gz")
#' 
#' demultiplex(reads, barcodes[1:10,], fwd="FWD", rev="REV",
#' name = name, mode = "dada", cores = 10, writeOut = "tests", stat = FALSE)
#' 
#' path <- '/run/user/1000/gvfs/smb-share:server=nas01,share=research/bart_lab/Previous_people/Mingsheng_Qi/research_dr/SINC/limiting_dilution/20210314/p2p11_rawseq'
#' barcodes = read.csv(paste0(path,"/barcodetable1.csv"))
#' barcodes$exp <- "exp"
#' reads <- ShortRead::readFastq(paste0(path, "/210314S05P10_merge.fq"))
#' 
#' x<-demultiplex(reads, barcodes, fwd="FWD_primer", rev="REV_primer",name = c("exp"), mode = "general",
#'                writeOut = TRUE, stat=TRUE, cores=10)
#' x
#' 
#' 
#' # benchmarking against:
#' # "plate","totalReads","groupedReads","ratio"
#' # "P10",234509,154507,0.6589
#' }
#' ## End not run
#' 
#' @export


demultiplex <- function(reads, barcodes, fwd="FWD", rev="REV",name = c("Name", "Sample_name"),
                        mode=c("dada", "general"), cores = getOption("mc.cores",1), writeOut = TRUE, stat=FALSE){
  
  if(is.character(reads)){
    reads <- ShortRead::readFastq(reads)
  }
  if(!is(reads, "ShortReadQ")){
    stop("reads must be a ShortReadQ object or a file path to a fastq file.")
  }
  if(is.character(barcodes)){
    barcodes <- utils::read.delim(barcodes)
  }
  if(!is(barcodes, "data.frame")){
    stop("barcodes must be a data.frame or a file path to a file readable by read.delim.")
  }
  if(!all(c(fwd, rev, name) %in% colnames(barcodes))){
    stop("fwd, rev, and name must all be column names in the barcodes dataframe.")
  }
  mode = mode[1]
  if(mode == "dada"){
    out<-.demultiplex_for_dada(reads, barcodes, fwd, rev, name, writeOut, stat, cores)
  } else if (mode=="general"){
    out<-.demultiplex_general(reads, barcodes, fwd, rev, name, writeOut, stat, cores)
  } else {
    stop("mode must be one of 'dada' or 'general'")
  }
  if(!is.logical(writeOut)){writeOut = TRUE}
  if(stat & writeOut){
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

.demultiplex_for_dada <- function(reads, barcodes, fwd="FWD", rev="REV",name = c("Name", "Sample_name"), writeOut = TRUE, stat=FALSE, cores = 1){
  # format barcodes to uppercase
  barcodes[[fwd]] <- toupper(barcodes[[fwd]])
  barcodes[[rev]] <- toupper(barcodes[[rev]])
  # per barcode write out a fq and optionally a summary table
  test<-parallel::mclapply(1:nrow(barcodes), function(i){
    # get a name for the sub sample
    name <- paste(barcodes[i, c(name)], collapse=".")
    # get barcodes
    fwd_bar <- barcodes[i,fwd]
    rev_bar <- barcodes[i,rev] 
    # filter for barcodes at beginning of sequence
    RowReads    <- reads[grepl(paste0("^",fwd_bar), ShortRead::sread(reads))]
    ColumnReads <- reads[grepl(paste0("^",rev_bar), ShortRead::sread(reads))]
    # make reverse complements of fwd and reverse barcodes
    fwd_bar_reverseComplement <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(fwd_bar)))
    rev_bar_reverseComplement <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(rev_bar)))
    # search for reverse complements of opposite barcode at end of sequences
    RowReads <- RowReads[grepl(paste0(rev_bar_reverseComplement, "$"), ShortRead::sread(RowReads))]
    ColumnReads <- ColumnReads[grepl(paste0(fwd_bar_reverseComplement, "$"), ShortRead::sread(ColumnReads))]
    # append row and column reads
    Clone <- append(RowReads, ColumnReads)
    # write fastq
    if(!is.logical(writeOut)){
      write = writeOut
      if(substr(write, nchar(write), nchar(write))!="/"){
        write <- paste0(write, "/")
      }
      writeOut = TRUE} else{
        write = NULL
    }
    if(writeOut){
      ShortRead::writeFastq(object = Clone, file =paste0(write, name, ".fq"),mode = "w", compress = T)
    }
    if(stat & writeOut){
      stats <- data.frame(name = name,
                         totalReads = length(reads),
                         groupedReads = ShortRead::countFastq(paste0(name, ".fq"))[1,1])
      stats$ratio <- signif(stats$groupedReads/stats$totalReads, digits = 4)
      write.csv(stats, file=paste0(name, ".csv"), row.names = FALSE)
      return(stats)
    }
  }, mc.cores = cores)
  if(!is.logical(writeOut)){writeOut = TRUE}
  if(stat & writeOut){
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

.demultiplex_general <- function(reads, barcodes, fwd="FWD", rev="REV",name = c("Name", "Sample_name"), writeOut = TRUE, stat=FALSE, cores=1){
  # format barcodes to uppercase
  barcodes[[fwd]] <- toupper(barcodes[[fwd]])
  barcodes[[rev]] <- toupper(barcodes[[rev]])
  # per barcode write out a fq and optionally a summary table
  test<-lapply(1:nrow(barcodes), function(i){
    # get a name for the sub sample
    name <- paste(barcodes[i, c(name)], collapse=".")
    # get barcodes
    fwd_bar <- barcodes[i,fwd]
    rev_bar <- barcodes[i,rev] 
    # filter for barcodes at beginning of sequence
    RowReads    <- reads[grepl(paste0("^",fwd_bar), ShortRead::sread(reads))]
    ColumnReads <- reads[grepl(paste0("^",rev_bar), ShortRead::sread(reads))]
    # make reverse complements of fwd and reverse barcodes
    fwd_bar_reverseComplement <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(fwd_bar)))
    rev_bar_reverseComplement <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(rev_bar)))
    # search for reverse complements of opposite barcode at end of sequences
    RowReads <- RowReads[grepl(paste0(rev_bar_reverseComplement, "$"), ShortRead::sread(RowReads))]
    ColumnReads <- ColumnReads[grepl(paste0(fwd_bar_reverseComplement, "$"), ShortRead::sread(ColumnReads))]
    # append row and column reads
    Clone <- append(RowReads, ColumnReads)
    # write fastq
    if(!is.logical(writeOut)){
      write = writeOut
      if(substr(write, nchar(write), nchar(write))!="/"){
        write <- paste0(write, "/")
      }
      writeOut = TRUE} else{
        write = NULL
      }
    if(writeOut){
      if(i==1){md = "w"}else{md="a"}
      ShortRead::writeFastq(object = Clone, file =paste0(write, name, ".fq"), mode = md, compress = T)
    }
  })
  if(stat & writeOut){
    stats <- data.frame(name = name,
                       totalReads = length(reads),
                       groupedReads = ShortRead::countFastq(paste0(name, ".fq"))[1,1])
    stats$ratio <- signif(stats$groupedReads/stats$totalReads, digits = 4)
    write.csv(stats, file=paste0(name, ".csv"), row.names = FALSE)
    return(stats)
  }
}
