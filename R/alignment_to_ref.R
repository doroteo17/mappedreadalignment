#' Align read to reference genome based on CIGAR string
#'
#' This function takes as input a specific read extracted from a BAM file and
#' a reference genome, if no reference genome is given in input,
#' a UCSC hg38 reference genome is used, and return the alignment of the read
#' to the reference genome, using the CIGAR of the specific read.
#' All the mapped reads, including those that have been mapped to the
#' reverse strand,in the BAM file are represented on the forward
#' genomic strand and also CIGAR string is reversed and thus recorded
#' consistently with the sequence bases as represented. Taking this into
#' consideration the alignment of a query read, that have been mapped to the
#' reverse strand, to the genome is also done on the forward genomic strand.
#'
#' @param read A specific read extracted from a BAM file.
#' @param genome An optional reference genome given in input by the user.
#' @return The function return a list object reporting the read ID, CIGAR string
#' and the alignment of the query read to the reference genome.
#' @importFrom Biostrings readDNAStringSet subseq
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @examples
#' example_read <- list(rname = "chr6",
#' pos = 32402000,
#' strand = "+",
#' cigar = "2S5M2D4M",
#' seq = "AAAGATCGACC",
#' qname = "read4")
#' alignment_to_ref(example_read)
#' # Returns: $Read_ID
#' # [1] "read4"
#' #
#' # $CIGAR
#' # [1] "2S5M2D4M"
#' #
#' # $Alignment
#' # $Alignment$Reference_genome
#' # [1] "  AGATCGAGACC"
#' #
#' # $Alignment$Query_read
#' # [1] "AAAGATC--GACC"
#' @export
alignment_to_ref <- function(read, genome = NULL) {
  if (is.null(genome)) {
    if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE))
      stop("The package 'BSgenome.Hsapiens.UCSC.hg38' is required but is not installed.")
    genome <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens }
  else {
    genome <- Biostrings::readDNAStringSet(genome, format = "fasta")}
  chr <- as.character(read$rname)
  cigar <- cigar_transform(read$cigar)
  pos <- read$pos
  seq <- as.character(read$seq)
  read_pos <- 1
  ref_read <- c()
  query_read <- c()
  genome <- genome[[chr]]
  for (operation in cigar) {
    if (operation == "M" || operation == "=" || operation == "X") {
      if (read_pos <= nchar(seq)) {
        ref_read <- c(ref_read, as.character(subseq(genome, start = pos,width = 1)))
        query_read <- c(query_read, substr(seq, read_pos,read_pos))
        read_pos <- read_pos + 1
        pos <- pos + 1 }
    } else if (operation == "I") {
      if (read_pos <= nchar(seq)) {
        ref_read <- c(ref_read, "-")
        query_read <- c(query_read, substr(seq, read_pos,read_pos))
        read_pos <- read_pos + 1 }
    } else if (operation == "D") {
      ref_read <- c(ref_read, as.character(subseq(genome, start = pos,width = 1)))
      query_read <- c(query_read, "-")
      pos <- pos + 1
    } else if (operation == "S") {
      query_read <- c(query_read, substr(seq, read_pos,read_pos))
      ref_read <- c(" ", ref_read)
      read_pos <- read_pos + 1
    } else if (operation == "N") {
      ref_read <- c(ref_read, as.character(subseq(genome, start = pos,width = 1)))
      query_read <- c(query_read, "_")
      pos <- pos + 1
    } else if (operation == "P") {
      ref_read <- c(ref_read, "*")
      query_read <- c(query_read, "*")
    }
  }
  ref_read <- paste(ref_read, collapse = "")
  query_read <- paste(query_read, collapse = "")
  result <- list(
    Read_ID = read$qname,
    CIGAR = read$cigar,
    Alignment = list(Reference_genome = ref_read,
                     Query_read = query_read))
  return(result)
}
