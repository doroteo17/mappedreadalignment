#' Access the specific entry in the BAM file
#'
#' This function takes a BAM file and the specific read ID as input and
#' return a list with the id, chromosome, sequence, CIGAR, position and
#' strand of the specific read.
#'
#'
#' @param file_bam Path to the input BAM file.
#' @param specific_read_id The ID of the read to extract.
#' @return A list containing the specific info of the read.
#' @importFrom Rsamtools scanBam ScanBamParam
#' @importFrom Biostrings DNAString
#' @export
specific_read <- function(file_bam, specific_read_id) {
  bam <- Rsamtools::scanBam(file_bam, param =
                              Rsamtools::ScanBamParam(what = c("qname", "seq", "cigar",
                                                               "rname","pos","strand")))
  reads <- bam[[1]]
  idx <- which(reads$qname == specific_read_id)

  if (length(idx) == 0) {
    stop("Read ID not found in BAM/SAM file.")
  }

  list(qname = as.character(reads$qname[idx]),
       rname = as.character(reads$rname[idx]),
       seq = Biostrings::DNAString(as.character(reads$seq[idx])),
       cigar = as.character(reads$cigar[idx]),
       pos = reads$pos[idx],
       strand = reads$strand[idx])
}
