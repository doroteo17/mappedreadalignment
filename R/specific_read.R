#' Access the specific entry in the BAM file
#'
#' This function takes a BAM file and the specific read ID as input and
#' return a list with the id, chromosome, sequence, CIGAR and
#' the position of the specific read. If the input file is not a BAM file,
#' or the specific read id is not found, an error is raised.
#'
#'
#' @param file_bam Path to the input BAM file.
#' @param specific_read_id The ID of the read to extract.
#' @return A list containing the specific info of the read.
#' @importFrom Rsamtools scanBam ScanBamParam scanBamHeader
#' @importFrom Biostrings DNAString
#' @examples
#' specific_read(system.file("extdata", "example.bam",
#' package = "mappedreadalignment"), "read4")
#' # Returns: $qname
#' # [1] "read4"
#' #
#' # $rname
#' # [1] "chr6"
#' #
#' # $seq
#' # 11-letter DNAString object
#' # seq: AAAGATCGACC
#' #
#' # $cigar
#' # [1] "2S5M2D4M"
#' #
#' # $pos
#' # [1] 32401000
#' @export
specific_read <- function(file_bam, specific_read_id) {
    if (!file.exists(file_bam)) {
        stop("File does not exist!")}

    tryCatch({
        Rsamtools::scanBamHeader(file_bam)},
        error = function(e) {
        stop("This is not a valid BAM file!")})

    bam <- Rsamtools::scanBam(file_bam, param = Rsamtools::ScanBamParam(
        what = c("qname", "seq", "cigar","rname","pos")))
    reads <- bam[[1]]
    idx <- which(reads$qname == specific_read_id)

    if (length(idx) == 0) {
        stop("Read ID not found in BAM/SAM file.")}

    list(qname = as.character(reads$qname[idx]),
        rname = as.character(reads$rname[idx]),
        seq = Biostrings::DNAString(as.character(reads$seq[idx])),
        cigar = as.character(reads$cigar[idx]),
        pos = reads$pos[idx])}
