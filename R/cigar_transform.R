#' Transform CIGAR into a character vector
#'
#' This function takes a CIGAR string as input and transforms
#' it into a character vector
#' by expanding each operation symbol according to the number of times it occurs.
#'
#' @param cigar A CIGAR string extracted from a BAM file.
#' @return  A character vector where each element represents an individual
#' operation from the CIGAR string.
#' @examples
#' cigar_transform("2M1I3M1D")
#' # Returns: [1] "M" "M" "I" "M" "M" "M" "D"
#' @export
cigar_transform <- function(cigar) {

  numbers <- as.integer(unlist(regmatches(cigar, gregexpr("\\d+", cigar))))
  code <- unlist(regmatches(cigar, gregexpr("[A-Z]", cigar)))

  result <- ""

  for (i in seq_along(code)) {
    result <- paste0(result, strrep(code[i], numbers[i]))
  }

  return(unlist(strsplit(result, NULL)))
}
