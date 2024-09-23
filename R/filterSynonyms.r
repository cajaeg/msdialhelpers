##' Filter chemical synonyms as obtained from querying Pubchem
##'
##' @param x character vector
##' @return character vector
##' @export
##'
filterSynonyms <- function(x) {
  if (length(x) > 1) {
    ## remove CAS strings
    flt <- !grepl("[0-9]+-[0-9]+-[0-9]+", x)
    x <- x[flt]
    if (length(x) > 1) {
      charwise <- strsplit(x, "")
      ## remove UPPER CASE-only names
      flt <- !sapply(charwise, function(x)
        all(x %in% LETTERS[1:26]))
      x <- x[flt]
      if (length(x) > 1) {
        charwise <- strsplit(x, "")
        ## check ratio of upper-to-lower case
        flt <- sapply(charwise, function(x)
          sum(x %in% LETTERS[1:26]) / sum(x %in% letters[1:26]) < 0.3)
        x <- x[flt]
        if (length(x) > 1) {
          ## remove entries containing purities e.g. '99.0%'
          flt <- !grepl("[0-9\\.]+\\%", x)
          x <- x[flt]
        }
        if (length(x) > 1) {
          charwise <- strsplit(x, "")
          ## check ratio of numbers to letters
          flt <- sapply(charwise, function(x)
            sum(x %in% as.character(0:9)) / sum(!x %in% as.character(0:9)) < 0.3)
          x <- x[flt]
        }
      }
    }
  }
  return(x)
}
