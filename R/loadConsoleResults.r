#' Import MS-DIAL console results
#' 
#' Import '.msdial' files created by MS-DIAL console app.
#' @param txtfile .msdial file for individual MS files. Use \code{\link{loadAlignmentResults}()} to import alignment results for multiple files.
#' @return data.frame
#' @export
#'
#' @examples
#' fp <- system.file("extdata/MSDconsole_single_result.msdial", package = "msdialhelpers")
#' dat <- loadConsoleResults(fp)
loadConsoleResults <- function(txtfile) {
  out <- utils::read.table(
    txtfile,
    header = TRUE,
    quote = "",
    sep = "\t",
    check.names = FALSE,
    as.is = TRUE,
    comment.char = ""
  )
  cn <- colnames(out)
  cn <- gsub("[^a-zA-Z0-9]+", "_", trimws(tolower(cn)))
  cn <- sub("_+$", "", cn)
  cn <- sub("^_+", "", cn)
  colnames(out) <- cn
  peakid <- if (nrow(out) > 0)
    1:nrow(out)
  else
    numeric()
  out <- data.frame(peakid = peakid, out, stringsAsFactors = FALSE)
  message(sprintf("peak data loaded for %d compounds", nrow(out)))
  return(out)
}
