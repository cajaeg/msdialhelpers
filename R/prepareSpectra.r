#' Prepare spectra for annotation process
#'
#' This function converts text-encoded spectra to data.frame objects, applies an intensity threshold and adds data.frame columns needed for annotation.
#' @param x MS-DIAL results table as obtained by \code{\link{loadAlignmentResults}()} or \code{\link{loadConsoleResults}()}
#' @param in_column name of column to be converted, e.g. "ms1" (default), "ms2", or "spectrum"
#' @param out_column name of column containing the prepared spectra
#' @param relthr relative intensity threshold to be applied to each spectrum
#' @return data.frame
#' @export
#'
#' @examples
#' fp <- system.file("extdata/MSDconsole_single_result.msdial", package = "msdialhelpers")
#' dat <- loadConsoleResults(fp)
#' colnames(dat) # spectrum column is 'spectrum'
#' dat <- prepareSpectra(dat, in_column = "spectrum")
#' head(dat$s[[1]])
#' plot(dat$s[[1]][,1:2], type = "h")
prepareSpectra <- function(x, in_column = "ms1", out_column = "s", relthr = 0.001) {
  stopifnot(inherits(x, c("matrix", "data.frame", "tbl", "tbl_df")))
  stopifnot(in_column %in% colnames(x))
  s_chr <- x[[in_column]]
  s_mat <- lapply(s_chr, str2spec)
  s_mat <- lapply(s_mat, function(x) {
    x[,"i"] <- ifelse(is.finite(x[,"i"]), x[,"i"], 0)
    return(x)
  })
  s_mat <- lapply(s_mat, function(x) {
    x[ x[,"i"] >= max(x[,"i"]) * relthr, ]
    return(x)
  })
  s_df <- lapply(s_mat, as.data.frame, stringsAsFactors = FALSE)
  s_df <- lapply(s_df, function(x) {
    x$isogr <- NA_integer_
    x$iso <- NA_integer_
    x$charge <- NA_integer_
    x$formula <- NA_character_
    x$exactmass <- NA_complex_
    x$error <- NA_character_
    x$label <- NA_character_
    x$peakcol <- 1
    x$labelcol <- 1
    return(x)
  })
  x[[out_column]] <- s_df
  return(x)
}
