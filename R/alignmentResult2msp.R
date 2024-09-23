#' Create a new MSP library from MS-DIAL alignment results
#'
#' @param x alignment result data.frame as imported by loadAlignmentResult()
#' @param intabs intensity threshold applied to mass spectra (PEAKS field in the resulting msp)
#' @param makeNamesUnique apply a make.unique() to the spectra names (recommended)
#' @param digits decimal digits to round to (m/z, intensity)
#' @param ... (currently unused)
#'
#' @return data.frame
#' @export
#'
#' @examples
#' fp <- system.file("extdata/MSDIAL_Alignment_result_GC-MS.txt", package = "msdialhelpers")
#' aligned <- loadAlignmentResults(fp)
#' flt <- aligned$metabolite_name != "Unknown" # restrict to identified compounds
#' aligned <- aligned[flt, ]
#' msp <- alignmentResult2msp(aligned)
#' \dontrun{
#' writeMSP(msp, file = "library.msp")
#' }
alignmentResult2msp <- function(x,
                                intabs = 100,
                                makeNamesUnique = TRUE,
                                digits = c(4, 0),
                                ...) {
  spectra <- lapply(x$ms1, function(x) {
    x <- str2spec(x)
    x <- if (any(x[, 2] >= intabs))
      x[x[, 2] >= intabs, , drop = FALSE]
    else
      x[which.max(x[, 2]), , drop = FALSE]
    x[, 1] <- round(x[, 1], digits = digits[1])
    x[, 2] <- round(x[, 2] / max(x[, 2], na.rm = TRUE) * 999, digits = digits[2])
    x
  })
  spectra_txt <- lapply(spectra, spec2str)
  out <- data.frame(matrix(nrow = nrow(x), ncol = 0))
  tmp <- x$metabolite_name
  if (any(duplicated(tmp)) && makeNamesUnique)
    tmp <- make.unique(tmp, sep = "_")
  out$Name <- tmp
  out$Formula <- x$formula
  out$RETENTIONTIME <- x$average_rt_min
  out$RETENTIONINDEX <- x$average_ri
  out$SMILES <- x$smiles
  out$INCHIKEY <- x$inchikey
  out$Ontology <- x$ontology
  out$`Num peaks` <- unlist(lapply(spectra, nrow))
  out$PEAKS <- unlist(spectra_txt)
  out
}
