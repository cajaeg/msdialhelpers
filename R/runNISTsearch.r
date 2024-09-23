##' Annotate MS-DIAL alignment results using MSPepSearch
##'
##' @title Annotate MS-DIAL spectra using NIST search
##' @param x MS-DIAL data.frame as returned by \code{\link{loadConsoleResults}()} or \code{\link{loadAlignmentResults}()}
##' @param in_column name of column containing prepared spectra (see \code{\link{prepareSpectra}()})
##' @param ri_column name of column containing RI values or 'NULL' to exclude RI
##' @param out_column name of column added to data.frame
##' @param ... further arguments passed to \code{\link{searchNIST}()}
##' @return data.frame
##' @export
runNISTsearch <- function(x,
                          in_column = "s",
                          ri_column = "average_ri",
                          out_column = "NISTres",
                          ...) {
    stopifnot(inherits(x, c("matrix", "data.frame", "tbl", "tbl_df")))
    stopifnot(in_column %in% colnames(x))
    spectra <- x[[ in_column ]]
    if(!is.null(ri_column)) {
        ri_obs <- x[[ ri_column ]]
    } else {
        ri_obs <- NULL
    }
    tmp <- searchNIST(spec = spectra, ri_obs = ri_obs, ...)
    tmp$query <- factor(tmp$query, levels = 1:length(spectra))
    x[[ out_column ]] <- split(tmp, tmp$query)
    return(x)
}
