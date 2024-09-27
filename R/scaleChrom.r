#' Scale chromatograms
#'
#' @param chrom chromatogram(s) as obtained by \code{\link{getChrom}()}
#' @param scope one of "global" (default), "by_file", "by_trace" or "none"
#'
#' @return matrix or list equivalent to original 'chrom'
#' @export
#'
#' @examples
#' # see getChrom()
scaleChrom <- function(chrom,
                       scope = c("global", "by_file", "by_trace", "none")[1]) {
  if (inherits(chrom, "list")) {
    scope <- match.arg(scope, choices = c("global", "by_file", "by_trace", "none"))
    return(lapply(chrom, scaleChrom, scope = switch(
      scope,
      global = max(1, unlist(chrom), na.rm = TRUE),
      by_file = "by_file",
      by_trace = "by_trace",
      none = "none"
    )))
  } else {
    old_attr <- attributes(chrom)
    out <- if (is.numeric(scope)) {
      chrom / scope * 100
    } else {
      scope <- match.arg(scope, choices = c("global", "by_file", "by_trace", "none"))
      switch(
        scope,
        global = chrom / max(1, chrom, na.rm = TRUE) * 100,
        by_file = chrom / max(1, chrom, na.rm = TRUE) * 100,
        by_trace = apply(chrom, 2, function(x)
          x / max(1, x, na.rm = TRUE) * 100),
        none = chrom
      )
    }
    mostattributes(out) <- old_attr
    return(out)
  }
}
