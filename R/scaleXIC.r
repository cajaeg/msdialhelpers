#' Scale chromatograms
#'
#' @param xic chromatogram(s) as obtained by \code{\link{getXIC}()}
#' @param scope one of "global" (default), "by_file", "by_trace" or "none"
#'
#' @return matrix or list equivalent to original 'xic'
#' @export
#'
#' @examples
#' # see getXIC()
scaleXIC <- function(xic,
                       scope = c("global", "by_file", "by_trace", "none")[1]) {
  if (inherits(xic, "list")) {
    scope <- match.arg(scope, choices = c("global", "by_file", "by_trace", "none"))
    return(lapply(xic, scaleXIC, scope = switch(
      scope,
      global = max(1, unlist(xic), na.rm = TRUE),
      by_file = "by_file",
      by_trace = "by_trace",
      none = "none"
    )))
  } else {
    old_attr <- attributes(xic)
    out <- if (is.numeric(scope)) {
      xic / scope * 100
    } else {
      scope <- match.arg(scope, choices = c("global", "by_file", "by_trace", "none"))
      switch(
        scope,
        global = xic / max(1, xic, na.rm = TRUE) * 100,
        by_file = xic / max(1, xic, na.rm = TRUE) * 100,
        by_trace = apply(xic, 2, function(x)
          x / max(1, x, na.rm = TRUE) * 100),
        none = xic
      )
    }
    mostattributes(out) <- old_attr
    return(out)
  }
}

