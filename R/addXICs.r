#' Add extracted ion chromatograms to alignment results
#'
#' @param x data.frame as created by \code{\link{loadAlignmentResults}()}
#' @param in_column name of column containing m/z's to extract: either a list of
#'   spectra or single m/z values
#' @param rt_column name of column containing retention times
#' @param out_column name of column where results are
#' @param xraw list of 'xcmsRaw' objects
#' @param mz_dev width of m/z window to extract (passed to
#'   \code{\link{getChrom}()})
#' @param rt_dev width of RT window to extract (passed to
#'   \code{\link{getChrom}()})
#' @param zeroVal value to replace NA's by (passed to \code{\link{getChrom}()})
#' @param smooth (passed to \code{\link{getChrom}()})
#' @param max_mz restrict BPC extraction to the 'max_mz' highest peaks (if
#'   'in_column' contains spectra)
#' @param EIC if 'TRUE' add (sum-based) EICs instead of (maximum-based) BPCs
#' @param .pivot_longer transform resulting BPC matrices to long format. Useful
#'   for plotting with 'ggplot'.
#' @return data.frame
#' @export
#'
#' @examples
#' # see ?loadAlignmentResults
#'
addXICs <- function(x,
                    in_column = "s",
                    rt_column = c("average_rt_min", "rt_min"),
                    out_column = "xic",
                    xraw,
                    mz_dev = 0.01,
                    rt_dev = 10,
                    zeroVal = NA,
                    smooth = 0,
                    max_mz = 5,
                    EIC = FALSE,
                    .progress = FALSE,
                    .pivot_longer = FALSE) {
  stopifnot(in_column %in% colnames(x))
  stopifnot(any(rt_column %in% colnames(x)))
  mz_col <- x[[in_column]]
  mz <- if (is.list(mz_col))
    # assume in_column holds a list of spectra
    lapply(mz_col, function(spec)
      spec[order(-spec[, 2]), ][1:min(max_mz, nrow(spec)), ][, 1])
  else
    mz <- as.numeric(unlist(mz_col))
  rt_column <- rt_column[which(rt_column %in% colnames(x))[1]]
  rt <- x[[rt_column]]
  out <- vector("list", nrow(x))
  if(.progress) {
    pb <- utils::txtProgressBar(min = 0,
                                max = length(out),
                                style = 3)
    on.exit(close(pb))
  }
  for (i in 1:length(out)) {
    if(.progress) {
      utils::setTxtProgressBar(pb, i)
    }
    tmp <- vector("list", length(xraw))
    flt <- sapply(xraw, inherits, "xcmsRaw")
    for (j in seq_along(xraw)[flt]) {
      tmp[[j]] <- getChrom(
        xraw[[j]],
        mz = mz[[i]],
        mz_dev = mz_dev,
        rt = rt[[i]] * 60,
        rt_dev = rt_dev,
        zeroVal = zeroVal,
        smooth = smooth,
        EIC = EIC
      )
    }
    names(tmp)[flt] <- sub("\\.mzML$", "", basename(sapply(xraw[flt], methods::slot, "filepath")))
    out[[i]] <- if (.pivot_longer) {
      plyr::ldply(tmp, function(x)
        data.frame(rt = attr(x, "rt"), x, check.names = F), .id = "fromFile") |>
        tidyr::pivot_longer(
          cols = !tidyselect::matches("fromFile|rt"),
          names_to = "mz",
          values_to = "i"
        ) |>
        dplyr::mutate(mz = factor(round(as.numeric(mz), 4), levels = round(attr(tmp[[1]], "mz"), 4)))
    } else
      tmp
  }
  x[[out_column]] <- out
  return(x)
}
