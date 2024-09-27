#' Extract ion chromatograms (BPC, EIC, TIC)
#'
#' Wrapper for \code{\link[HiResTEC]{getMultipleBPC}()} accepting lists of
#' 'xcmsRaw' objects.
#' @param xraw 'xcmsRaw' object or list of 'xcmsRaw' objects
#' @param mz mass vector or NULL to extract the full mass range
#' @param mz_dev allowed m/z deviations, see
#'   \code{\link[HiResTEC]{getMultipleBPC}()} for details
#' @param rt retention time point in seconds or NULL to extract the full RT
#'   range
#' @param rt_dev RT window
#' @param zeroVal value to use for intensities <= 0 (typically NA or 0)
#' @param smooth window size for moving average smoother, 0 = no smoothing
#' @param EIC if 'TRUE' return (sum-based) EIC instead of (maximum-based) BPC
#'
#' @return A matrix with scan wise (rows) intensities for all requested m/z's
#'   (columns), or a list of such matrices.
#' @export
#'
#' @examples
#' mzml_files <- list.files(system.file("extdata", package = "msdialhelpers"),
#' "GCQTOF_PAH.*\\.mzML", full.names = TRUE)
#' xraw <- loadRaws(mzml_files)
#'
#' tic <- getXIC(xraw, EIC = TRUE) # total ion chromatogram (TIC)
#' plotXIC(tic)
#'
#' bpc <- getXIC(xraw, EIC = FALSE) # full-range base peak chromatogram (BPC)
#' plotXIC(bpc)
#'
#' xic <- getXIC(xraw, mz = c(202.078, 101.042), mz_dev = 0.005, rt = 20*60, rt_dev = 1*60)
#' plotXIC(xic)
#'
#' par(mfrow = c(1,2))
#' plotXIC(xic, set.mfrow = FALSE) # override default multi-figure layout
#' 
#' xic_scaled <- scaleXIC(xic, scope = "by_file")
#' plotXIC(xic_scaled)
getXIC <- function(xraw,
                   mz = NULL,
                   mz_dev = 0.01,
                   rt = NULL,
                   rt_dev = 2,
                   zeroVal = NA,
                   smooth = 0,
                   EIC = FALSE) {
  checkFinite <- function(x)
    if (!is.null(x) && any(!is.finite(x)))
      NULL
  else
    x
  mz <- checkFinite(mz)
  rt <- checkFinite(rt)
  if (inherits(xraw, "list")) {
    .n <- length(xraw)
    flt <- sapply(xraw, inherits, "xcmsRaw")
    out <- vector("list", .n)
    for (i in which(flt)) {
      out[[i]] <- getXIC(
        xraw[[i]],
        mz = mz,
        mz_dev = mz_dev,
        rt = rt,
        rt_dev = rt_dev,
        zeroVal = zeroVal,
        smooth = smooth,
        EIC = EIC
      )
    }
    return(out)
  } else {
    out <- HiResTEC::getMultipleBPC(
      xraw,
      mz = mz,
      mz_dev = mz_dev,
      rt = rt,
      rt_dev = rt_dev,
      zeroVal = zeroVal,
      smooth = smooth,
      returnEIC = EIC
    )
    if (is.null(out)) {
      out <- matrix(nrow = 0, ncol = length(mz))
      attr(out, "mz") <- mz
      attr(out, "mz_dev") <- mz_dev
    }
    attr(out, "from_file") <- basename(xraw@filepath[])
    return(out)
  }
}
