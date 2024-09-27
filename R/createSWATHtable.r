#' Create an "MS-DIAL experiment information" file required for SWATH/MS^E data
#'
#' @param mzrange Method m/z range (scan range)
#' @param winwidth SWATH window width (Da)
#' @param ms2range m/z range to create SWATH windows for, often identical to
#'   \code{mzrange}
#' @param spectra_rate Only used to calculate accumulation time (\code{attr(out,
#'   "acc_time")})
#' @param ms1_acq_factor Increase accumulation time for MS1 scan by this factor
#' @param outfile If non-NULL, divert output to text file that can be used as
#'   MS-DIAL experiment information file
#'
#' @return data.frame, with additional attributes "cycle_time", "acc_time"
#'   (accumulation time) and "bruker_mrm_table". The latter can serve as
#'   template in Bruker DataAcquisition.
#' @export
#'
#' @examples
#' createSWATHtable()
#' \dontrun{
#' createSWATHtable(outfile = "MSDIAL_exp_info.txt")
#' }
#' out <- createSWATHtable(c(50, 600))
#' out
#' attr(out, "cycle_time")
#' attr(out, "acc_time")
#' attr(out, "bruker_mrm_table")
#' createSWATHtable(c(50, 500),
#'   winwidth = 25,
#'   ms1_acq_factor = 1.67,
#'   spectra_rate = 30) # reproduce HILIC SWATH scheme from Tsugawa et al. 2015 (cycle time 640 ms)
#' createSWATHtable(c(100, 1250),
#'   winwidth = 21,
#'   ms1_acq_factor = 10,
#'   spectra_rate = 90) # reproduce lipid SWATH scheme (cycle time 730 ms)
#' createSWATHtable(c(100, 1250),
#'   winwidth = 25,
#'   ms2range = c(350, 1250),
#'   ms1_acq_factor = 5,
#'   spectra_rate = 50) # compromise for Impact II (max scan rate=50)
createSWATHtable <- function(mzrange = c(50, 500),
                             winwidth = 25,
                             ms2range = NULL,
                             spectra_rate = 50,
                             ms1_acq_factor = 5,
                             outfile = NULL) {
  op <- options(stringsAsFactors = F)
  on.exit(options(op))
  if (is.null(ms2range))
    ms2range <- mzrange
  win_lower <- seq(from = ms2range[1],
                   to = ms2range[2] - winwidth,
                   by = winwidth)
  win_upper <- seq(from = ms2range[1] + winwidth,
                   to = ms2range[2],
                   by = winwidth)
  win_center <- seq(
    from = ms2range[1] + winwidth / 2,
    to = ms2range[2] - winwidth / 2,
    by = winwidth
  )
  nwin <- length(win_center)
  acq_factor <- c(ms1_acq_factor, rep(1, nwin))
  out_colnames <- c("Experiment", "MS Type", "Min m/z", "Max m/z")
  out <- data.frame(matrix(
    nrow = nwin + 1,
    ncol = length(out_colnames),
    dimnames = list(NULL, out_colnames)
  ),
  check.names = F)
  out[, 1] <- (1:nrow(out)) - 1
  out[, 2] <- c("SCAN", rep("SWATH", nwin))
  out[, 3] <- c(mzrange[1], win_lower)
  out[, 4] <- c(mzrange[2], win_upper)
  attr(out, "cycle_time") <- sum(acq_factor) / spectra_rate
  attr(out, "acc_time") <- acq_factor / spectra_rate
  attr(out, "bruker_mrm_table") <- data.frame(
    Mass = c(mzrange[1], win_center),
    Width = c(0, rep(winwidth, nwin)),
    Collision = c(8, rep(40, nwin)),
    xAcq = acq_factor
  )
  if (is.null(outfile)) {
    cat("Cycle time:", attr(out, "cycle_time"), "\n")
    cat("Accumulation times:", unique(attr(out, "acc_time")), "\n")
    return(out)
  } else {
    utils::write.table(
      out,
      file = outfile,
      sep = "\t",
      quote = F,
      col.names = T,
      row.names = F
    )
  }
}
