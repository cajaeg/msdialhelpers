#' Export MSP library file
#'
#' @param x data.frame to be exported as MSP
#' @param file MSP file to be written
#' @param append append to existing file?
#' @param progress display progress bar?
#'
#' @return None
#' @export
#'
#' @examples
#' msp <- data.frame(
#'   Name = "Test_spectrum",
#'   `Num peaks` = 3,
#'   PEAKS = "100:999 101:122 102:33",
#'   check.names = FALSE
#' )
#' msp
#' \dontrun{
#' writeMSP(msp, file = "test.msp")
#' }
writeMSP <- function (x,
                      file,
                      append = FALSE,
                      progress = FALSE)
{
  format_section <- function(x, normal_columns, peak_columns, nr) {
    utils::write.table(
      paste(colnames(x)[c(normal_columns, peak_columns[1])], x[, c(normal_columns, peak_columns[1])], sep = ": "),
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE,
      file = ""
    )
    utils::write.table(
      str2spec(x[peak_columns[2]]),
      row.names = FALSE,
      col.names = FALSE,
      file = ""
    )
    if (nr > 1)
      cat("\n", file = "")
  }
  peak_columns <- match(c("num peaks", "peaks"), tolower(colnames(x)))
  normal_columns <- (1:ncol(x))[-peak_columns]
  nr <- nrow(x)
  out <- vector("list", nr)
  if (progress)
    pb <- utils::txtProgressBar(0, nr, style = 3)
  for (i in 1:nr) {
    out[[i]] <- utils::capture.output(format_section(x[i, ], normal_columns, peak_columns, nr))
    if (progress)
      utils::setTxtProgressBar(pb, i)
  }
  if (progress)
    close(pb)
  cat(unlist(out),
      sep = "\n",
      file = file,
      append = append)
}
