#' Import alkane definition files
#'
#' 
#' @param txtfile text file with two tab-separated columns (Num|RT(min))
#' @param plot show diagnostic plot?
#' @return data.frame
#' @export
#'
#' @examples
#' fp <- system.file("extdata/Alkane.txt", package = "msdialhelpers")
#' loadAlkanes(fp, plot = TRUE)
loadAlkanes <- function(txtfile, plot = FALSE) {
  stopifnot(file.exists(txtfile))
  out <- utils::read.csv(txtfile, sep = "\t", as.is = TRUE)
  stopifnot(ncol(out) == 2)
  out$RI <- out[, 1] * 100
  colnames(out)[1:2] <- c("Num", "RT")
  if (plot)
    graphics::plot(out[, c("RI", "RT")])
  return(out)
}
