#' Numeric equivalent to 'match'
#' 
#' @param x numeric vector
#' @param table numeric vector against which to match x
#' @param delta_x allowed absolute difference
#' @return array index
#' @export
#'
#' @examples
#' x <- c(3:5)+0.01
#' y <- 1:10
#' nummatch(x, y, delta_x = 0.1)
nummatch <- function(x, table, delta_x = 0.005) {
  do_match <- function(x, table, delta_x) {
    which(abs(table - x) < delta_x)
  }
  pidx <- mapply(do_match, x, MoreArgs = list(table=table, delta_x=delta_x))
  return(cbind(
      x = rep(seq_along(pidx), sapply(pidx, length)),
      table = unlist(pidx)
  ))
}
