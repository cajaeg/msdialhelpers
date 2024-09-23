#' Calculate a simple linear score from RI deviation
#'
#' @param ri_dev deviation of observed RI from library RI
#' @param max_ri_dev deviation greater than this receive the min score; NA values also receive the min score
#' @param minscore minimum score
#' @param maxscore maximum score
#' @return vector
#' @export
#'
#' @examples
#' ri_obs <- 1225
#' ri_lib <- 1200
#' scoreRI(ri_obs - ri_lib)
scoreRI <- function(ri_dev, max_ri_dev = 100, minscore = 0.75, maxscore = 1) {
    ri_dev[!is.finite(ri_dev)] <- max_ri_dev
    ri_dev <- abs(ri_dev)
    ((pmin(ri_dev, max_ri_dev) - 0) / (max_ri_dev - 0)) * (minscore - maxscore) + maxscore
}
