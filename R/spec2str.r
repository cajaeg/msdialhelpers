#' Convert mass spectrum from character to matrix representation
#'
#' @param x character string encoding mass spectra, e.g. "100:999 101:121"
#' @param colnames column names in resulting data.frame (default \code{c("mz",
#'   "i")})
#'
#' @return data.frame
#' @export
#'
#' @examples
#' str2spec("100:999 101:121")
str2spec <- function (x, colnames = c("mz", "i"))
{
  tmp <- trimws(gsub("[^0-9\\.]*([0-9\\.]+)[^0-9\\.]+([0-9\\.]+)[^0-9\\.]*", "\\1 \\2 ", trimws(x)))
  tmp <- trimws(unlist(strsplit(tmp, " ")))
  tmp <- as.numeric(tmp)
  is_mz <- 1:length(tmp) %% 2 == 1
  out <- data.frame(mz = tmp[is_mz], i = tmp[!is_mz])
  base::colnames(out) <- colnames
  return(out)
}

#' Convert mass spectrum from matrix to character representation
#'
#' @param x two-column matrix or data.frame
#' @param sep separator to be used between 'mz' and 'intensity' (default ':')
#' @param collapse separator to be used between (m/z, i) blocks (default ' ')
#' @param round Number of decimal places to be used for mz and intensity,
#'   respectively (default \code{c(4,0)})
#'
#' @return character string
#' @export
#'
#' @examples
#' s <- data.frame(mz = c(100, 101), i = c(999, 121))
#' spec2str(s)
spec2str <- function(x, sep = ":", collapse = " ", round = c(4,0)) {
  if(!any(inherits(x, c("matrix", "data.frame"), which = TRUE) == 1)) {
    stop("no matrix or data.frame provided")
  }
  if(!is.null(round)) {
    x[,1] <- base::round(x[,1], round[1])
    x[,2] <- base::round(x[,2], round[2])
  }
  character_columns <- which(apply(x, 2, function(column) {
    methods::is(tryCatch(as.numeric(column), warning = function(w) w), "warning")
  }))
  if(length(character_columns) > 0) {
    for(i in character_columns) {
      ## quote character columns
      x[,i] <- sprintf("\"%s\"", x[,1])
    }
  }
  paste(apply(x, 1, function(x) {
    paste(trimws(x), sep="", collapse=":")
  }), sep="", collapse=" ")
}
