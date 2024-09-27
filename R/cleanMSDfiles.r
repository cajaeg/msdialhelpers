#' Remove text files created by an MS-DIAL run
#'
#' @param msdir MS-DIAL working directory
#' @param pat file pattern passed to \code{list.files()}. The default
#'   pattern selects .pai2, .dcl, .msdial and .aef files.
#' @param ask ask before deleting
#' @return (invisibly) character vector of file names
#' @export
#'
cleanMSDfiles <- function(msdir, pat = "\\.pai2$|\\.dcl$|\\.msdial$|\\.aef$", ask = TRUE) {
  fls <- list.files(msdir, pat, full.names = TRUE)
  if (length(fls) > 0) {
    print(fls)
    if (ask) {
      if (utils::askYesNo("Delete these files?")) {
        file.remove(fls)
      }
    } else {
      file.remove(fls)
    }
  } else {
    message("no files found")
  }
  invisible(fls)
}
