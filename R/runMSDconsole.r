#' Trigger MS-DIAL console process
#'
#' @param msdir directory holding MS files supported by MS-DIAL (.abf, .mzML, .mzXML, .cdf etc.)
#' @param analysis_type one of "gcms", "lcmsdda", "lcmsdia", "lcimmsdda", "lcimmsdia"
#' @param method list of parameters as returned by \code{\link{createMSDCmethod}()} or name of MS-DIAL method file (.txt)
#' @param MSDCexe full path of MS-DIAL console executable (MsdialConsoleApp.exe)
#' @param readResults read .msdial using \code{\link{loadConsoleResults}()} or \code{\link{loadAlignmentResults}()}
#' @param skipRun don't run MS-DIAL console, just read existing .msdial files
#' @param verbose show output of MS-DIAL console
#' @return list of data.frame's containing MS-DIAL results or vector of .msdial files, depending on 'readResults'
#' @export
#'
#' @examples
#' \dontrun{
#' mzml_files <- list.files(system.file("extdata", package = "msdialhelpers"), 
#'   "GC-QTOF.*\\.mzML", 
#'   full.names = TRUE)
#' msdir <- "~/MSDCtest"
#' dir.create(msdir)
#' file.copy(mzml_files, msdir)
#' (resfiles <- runMSDconsole(msdir))
#' peaks <- loadConsoleResults(resfiles[2])
#' aligned <- loadAlignmentResults(resfiles[1])
#' }
#'
#' ## use a custom method file
#' \dontrun{
#' (resfiles <- runMSDconsole(msdir, method = createMSDCmethod(minimum_peak_height = 5000)))
#' peaks <- loadConsoleResults(resfiles[2])
#' }
runMSDconsole <- function(msdir,
                          analysis_type = c("gcms", "lcmsdda", "lcmsdia", "lcimmsdda", "lcimmsdia")[1],
                          method = createMSDCmethod(),
                          MSDCexe = getOption("msdialhelpers_msdialconsole_path"),
                          readResults = TRUE,
                          skipRun = FALSE,
                          verbose = TRUE
                          ) {
  MSDCexe <- normalizePath(getOption("msdialhelpers_msdialconsole_path"))
  if(is.null(MSDCexe) || !file.exists(MSDCexe))
    stop("Set option 'msdialhelpers_msdialconsole_path' to the full path of MSD console executable.")
  msdir <- normalizePath(msdir)
  stopifnot(dir.exists(msdir))
  if(!skipRun) {
      if(is.list(method)) {
          if(!is.null(attr(method, "ri_index_file_paths"))) {
            mzml_files <- normalizePath(list.files(msdir, "\\.mzml", full.names = TRUE, ignore.case = TRUE))
            ri_index_file_path_outfile <- normalizePath(file.path(msdir, basename(attr(method, "ri_index_file_paths"))))
            cat(sprintf("%s\t%s", mzml_files, attr(method, "ri_index_file")),
                file = ri_index_file_path_outfile)
            method[["RI index file pathes"]] <- ri_index_file_path_outfile
          }
          method_file <- file.path(msdir, "MSDCparams.txt")
          cat(sprintf("%s: %s", names(method), unlist(method)), sep = "\n", file = method_file)

      } else if(is.character(method) && file.exists(method)) {
          method_file <- method
      }
      method_file <- normalizePath(method_file)
      stopifnot(file.exists(method_file))
      system2(MSDCexe, args = c(analysis_type,
                                "-i", utils::shortPathName(msdir),
                                "-o", utils::shortPathName(msdir),
                                "-m", utils::shortPathName(method_file)),
              stdout = if(verbose) "" else FALSE,
              stderr = if(verbose) "" else FALSE)
  }
  msdial_files <- list.files(msdir, "\\.msdial$", full.names = TRUE)
  if(readResults) {
      tryCatch({
          out <- vector("list", 2)
          names(out) <- c("aligned", "peaks")
          alignment_file <- grepl("AlignResult", msdial_files)
          if(any(alignment_file)) {
              fls <- msdial_files[which(alignment_file)]
              alignment_results <- lapply(fls, loadAlignmentResults, raw_folder = msdir)
              names(alignment_results) <- basename(sub("\\.msdial$", "", fls))
              out[[1]] <- alignment_results
          }
          single_files <- !alignment_file
          if(any(single_files)) {
              fls <- msdial_files[which(single_files)]
              single_results <-lapply(fls, loadConsoleResults)
              names(single_results) <- basename(sub("\\.msdial$", "", fls))
              out[[2]] <- single_results
          }
          return(out)
      }, error =
             function(e) {
                 message("Could not read results, return file names only")
                 return(msdial_files)
             })
  } else {
      return(msdial_files)
  }
}
