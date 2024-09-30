#' Trigger MS-DIAL console process
#'
#' Add \code{options(msdialhelpers_msdialconsole_path = "~/LocalApps/MSDIAL
#' ver.4.9.221218 Windowsx64/MsdialConsoleApp.exe")} (adjusted accordingly to your system) to your .Rprofile file. Use
#' \code{Sys.getenv("R_USER")} to check where this file is located on your
#' system.
#' @param msdir directory holding MS files supported by MS-DIAL (.abf, .mzML,
#'   .mzXML, .cdf etc.)
#' @param analysis_type one of "gcms" (default), "lcmsdda", "lcmsdia",
#'   "lcimmsdda", "lcimmsdia"
#' @param method list of parameters as returned by
#'   \code{\link{createMSDCmethod}()} or, alternatively, name of existing
#'   MS-DIAL method file (.txt).
#' @param readResults import results files after processing
#' @param skipRun don't run MS-DIAL console, just read existing .msdial files
#' @param verbose show output of MS-DIAL console
#' @return List of data.frame's containing MS-DIAL results, or a vector of file
#'   names for 'readResults' == FALSE
#' @export
#'
#' @examples
#' \dontrun{
#' mzml_files <- list.files(system.file("extdata", package = "msdialhelpers"),
#'   "GCQTOF_PAH.*\\.mzML",
#'   full.names = TRUE)
#' workdir <- tempfile()
#' dir.create(workdir)
#' file.copy(mzml_files, workdir)
#' full_results <- runMSDconsole(workdir)
#' alignment_results <- full_results[[1]]
#' # check ?loadAlignmentsResults for examples of further processing these results
#' cleanMSDfiles(workdir, ask = FALSE)
#' }
#'
#' ## use a custom method file
#' \dontrun{
#' msdc_method <- createMSDCmethod(minimum_peak_height = 10000)
#' (result_files <- runMSDconsole(workdir, method = msdc_method, readResults = FALSE)
#' peaks <- loadConsoleResults(resfiles[[2]])
#' }
#'
#' ## clean up temporary folder
#' \dontrun{
#' unlink(workdir, recursive = TRUE, force = TRUE)
#' }
runMSDconsole <- function(msdir,
                          analysis_type = c("gcms", "lcmsdda", "lcmsdia", "lcimmsdda", "lcimmsdia")[1],
                          method = createMSDCmethod(),
                          readResults = TRUE,
                          skipRun = FALSE,
                          verbose = TRUE) {
  MSDCexe <- normalizePath(getOption("msdialhelpers_msdialconsole_path"))
  if (is.null(MSDCexe) || !file.exists(MSDCexe))
    stop(
      "Set option 'msdialhelpers_msdialconsole_path' to the full path of MSD console executable."
    )
  msdir <- normalizePath(msdir)
  stopifnot(dir.exists(msdir))
  if (!skipRun) {
    if (is.list(method)) {
      if (!is.null(attr(method, "ri_index_file_paths"))) {
        mzml_files <- normalizePath(list.files(
          msdir,
          "\\.mzml",
          full.names = TRUE,
          ignore.case = TRUE
        ))
        ri_index_file_path_outfile <- normalizePath(file.path(msdir, basename(
          attr(method, "ri_index_file_paths")
        )))
        cat(sprintf("%s\t%s", mzml_files, attr(method, "ri_index_file")), file = ri_index_file_path_outfile)
        method[["RI index file pathes"]] <- ri_index_file_path_outfile
      }
      method_file <- file.path(msdir, "MSDCparams.txt")
      cat(sprintf("%s: %s", names(method), unlist(method)),
          sep = "\n",
          file = method_file)
      
    } else if (is.character(method) && file.exists(method)) {
      method_file <- method
    }
    method_file <- normalizePath(method_file)
    stopifnot(file.exists(method_file))
    system2(
      MSDCexe,
      args = c(
        analysis_type,
        "-i",
        utils::shortPathName(msdir),
        "-o",
        utils::shortPathName(msdir),
        "-m",
        utils::shortPathName(method_file)
      ),
      stdout = if (verbose)
        ""
      else
        FALSE,
      stderr = if (verbose)
        ""
      else
        FALSE
    )
  }
  msdial_files <- list.files(msdir, "\\.msdial$", full.names = TRUE)
  if (readResults) {
    tryCatch({
      out <- vector("list", 2)
      names(out) <- c("aligned", "peaks")
      alignment_file <- grepl("AlignResult", msdial_files)
      if (any(alignment_file)) {
        fls <- msdial_files[which(alignment_file)]
        alignment_results <- lapply(fls, loadAlignmentResults, raw_folder = msdir)
        names(alignment_results) <- basename(sub("\\.msdial$", "", fls))
        out[[1]] <- alignment_results
      }
      single_files <- !alignment_file
      if (any(single_files)) {
        fls <- msdial_files[which(single_files)]
        single_results <- lapply(fls, loadConsoleResults)
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
