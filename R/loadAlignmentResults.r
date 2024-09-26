#' Import alignment results from MS-DIAL
#'
#' Read "Height.txt" or "Area.txt" files obtained by running "Export ->
#' Alignment result" in MS-DIAL. Supports file format created by MS-DIAL
#' versions >= 4.0.
#' @param txt Alignment results file, i.e. "Height.txt" or "Area.txt"
#' @param raw_folder Location of raw data. If given, file paths of raw data will
#'   be attached to \code{attr(x, "msdial_sam")}.
#' @param raw_suffix File type of raw data. Adjust this when file type differs
#'   from 'mzML'. Only relevant when 'raw_folder' is given.
#' @param sample_names Optional character vector to use instead of original
#'   column names. No check for proper order is performed.
#'
#' @return data.frame with two attributes attached: "msdial_sam" and
#'   "intensity_columns".
#'
#'   \code{attr(x, "msdial_sam")} is derived from the file header and contains
#'   sample class information as defined in MS-DIAL.
#'
#'   \code{attr(x, "intensity_columns")} is an index of columns holding MS
#'   intensities (sample columns).
#' @export
#'
#' @examples
#' # basic use
#' fp <- system.file("extdata", package = "msdialhelpers")
#' txt_file <- file.path(fp, "MSDIAL_Alignment_result_GCQTOF.txt")
#' aligned <- loadAlignmentResults(txt_file)
#' head(aligned)
#'
#' # locate peak intensity columns
#' int_columns <- attr(aligned, "intensity_columns")
#' int_med <- apply(aligned[, int_columns], 1, median)
#' hist(log10(int_med))
#' 
#' # match raw mzML files to sample columns
#' aligned <- loadAlignmentResults(txt_file, raw_folder = fp)
#' msdial_sam <- attr(aligned, "msdial_sam")
#' colnames(aligned)[int_columns] == sub(".mzML", "", basename(msdial_sam$raw_file))
#' 
#' # add XICs for all peaks (e.g. for QC)
#' xraw <- loadRaws(msdial_sam$raw_file)
#' library(dplyr)
#' aligned2 <- aligned %>%
#'   as_tibble() %>%
#'   prepareSpectra() %>%
#'   addXICs(xraw = xraw)
#' 
#' flt <- 2 # plot compound of interest
#' plotXIC(aligned2$bpc[[flt]], main = aligned2$metabolite_name[flt])
#' plotXIC(aligned2$bpc[[flt]] %>% scaleXIC("by_file"), main = aligned2$metabolite_name[flt])
#' 
#' par(mfrow = n2mfrow(nrow(aligned2)), mar = c(2,2,1,0)+0.1, mgp = c(1.1, 0.1, 0), tck = 0.01)
#' plotXIC(lapply(aligned2$bpc, "[[", 1)) # QC panel for all compounds (file 1)
#'
#' # send spectra to MSPepSearch for identification
#' # needs MSPepSearch to be correctly configured, see ?MSPepSearch
#' \dontrun{
#' aligned3 <- aligned2 %>% 
#'   runNISTsearch(ri_column = NULL)
#' }
#'   
#' # dto., with retention index data
#' alkane_dict <- loadAlkanes(file.path(fp, "Alkane.txt"))
#' aligned4 <- aligned3 %>% 
#'   mutate(average_ri = rt2ri(average_rt_min, ref = alkane_dict)) %>% 
#'   runNISTsearch(ri_column = "average_ri")
#'
#' # compile spectral match results into new table
#' identified <- aligned4 %>% 
#'   group_by(alignment_id) %>% 
#'   reframe(NISTres[[1]])
#' head(identified)
loadAlignmentResults <-
  function(txt,
           raw_folder = NULL,
           raw_suffix = "mzml",
           sample_names = NULL) {
    stopifnot(file.exists(txt))
    hdr_nrow <- 4 # number of rows occupied by header
    ##
    ## evaluate header
    ##
    hdr <- utils::read.table(
      txt,
      sep = "\t",
      header = F,
      skip = 0,
      as.is = T,
      check.names = F,
      comment.char = "",
      na.strings = c("NA", "No record", "null"),
      quote = "",
      nrows = hdr_nrow + 1
    )
    ok <- "Class" %in% hdr[1, ]
    if (!ok) {
      stop("\"Class\" keyword not found in header, check file format")
    }
    intensity_columns_start <- grep("^Class$", hdr[1, ])[1] + 1
    intensity_columns_end <- max(which(!is.na(hdr[2, ]) &
                                         !duplicated(unlist(hdr[5, ]))))
    intensity_columns <- intensity_columns_start:intensity_columns_end
    msdial_sam <- data.frame(t(hdr[, intensity_columns]), stringsAsFactors = FALSE)
    cn <- hdr[, grep("Class", hdr[1, ])[1]]
    cn <- gsub("[^a-zA-Z0-9]+", "_", tolower(trimws(cn)))
    cn <- sub("_$", "", cn)
    cn <- sub("^_", "", cn)
    cn[5] <- "name"
    colnames(msdial_sam) <- cn
    rownames(msdial_sam) <- NULL
    msdial_sam <- msdial_sam[, c(grep("^name$", colnames(msdial_sam)), which(!grepl("^name$", colnames(msdial_sam))))]
    other_columns <- 1:(intensity_columns[1] - 1)
    ##
    ## evaluate body
    ##
    out <- utils::read.table(
      txt,
      sep = "\t",
      header = T,
      skip = hdr_nrow,
      as.is = T,
      check.names = F,
      comment.char = "",
      na.strings = c("NA", "No record", "null"),
      quote = ""
    )
    cn <- colnames(out)[other_columns]
    cn <- make.names(cn)
    cn <- gsub("[^a-zA-Z0-9]+", "_", tolower(trimws(cn))) # fix column names
    cn <- sub("_$", "", cn)
    cn <- sub("^_", "", cn)
    msd_mode <- if ("ei_spectrum" %in% cn)
      "GCMS"
    else
      "LCMS"
    if (msd_mode == "GCMS") {
      cn[cn == "ei_spectrum"] <- "ms1"
    } else {
      cn[cn == "ms1_isotopic_spectrum"] <- "ms1"
      cn[cn == "ms_ms_spectrum"] <- "ms2"
    }
    colnames(out)[other_columns] <- cn
    ## msdial_sam$name <- colnames(out)[intensity_columns] # include sample column names in 'msdial_sam'
    ## msdial_sam <- msdial_sam[,c(5,1:4)]
    for (i in which(apply(out, 2, function(x)
      all(x %in% c("False", "True"))))) {
      out[, i] <- as.logical(out[, i]) # encode as logical where appropriate
    }
    if (!is.null(sample_names)) {
      if (length(sample_names) == length(intensity_columns)) {
        colnames(out)[intensity_columns] <- sample_names
      } else {
        stop("supplied sample names do not match number of sample columns")
      }
    }
    if (ncol(out) > max(intensity_columns)) {
      out <- out[, -seq(max(intensity_columns) + 1, ncol(out))] # remove Average & SD columns
    }
    if (!is.null(raw_folder) && dir.exists(file.path(raw_folder))) {
      raw_suffix_pat <- sprintf("\\.%s$", raw_suffix)
      raw_files <- list.files(raw_folder,
                              raw_suffix_pat,
                              ignore.case = TRUE,
                              full.names = TRUE)
      pat <- sub(raw_suffix_pat, "", basename(raw_files), ignore.case = TRUE)
      o <- match(msdial_sam$name, pat)
      raw_files <- raw_files[o]
      msdial_sam$raw_file <- raw_files
    } else {
      msdial_sam$raw_file <- NA_character_
    }
    attr(out, "msdial_sam") <- msdial_sam
    attr(out, "intensity_columns") <- intensity_columns
    message(sprintf(
      "%s data loaded for %d compounds x %d samples",
      msd_mode,
      nrow(out),
      length(intensity_columns)
    ))
    return(out)
  }
