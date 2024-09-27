#' Create a method file required by MS-DIAL console
#'
#' MS-DIAL parameters in console are set via a text file (option -m). Use this function to create that text file.
#' @param outfile name of text file to be created or 'NULL' to return method contents as list
#' @param data_type "Centroid" or "Profile"
#' @param ion_mode "Positive" or "Negative"
#' @param accuracy_type "IsAccurate" or "IsNominal"
#' @param retention_time_begin RT [min]
#' @param retention_time_end RT [min]
#' @param mass_range_begin m/z
#' @param mass_range_end m/z
#' @param number_of_threads 1
#' @param smoothing_method "LinearWeightedMovingAverage"
#' @param smoothing_level 2
#' @param average_peak_width 30
#' @param minimum_peak_height 1000
#' @param mass_slice_width 0.1 (for accurate mass)
#' @param mass_accuracy 0.025 (for accurate mass)
#' @param sigma_window_value 0.5
#' @param amplitude_cut_off 10
#' @param msp_file MSP file
#' @param ri_index_file ""
#' @param retention_type "RT" or "RI"
#' @param ri_compound "Alkanes" or "Fames"
#' @param retention_time_tolerance_for_identification 0.5
#' @param retention_index_tolerance_for_identification 20
#' @param ei_similarity_tolerance_for_identification 70
#' @param identification_score_cut_off 70
#' @param use_retention_information_for_scoring "True"
#' @param use_retention_information_for_filtering "False"
#' @param reference_file reference file for alignment
#' @param retention_time "RT" or "RI"
#' @param retention_time_tolerance_for_alignment 0.075
#' @param retention_index_tolerance_for_alignment 20
#' @param ei_similarity_tolerance 70
#' @param retention_time_factor 0.5
#' @param ei_similarity_factor 0.5
#' @param peak_count_filter 0
#' @param qc_at_least_filter "True"
#'
#' @return list or NULL depending on 'outfile'
#' @export
#'
#' @examples
#' createMSDCmethod() # print method contents to stdout
#' \dontrun{
#' createMSDCmethod(outfile = "test.txt") # write tag=value pairs
#' }
createMSDCmethod <- function(outfile = NULL,
                             ## Data type
                             data_type = "Centroid",
                             ion_mode = c("Positive", "Negative")[1],
                             accuracy_type = c("IsAccurate", "IsNominal")[1],
                             ## Data collection
                             retention_time_begin = 0,
                             retention_time_end = 100,
                             mass_range_begin = 0,
                             mass_range_end = 1000,
                             ## Data processing
                             number_of_threads = 1,
                             ## Peak detection
                             smoothing_method = "LinearWeightedMovingAverage",
                             smoothing_level = 2,
                             average_peak_width = 30,
                             minimum_peak_height = 1000,
                             mass_slice_width = 0.1,
                             mass_accuracy = 0.025,
                             ## MS1Dec parameters
                             sigma_window_value = 0.5,
                             amplitude_cut_off = 100,
                             ## Identification
                             msp_file = "",
                             ri_index_file = "",
                             retention_type = c("RT", "RI")[1],
                             ri_compound = c("Alkanes", "Fames")[1],
                             retention_time_tolerance_for_identification = 0.5,
                             retention_index_tolerance_for_identification = 20,
                             ei_similarity_tolerance_for_identification = 70,
                             identification_score_cut_off = 70,
                             use_retention_information_for_scoring = "True",
                             use_retention_information_for_filtering = "False",
                             ## Alignment
                             reference_file = "",
                             retention_time = c("RT", "RI")[1],
                             retention_time_tolerance_for_alignment = 0.15,
                             retention_index_tolerance_for_alignment = 200,
                             ei_similarity_tolerance = 70,
                             retention_time_factor = 0.5,
                             ei_similarity_factor = 0.5,
                             ## Filtering
                             peak_count_filter = 0,
                             qc_at_least_filter = "True") {
  method_contents <- list(
    ## Data type
    `Data type` = data_type,
    `Ion mode` = ion_mode,
    `Accuracy type` = accuracy_type,
    ## Data collection
    `Retention time begin` = retention_time_begin,
    `Retention time end` = retention_time_end,
    `Mass range begin` = mass_range_begin,
    `Mass range end` = mass_range_end,
    ## Data processing
    `Number of threads` = number_of_threads,
    ## Peak detection
    `Smoothing method` = smoothing_method,
    `Smoothing level` = smoothing_level,
    `Average peak width` = average_peak_width,
    `Minimum peak height` = minimum_peak_height,
    `Mass slice width` = mass_slice_width,
    `Mass accuracy` = mass_accuracy,
    ## MS1Dec parameters
    `Sigma window value` = sigma_window_value,
    `Amplitude cut off` = amplitude_cut_off,
    ## Identification
    `MSP file` = msp_file,
    `RI index file pathes` = if (identical(ri_index_file, ""))
      ""
    else
      "RI_index_file_paths.txt",
    `Retention type` = retention_type,
    `RI compound` = ri_compound,
    `Retention index tolerance for identification` = retention_index_tolerance_for_identification,
    `Retention time tolerance for identification` = retention_time_tolerance_for_identification,
    `EI similarity tolerance for identification` = ei_similarity_tolerance_for_identification,
    `Identification score cut off` = identification_score_cut_off,
    `Use retention information for scoring` = use_retention_information_for_scoring,
    `Use retention information for filtering` = use_retention_information_for_filtering,
    ## Alignment
    `Reference file` = reference_file,
    `Retention time` = retention_time,
    `Retention time tolerance for alignment` = retention_time_tolerance_for_alignment,
    `Retention index tolerance for alignment` = retention_index_tolerance_for_alignment,
    `EI similarity tolerance` = ei_similarity_tolerance,
    `Retention time factor` = retention_time_factor,
    `EI similarity factor` = ei_similarity_factor,
    ## Filtering
    `Peak count filter` = peak_count_filter,
    `QC at least filter` = qc_at_least_filter
  )
  if (!identical(ri_index_file, "")) {
    attr(method_contents, "ri_index_file_paths") <- "RI_index_file_paths.txt"
    attr(method_contents, "ri_index_file") <- ri_index_file
  }
  if (is.null(outfile)) {
    return(method_contents)
  } else {
    cat(sprintf("%s: %s", names(method_contents), unlist(method_contents)),
        sep = "\n",
        file = outfile)
    invisible(NULL)
  }
}
