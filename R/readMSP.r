#' Read- MSP library files
#'
#' @param mspfile file name
#' @param GMDsynonfix Fix for GMD MSP file (remove "Synon: " prefix)
#' @param .progress Show a progress bar, see \code{\link[plyr]{ldply}()}.
#'
#' @return data.frame
#' @export
#'
#' @examples
#' msp_file <- system.file("extdata/MSP_demo.msp", package = "msdialhelpers")
#' msp_lib <- readMSP(msp_file)
#' head(msp_lib)
readMSP <-
  function(mspfile,
           GMDsynonfix = FALSE,
           .progress = c("text", "none")[1])
  {
    if (!file.exists(mspfile))
      stop("Input file not found")
    txt <- readLines(mspfile)
    sec_start <- grep("^NAME: ", txt, ignore.case = T)
    sec_split <- rep(seq_along(sec_start), diff(c(sec_start, length(txt) +
                                                    1)))
    tmp <- split(txt, sec_split)
    .progress <- if (length(tmp) <= 50)
      "none"
    else
      .progress # disable PB for short MSP's
    plyr::ldply(tmp, function(x) {
      ## parse sections
      flt <- nchar(x) > 0 # remove blank lines
      x <- x[flt]
      is_metadata <- 1:length(x) <= grep("^num peaks:", x, ignore.case = T) # it's safe to assume that "Num Peaks:" comes last
      if (GMDsynonfix) {
        x[is_metadata] <- sub("^Synon: ", "", x[is_metadata])
      }
      metadata <- strsplit(x[is_metadata], ": ")
      df <- data.frame(
        matrix(
          sapply(metadata, "[", 2),
          nrow = 1,
          ncol = sum(is_metadata),
          dimnames = list(NULL, c(sapply(metadata, "[", 1)))
        ),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      convertsToNumeric <- which(!apply(df, 2, function(x)
        methods::is(tryCatch(
          as.numeric(x),
          warning = function(w)
            w
        ), "warning")))
      for (i in convertsToNumeric)
        df[, i] <- as.numeric(df[, i])
      pks <- x[!is_metadata]
      pks_str <- paste(gsub("\\t", ":", x[!is_metadata]), collapse = " ")
      df$PEAKS <- pks_str
      df
    }, .id = NULL, .progress = .progress)
  }
