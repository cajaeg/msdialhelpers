#' Send mass spectrum to MSPepSearch
#'
#' Use "MSPepSearch64.exe" to search NIST library.
#' Download zip file from NIST site, copy all files to "MSSEARCH" directory and add the following line to your .Rprofile:
#' \code{options(msdialhelpers_mspepsearch_path = "/path/to/MSPepSearch64.exe")}
#'
#' Tested with MSPepSearch64 v2019_02_22.
#' @param spec mass spectrum: matrix or dataframe with m/z and intensity in the first two columns, or a list of matrices
#' @param searchtype one of "I" (identity, default), "Q" (quick identity), "S" (simple similarity),
#' "H" (hybrid similarity), "L" (neutral loss similarity), "M" (ms/ms in EI library),
#' "P" (peptide ms/ms search in a ms/ms library (use peak annotations and weighting), "G" (generic ms/ms search in a ms/ms library (no peak annotations or weighting used))
#' @param presearchtype one of "d" (standard (use pre-indexed peaks); default), "s" (sequential (compare to all library spectra)),
#' "f" (fast presearch (use pre-indexed peaks)), "m" (precursor ion m/z within precursor ion uncertainty - only for MS/MS search (P or G))
#' @param hits how many hits to return
#' @param extralibs libraries to use in addition to 'mainlib' and 'replib'; check "MSSEARCH" directory for further libraries that may be installed
#' @param ri_obs measured retention index (RI)
#' @param max_ri_dev passed to \code{scoreRI()}
#' @param ri_minscore minimum RI score; passed to \code{scoreRI()}
#' @param mzlimits set m/z limits for spectral comparison
#' @param localRdata attach meta data from local .rda (path to .rda or object already loaded in workspace)
#' @param .debug if TRUE keep temporary files created by MSPepSearch  
#' @return data.frame
#' @export
#'
#' @examples
#' \dontrun{
#' spec <- str2spec("135:999 107:275 150:262 41:150 95:142 136:103 77:84 39:74 91:72 65:49")
#' searchNIST(spec)
#' searchNIST(spec, ri_obs = 1300)
#' }
searchNIST <- function(spec,
                       searchtype = "I",
                       presearchtype = "d",
                       hits = 5,
                       extralibs = c(""),
                       ri_obs = NULL,
                       max_ri_dev = 100,
                       ri_minscore = 0.9,
                       mzlimits = c(-1, -1),
                       localRdata = getOption("msdialhelpers_localNISTrda"),
                       .debug = FALSE) {
  ri_dev = 5 # argument removed from argument list as it doesn't affect MSPepSearch scores; use max_ri_dev instead
  # previous RD text: allowed RI tolerance (1-32767)
  ri_penalty = 5 # argument removed from argument list as it doesn't have affect MSPepSearch scores; use ri_minscore instead
  # previous RD text: penalty rate (1-65535); NIST presets: Very weak(10), WK=Weak(20), AV=Average(50), ST=Strong(100), VS=Very strong(200), IN=Infinite(65535)
  mspepsearch <- getOption("msdialhelpers_mspepsearch_path")
  if (is.null(mspepsearch) || !file.exists(mspepsearch)) {
    stop(
      "Option \"msdialhelpers_mspepsearch_path\" not set or tool not found. See help file."
    )
  }
  nistmssearchdir <- dirname(mspepsearch)
  format_cas <- function(CAS) {
    CAS <- as.character(CAS)
    if (CAS == "0")
      return(NA_character_)
    sprintf(
      "%s-%s-%s",
      substr(CAS, 1, nchar(CAS) - 3),
      substr(CAS, nchar(CAS) - 2, nchar(CAS) - 1),
      substr(CAS, nchar(CAS), nchar(CAS))
    )
  }
  if (missing(spec)) {
    message("No spectrum provided, displaying MSPepSearch help page.")
    owd <- setwd(dirname(mspepsearch))
    on.exit(setwd(owd), add = TRUE)
    system2(command = basename(mspepsearch))
    return(invisible(NULL))
  }
  if (is.data.frame(spec) || is.matrix(spec)) {
    spec <- list(spec)
  }
  if (!is.null(ri_obs) && length(ri_obs) != length(spec))
    ri_obs <- rep(ri_obs[1], length(spec))
  stopifnot(extralibs == "" ||
              all(extralibs %in% tolower(basename(
                list.dirs(nistmssearchdir, recursive = FALSE)
              ))))
  stopifnot(nchar(searchtype) == 1 &&
              searchtype %in% c("I", "Q", "S", "H", "L", "M", "P", "G"))
  stopifnot(nchar(presearchtype) == 1 &&
              presearchtype %in% c("d", "s", "f", "m"))
  tmpdir <- tempdir()
  tmpdir_work <- file.path(tmpdir, "searchNIST")
  if (!.debug) {
    on.exit(unlink(tmpdir_work, recursive = TRUE, force = TRUE))
  }
  dir.create(tmpdir_work, recursive = TRUE)
  mspfile <- file.path(tmpdir_work, "tmp.msp")
  outfile <- file.path(tmpdir_work, "tmp.txt")
  if (file.exists(mspfile))
    unlink(mspfile, recursive = F, force = T)
  n_spec <- length(spec)
  mspNames <- sprintf("RQuery_%06d", 1:n_spec)
  showProgress <- n_spec > 100
  if (showProgress) {
    message("Exporting spectra to msp file ...")
    pb <- utils::txtProgressBar(min = 0,
                                max = length(mspNames),
                                style = 3)
  }
  zz_len <- 1e6 # size of pre-allocated character vector
  zz <- vector("character", zz_len)
  j = 1
  for (i in 1:length(spec)) {
    zz[j] <- sprintf("Name: %s", mspNames[i])
    j = j + 1
    if (!is.null(ri_obs))
      zz[j] <- sprintf("Retention_index: SemiStdNP=%d/1/4 StdNP=%d",
                       round(ri_obs[i]),
                       round(ri_obs[i]))
    j = j + 1
    m <- as.matrix(spec[[i]][, 1:2])
    zz[j] <- sprintf("Num Peaks: %d", nrow(m))
    j = j + 1
    tmp <- trimws(apply(m, 1, paste, collapse = " "))
    for (k in 1:length(tmp)) {
      zz[j] <- tmp[k]
      j = j + 1
    }
    zz[j] <- ""
    j = j + 1
    if (j > zz_len * 0.9) {
      cat(zz[1:j],
          file = mspfile,
          sep = "\n",
          append = file.exists(mspfile)) # space gets tight, write batch to disk
      zz <- vector("character", zz_len)
      j = 1
    }
    if (showProgress)
      utils::setTxtProgressBar(pb, i)
  }
  if (showProgress)
    close(pb)
  cat(zz[1:j],
      file = mspfile,
      sep = "\n",
      append = file.exists(mspfile))
  args <- c(
    sprintf("%s%s%s", presearchtype, "v", searchtype),
    # "v" allows RMF output with /OUTTAB
    "/MAIN mainlib",
    "/REPL replib",
    if (!is.null(ri_obs))
      sprintf("/RI sut%dr%d", ri_dev, ri_penalty),
    if (extralibs[1] != "")
      paste(paste("/LIB", extralibs), collapse = " "),
    sprintf("/INP %s", normalizePath(mspfile)),
    sprintf("/HITS %d", hits),
    sprintf("/MzLimits %d %d", mzlimits[1], mzlimits[2]),
    sprintf("/COL mw,nn,nm,cf,cn,ik"),
    sprintf("/OUTTAB %s", suppressWarnings(normalizePath(outfile)))
  )
  owd <- setwd(dirname(mspepsearch))
  on.exit(setwd(owd), add = TRUE)
  if (showProgress)
    message("Calling MSPepSearch ...")
  rv0 <- system2(
    command = basename(mspepsearch),
    args = args,
    stdout = TRUE,
    stderr = FALSE
  )
  rv <- iconv(rv0, "IBM860", "UTF-8")
  if (showProgress)
    message("Reading MSPepSearch results ...")
  txt0 <- readLines(outfile)
  ## browser()
  lidx <- grep("^>$", txt0) + c(2, -1)
  nr <- diff(lidx) + 1
  outcols <- c(
    "query",
    "rank",
    "name",
    "formula",
    "mf",
    "rmf",
    "prob",
    "ri_obs",
    "ri_lib",
    "ri_score",
    "total_score",
    "column_type",
    "match_no_ri",
    "cas",
    "mw",
    "exactmass",
    "lib",
    "lib_id",
    "nistrn",
    "inchikey",
    "peaks"
  )
  out <- data.frame(matrix(
    nrow = max(0, nr),
    ncol = length(outcols),
    dimnames = list(NULL, outcols)
  ),
  stringsAsFactors = FALSE)
  NISTcn <- strsplit(txt0[lidx[1] - 1], "\t")[[1]] # column names in NIST output
  hasRI <- length(grep("^RI$", NISTcn)) > 0
  if (nr > 0) {
    txt <- txt0[seq(lidx[1], lidx[2])]
    txt <- strsplit(txt, "\t")
    out$query <- as.numeric(sub("RQuery_([0-9]+).*", "\\1", sapply(txt, "[[", 1)))
    out$rank <- as.numeric(sapply(txt, "[[", grep("^Rank$", NISTcn)))
    out$name <- trimws(sapply(txt, "[[", grep("^Name$", NISTcn)))
    out$formula <- trimws(sapply(txt, "[[", grep("^Formula$", NISTcn)))
    out$mf <- as.numeric(sapply(txt, "[[", grep("^MF$", NISTcn)))
    out$rmf <- as.numeric(sapply(txt, "[[", grep("^R.Match$", NISTcn)))
    out$prob <- as.numeric(sapply(txt, "[[", grep("^Prob\\(", NISTcn)))
    if (!is.null(ri_obs))
      out$ri_obs <- round(rep(ri_obs, tapply(out$query, out$query, length)))
    if (hasRI) {
      out$ri_lib <- as.numeric(gsub("[^0-9]", "", sapply(txt, "[[", grep(
        "^RI$", NISTcn
      ))))
      out$ri_score <- round(
        scoreRI(
          out$ri_obs - out$ri_lib,
          max_ri_dev = max_ri_dev,
          minscore = ri_minscore
        ),
        2
      )
      out$total_score <- round(out$mf * out$ri_score)
      out$column_type <- sub(".*([A-Z])", "\\1", sapply(txt, "[[", grep("^RI$", NISTcn)))
      out$match_no_ri <- as.numeric(gsub("[^0-9]", "", sapply(
        txt, "[[", grep("^Match no RI$", NISTcn)
      )))
    }
    out$cas <- sapply(txt, function(x)
      format_cas(x[[grep("^CAS$", NISTcn)]]))
    out$mw <- as.numeric(sapply(txt, "[[", grep("^Lib MW$", NISTcn)))
    out$exactmass <- as.numeric(sapply(txt, "[[", grep("^Mass$", NISTcn)))
    out$lib <- trimws(sapply(txt, "[[", grep("^Library$", NISTcn)))
    out$lib_id <- as.numeric(sapply(txt, "[[", grep("^Id$", NISTcn)))
    out$nistrn <- as.numeric(sapply(txt, "[[", grep("^NIST r.n.$", NISTcn)))
    out$inchikey <- trimws(sapply(txt, "[[", grep("^InChIKey$", NISTcn)))
  }
  out <- out[order(out$query, if (hasRI)
    - out$total_score
    else
      - out$mf), ]
  out$rank <- as.numeric(unlist(tapply(out$query, out$query, seq_along)))
  if (!is.null(localRdata)) {
    if (inherits(localRdata, "data.frame")) {
      nistdb <- localRdata
    } else if (file.exists(localRdata)) {
      load(localRdata) # "nistdb"
    }
    if (exists("nistdb")) {
      idx <- match(out$nistrn, nistdb$nistrn)
      if (length(idx[is.finite(idx)]) > 0) {
        out$peaks[which(is.finite(idx))] <- nistdb$peaks[idx[is.finite(idx)]]
      }
    }
  }
  attr(out, "timeinfo") <- txt0[length(txt0)]
  return(out)
}
