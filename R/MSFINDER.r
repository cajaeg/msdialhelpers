#' Create MSFINDER query
#'
#' Export mass spectrum to 'MAT' file that can be opened in MSFINDER.
#' @param name Give the spectrum a name
#' @param precursormz Precursor m/z
#' @param precursortype Precursor type
#' @param ionmode "Positive" or "Negative"
#' @param formula Sum formula (required for fragment annotation)
#' @param SMILES SMILES string (required for fragment annotation)
#' @param collision_energy Collision energy (optional)
#' @param rt Retention time of the spectrum that will be used by MSFINDER for refined prediction (optional)
#' @param ms1spec MS1 spectrum (optional)
#' @param ms2spec MS2 spectrum (required)

#' @param outfile Name of MAT file, or \code{NULL} for 'stdout'.
#'
#' @return NULL
#' @export
#' @keywords internal
writeMAT <- function(name = "unknown",
                     precursormz = NULL,
                     precursortype = "[M+H]+",
                     ionmode = c("Positive", "Negative", "Both")[1],
                     formula = "",
                     SMILES = "",
                     collision_energy = 40,
                     rt = 0,
                     ms1spec = NULL,
                     ms2spec = NULL,
                     outfile = NULL) {
  if (!is.null(outfile)) {
    outdir <- dirname(outfile)
    if (!dir.exists(outdir)) {
      dir.create(outdir, recursive = TRUE)
    }
  }
  "%&%" <- function(x, y)
    paste0(x, y)
  spectra <- ""
  if (!is.null(ms1spec) &&
      nrow(ms1spec) > 0 &&
      any(ms1spec[, 2] > 0) && any(ms1spec[, 2] > 0)) {
    spectra <- spectra %&% "MSTYPE: MS1\n"
    spectra <- spectra %&% "Num Peaks: " %&% nrow(ms1spec[ms1spec[, 2] >
                                                            0, , drop = F]) %&% "\n"
    for (i in which(ms1spec[, 2] > 0)) {
      spectra <- spectra %&% ms1spec[i, 1] %&% " " %&% ms1spec[i, 2] %&% "\n"
    }
  }
  if (nrow(ms2spec) > 0 &&
      any(!is.na(ms2spec[, 2])) && any(ms2spec[, 2] > 0)) {
    spectra <- spectra %&% "MSTYPE: MS2\n"
    spectra <- spectra %&% "Num Peaks: " %&% nrow(ms2spec[ms2spec[, 2] >
                                                            0, , drop = F]) %&% "\n"
    for (i in which(ms2spec[, 2] > 0)) {
      spectra <- spectra %&% ms2spec[i, 1] %&% " " %&% ms2spec[i, 2] %&% "\n"
    }
  }
  out <-
    sprintf(
      "NAME: %s\nPRECURSORMZ: %f\nPRECURSORTYPE: %s\nIONMODE: %s\nFORMULA: %s\nSMILES: %s\nCOLLISIONENERGY: %d\nRETENTIONTIME: %f\n%s",
      as.character(name),
      precursormz,
      as.character(precursortype),
      as.character(ionmode),
      if (is.null(formula))
        ""
      else
        as.character(formula),
      if (is.null(SMILES))
        ""
      else
        as.character(SMILES),
      as.numeric(collision_energy),
      as.numeric(rt),
      spectra
    )
  # browser()
  cat(out, file = if (is.null(outfile))
    ""
    else
      outfile)
}


#' Trigger MS-FINDER process (desktop app)
#'
#' @param matpath Path containing .mat files.
#' @param MSFINDER_path Full path of 'MSFINDER.exe' file. Uses option
#'   'msdialhelpers_msfinder_exe' if set. Add
#'   \code{options(msdialhelpers_msfinder_exe = "/path/to/MSFINDER.exe")} or
#'   similar to '~/.Rprofile' to customize. Use \code{Sys.getenv("R_USER")} to
#'   check where this file is located on Windows.
#'
#' @return \code{TRUE} if process was spawned successfully
#' @export
#'
#' @examples
#' \dontrun{
#' matpath <- "~/MATtmp"
#' dir.create(matpath)
#' file.copy(list.files(system.file("extdata/MSFINDER", package = "msdialhelpers"),
#'   "*.mat",
#'   full = TRUE), matpath)
#' triggerMSF(matpath) # works only if MSFINDER_path is set correctly
#' }
triggerMSF <- function(matpath,
                       MSFINDER_path = getOption("msdialhelpers_msfinder_exe")) {
  if (is.null(MSFINDER_path) || !file.exists(MSFINDER_path))
    stop("MSFINDER executable not found")
  rv <- system2(normalizePath(MSFINDER_path),
                normalizePath(matpath),
                wait =
                  FALSE)
  cat("Now start 'Compound annotation (batch job)' from MSFINDER menu\n")
  invisible(rv)
}


#' Generate .ini file needed by MSFINDER console app
#'
#' @return list
#' @export
#'
#' @examples
#' ini <- getMSFINDERini()
#' ini$Ms1Tolerance <- 0.01 # customize parameters
#' \dontrun{
#' cat(paste(names(ini), ini, sep="=", collapse="\n"), file="MSFINDER.INI")
#' }
getMSFINDERini <- function() {
  list(
    #Formula finder parameters
    LewisAndSeniorCheck = TRUE,
    Ms1Tolerance = 0.01,
    IsotopicAbundanceTolerance = 20,
    MassToleranceType = "Da",
    CommonRange = TRUE,
    ExtendedRange = FALSE,
    ExtremeRange = FALSE,
    ElementProbabilityCheck = FALSE,
    Ocheck = TRUE,
    Ncheck = TRUE,
    Pcheck = TRUE,
    Scheck = TRUE,
    Fcheck = FALSE,
    ClCheck = FALSE,
    BrCheck = FALSE,
    Icheck = FALSE,
    SiCheck = FALSE,
    IsTmsMeoxDerivative = FALSE,
    MinimumTmsCount = 1,
    MinimumMeoxCount = 0,
    CanExcuteMS1AdductSearch = FALSE,
    CanExcuteMS2AdductSearch = FALSE,
    FormulaMaximumReportNumber = 100,
    #Structure finder parameters
    TreeDepth = 2,
    Ms2Tolerance = 0.05,
    RelativeAbundanceCutOff = 1,
    StructureMaximumReportNumber = 100,
    IsUseEiFragmentDB = FALSE,
    #Data source
    MinesNeverUse = TRUE,
    MinesOnlyUseForNecessary = FALSE,
    MinesAllways = FALSE,
    PubChemNeverUse = TRUE,
    PubChemOnlyUseForNecessary = FALSE,
    PubChemAllways = FALSE,
    HMDB = TRUE,
    YMDB = TRUE,
    PubChem = TRUE,
    SMPDB = TRUE,
    UNPD = TRUE,
    ChEBI = TRUE,
    PlantCyc = TRUE,
    KNApSAcK = TRUE,
    BMDB = TRUE,
    FooDB = TRUE,
    ECMDB = TRUE,
    DrugBank = TRUE,
    T3DB = TRUE,
    STOFF = TRUE,
    NANPDB = TRUE,
    LipidMAPS = TRUE,
    Urine = TRUE,
    Saliva = TRUE,
    Feces = TRUE,
    Serum = TRUE,
    Csf = TRUE,
    IsUserDefinedDB = TRUE,
    UserDefinedDbFilePath = "",
    #Spectral database search
    IsRunSpectralDbSearch = FALSE,
    IsRunInSilicoFragmenterSearch = TRUE,
    IsPrecursorOrientedSearch = TRUE,
    IsUseInternalExperimentalSpectralDb = FALSE,
    IsUseInSilicoSpectralDbForLipids = FALSE,
    IsUseUserDefinedSpectralDb = FALSE,
    UserDefinedSpectralDbFilePath = "",
    SolventType = "CH3COONH4",
    MassRangeMin = 0,
    MassRangeMax = 2000,
    #Retention time setting for structure elucidation
    IsUsePredictedRtForStructureElucidation = FALSE,
    IsUseRtInchikeyLibrary = TRUE,
    IsUseXlogpPrediction = FALSE,
    RtInChIKeyDictionaryFilepath = "",
    RtSmilesDictionaryFilepath = "",
    Coeff_RtPrediction = -1,
    Intercept_RtPrediction = -1,
    RtToleranceForStructureElucidation = 2.5,
    RtPredictionSummaryReport = "",
    IsUseRtForFilteringCandidates = FALSE,
    #Retention time setting for spectral searching
    IsUseExperimentalRtForSpectralSearching = FALSE,
    RetentionType = "RT",
    RtToleranceForSpectralSearching = 0.5,
    #Batch job
    AllProcess = TRUE,
    FormulaFinder = TRUE,
    StructureFinder = TRUE,
    TryTopNMolecularFormulaSearch = 5,
    #FSEA parameter
    FseaRelativeAbundanceCutOff = 5,
    FseanonsignificantDef = "OntologySpace",
    FseaPvalueCutOff = 1,
    #Msfinder molecular networking (mmn)
    IsMmnLocalCytoscape = TRUE,
    IsMmnMsdialOutput = FALSE,
    IsMmnFormulaBioreaction = FALSE,
    IsMmnRetentionRestrictionUsed = FALSE,
    IsMmnOntologySimilarityUsed = TRUE,
    MmnMassTolerance = 0.025,
    MmnRelativeCutoff = 1,
    MmnMassSimilarityCutOff = 75,
    MmnRtTolerance = 100,
    MmnOntologySimilarityCutOff = 90,
    MmnOutputFolderPath = "",
    #Time out parameter
    FormulaPredictionTimeOut = -1,
    StructurePredictionTimeOut = -1
  )
}


#' Trigger MSFINDER console process
#'
#' @param inputFolder Location of MAT/MSP files to be processed. Required for 'predict' and 'annotate'.
#' @param inputFile Name of text file containing SMILES codes (required for 'generate').
#' @param outputFolder Where to place results files. By default, the same as 'inputFolder' or base folder of 'inputFile', respectively.
#' @param analysisType MSFINDER mode: 'predict', 'annotate' or 'generate'.
#' @param MSFINDER_console Full path of 'MsfinderConsoleApp.exe' file. Add \code{options(msdialhelpers_msfinder_console = "/path/to/MsfinderConsoleApp.exe")} to '~/.Rprofile' to customize. Use \code{Sys.getenv("R_USER")} to check where this file is located on Windows.
#' @param MSFINDER_ini List representation of MSFINDER.INI as returned by \code{\link{getMSFINDERini}()}
#' @param readResults read back MSFINDER results (using \code{readFGT()})
#' @param verbose show MSFINDER messages
#'
#' @return list of MSFINDER results or NULL depending on \code{results}
#' @export
#'
#' @examples
#' \dontrun{
#' if(!is.null(getOption("msdialhelpers_msfinder_console"))) {
#'   ##
#'   ## analysisType == "predict"
#'   ##
#'   dir.create(inputFolder <- tempfile(pattern = "RMSFINDER"))
#'   file.copy(list.files(system.file("extdata/MSFINDER", package = "msdialhelpers"), 
#'     "*.mat", 
#'     full.names = T), inputFolder)
#'   triggerMSFconsole(inputFolder)
#'   fgtfiles <- list.files(inputFolder, "*.fgt", full.names = T) # locate files produced by MSFINDER
#'   res <- readFGT(fgtfiles[1])
#'   summarizeMSFresults(res)
#'   unlink(inputFolder, recursive = T)
#'   ##
#'   ## analysisType == "annotate"
#'   ##
#'   dir.create(inputFolder <- tempfile(pattern = "RMSFINDER"))
#'   file.copy(list.files(system.file("extdata/MSFINDER", package = "msdialhelpers"), 
#'     "*.mat", 
#'     full.names = T), inputFolder)
#'   triggerMSFconsole(inputFolder, analysisType = "annotate")
#'   unlink(inputFolder, recursive = T)
#'   ##
#'   ## analysisType == "generate"
#'   ##
#'   dir.create(inputFolder <- tempfile(pattern = "RMSFINDER"))
#'   smiles <- c("N[C@@H](CCC(N)=O)C(O)=O", "NCCS(O)=O", 
#'     "O[C@@H]1[C@@H](COP(O)(O)=O)O[C@H]([C@@H]1O)N1C=NC2=C1N=CNC2=O")
#'   cat(smiles, file = inputFile <- file.path(inputFolder, "smiles.txt"), sep = "\n")
#'   triggerMSFconsole(inputFile = inputFile, analysisType = "generate")
#'   list.files(dirname(inputFile))
#'   unlink(inputFolder, recursive = TRUE)
#' }
#' }
triggerMSFconsole <- function(inputFolder = NULL,
                              inputFile = NULL,
                              outputFolder = NULL,
                              analysisType = c("predict", "annotate", "generate")[1],
                              MSFINDER_console = getOption("msdialhelpers_msfinder_console"),
                              MSFINDER_ini = getMSFINDERini(),
                              readResults = TRUE,
                              verbose = TRUE) {
  ## check arguments
  if (is.null(MSFINDER_console) || !file.exists(MSFINDER_console))
    stop("MS-FINDER console executable not found")
  
  analysisType <- match.arg(analysisType, choices = c("predict", "annotate", "generate"))
  if (analysisType %in% c("predict", "annotate")) {
    if (is.null(inputFolder) || !dir.exists(inputFolder)) {
      stop("'inputFolder' required for this analysis type")
    }
    inputFolder <- normalizePath(inputFolder)
    if (is.null(outputFolder))
      outputFolder <- inputFolder
    else
      dir.create(outputFolder)
  } else if (analysisType == "generate") {
    if (is.null(inputFile) || !file.exists(inputFile)) {
      stop("'inputFile' required for this analysis type")
    }
    inputFolder <- normalizePath(dirname(inputFile))
    if (is.null(outputFolder))
      outputFolder <- inputFolder
    else
      dir.create(outputFolder)
    outputFile <- file.path(outputFolder,
                            sub(
                              "(.*)\\.([A-Za-z]+)$",
                              "\\1_out.\\2",
                              basename(inputFile)
                            ))
  }
  
  msFinderConsoleApp <- normalizePath(MSFINDER_console)
  methodFile <- file.path(inputFolder, "MSFINDER.INI")
  cat(paste(
    names(MSFINDER_ini),
    MSFINDER_ini,
    sep = "=",
    collapse = "\n"
  ), file = methodFile)
  
  ## batchScriptContents <- switch(analysisType,
  ##    predict = sprintf('"%s" %s -i "%s" -o "%s" -m "%s"', msFinderConsoleApp, analysisType, inputFolder, outputFolder, normalizePath(methodFile)),
  ##    annotate = sprintf('"%s" %s -i "%s" -m "%s"', msFinderConsoleApp, analysisType, inputFolder, normalizePath(methodFile)),
  ##    generate = sprintf('"%s" %s -i "%s" -o "%s"', msFinderConsoleApp, analysisType, normalizePath(inputFile), suppressWarnings(normalizePath(outputFile)))
  ## )
  ## batchScriptFile <- file.path(inputFolder, "run.bat")
  ## cat(batchScriptContents, file=batchScriptFile)
  ## shell.exec(batchScriptFile)
  
  system2(
    utils::shortPathName(msFinderConsoleApp),
    args = switch(
      analysisType,
      predict = c(
        analysisType,
        paste0("-i ", utils::shortPathName(inputFolder)),
        paste0("-o ", utils::shortPathName(outputFolder)),
        paste0("-m ", utils::shortPathName(normalizePath(methodFile)))
      ),
      annotate = c(
        analysisType,
        paste0("-i ", utils::shortPathName(inputFolder)),
        paste0("-m ", utils::shortPathName(normalizePath(methodFile)))
      ),
      generate = c(
        analysisType,
        paste0("-i ", utils::shortPathName(normalizePath(inputFile))),
        paste0("-o ", utils::shortPathName((
          suppressWarnings(normalizePath(outputFile))
        )))
      )
    ),
    stdout = if (verbose)
      ""
    else
      FALSE
  )
  
  if (readResults) {
    if (analysisType == "predict") {
      fgt <- list.files(outputFolder, ".*\\.fgt", full.names = TRUE)
      res <- lapply(fgt, readFGT)
      names(res) <- basename(fgt)
    } else if (analysisType == "annotate") {
      res <- NULL ## TBD
    } else if (analysisType == "") {
      res <- NULL ## TBD
    }
  } else {
    res <- NULL
  }
  
  invisible(res)
}


#' Read MSFINDER result files
#'
#' @param file path to FGT files which MSFINDER puts in same folder as input MAT files
#'
#' @return List containing MSFINDER results. First element is match 1, second is match 2 etc. Each element contains data.frames 'cmpd', 'production' and 'neutralloss'. In case structures were found, it further contains further list elements named 'structure*'. Use \code{summarizeMSFresults} to extract the essential parts.
#' @export
#'
#' @examples
#' fgt <- list.files(system.file("extdata/MSFINDER", package="msdialhelpers"), 
#'   "*.fgt", 
#'   full.names = TRUE)
#' res <- lapply(fgt, readFGT)
readFGT <-
  function (file) {
    get.text.value <- function(x, field) {
      l <- strsplit(x, field, fixed = TRUE)
      sapply(l, function(x) {
        if (length(x) == 2) {
          gsub("^ +", "", x[2])
        } else {
          NA
        }
      })
    }
    convert.to.numeric <- function(x) {
      numericFields <- if (methods::is(x, "data.frame")) {
        which(apply(x, 2, function(y)
          ! any(is.na(
            suppressWarnings(as.numeric(y))
          ))))
      } else {
        which(sapply(x, function(y)
          ! any(is.na(
            suppressWarnings(as.numeric(y))
          ))))
      }
      for (i in numericFields)
        x[[i]] <- as.numeric(x[[i]]) # convert
      return(x)
    }
    read.compound <- function(strs) {
      # fields.idx <- grep("^[[:alnum:] \\(\\)]+:.*", strs)
      fields.idx <- grep("^[[:alnum:]][^:]+:.*", strs)
      fields <- sapply(strsplit(strs[fields.idx], ":"), "[[", 1)
      pk.idx <- grep("^Num ", fields) # MSP-like peak sections (ProductIon, NeutralLoss)
      npeaks <- as.numeric(sapply(pk.idx, function(i)
        get.text.value(strs[fields.idx[i]], paste0(fields[i], ":")))) # number of peaks in the different peak sections
      npkSections <- length(pk.idx) # number of peak sections
      peaks.idx <- unlist(sapply(1:npkSections, function(i)
        (fields.idx[pk.idx[i]] + 1):(fields.idx[pk.idx[i]] + npeaks[i]))) # all peak lines (don't process here)
      cmpd <- lapply(fields.idx[-peaks.idx], function(i)
        get.text.value(strs[i], paste0(fields[i], ":")))
      names(cmpd) <- tolower(fields[-peaks.idx])
      cmpd <- convert.to.numeric(cmpd)
      pk <- lapply(1:npkSections, function(i) {
        # parse peak sections
        ColNames <- tolower(unlist(strsplit(
          sub(".*\\(([A-Za-z\\/ ]+)\\).*", "\\1", fields[pk.idx[i]]),
          " ",
          TRUE
        )))
        out <- data.frame(matrix(
          ncol = length(ColNames),
          nrow = npeaks[i],
          dimnames = list(NULL, ColNames)
        ))
        if (npeaks[i] > 0) {
          peaks.idx <- (fields.idx[pk.idx[i]] + 1):(fields.idx[pk.idx[i]] + npeaks[i])
          peakrows <- strsplit(strs[peaks.idx], "\t")
          peakrows <- lapply(peakrows, "[", 1:length(ColNames)) # ensure equal lengths
          for (j in 1:length(peakrows))
            out[j, ] <- peakrows[[j]]
          out <- convert.to.numeric(out)
        }
        return(out)
      })
      names(pk) <- tolower(sub("Num ([A-Za-z]+).*", "\\1", fields[pk.idx]))
      c(list(cmpd = as.data.frame(cmpd, stringsAsFactors = FALSE)), pk)
    }
    txt <- scan(file,
                what = "",
                sep = "\n",
                quiet = TRUE)
    starts <- which(regexpr("^NAME: ", txt, ignore.case = TRUE) == 1)
    ends <- c(starts[-1] - 1, length(txt))
    if (length(starts) > 0) {
      out <- lapply(1:length(starts), function(i) {
        read.compound(txt[starts[i]:ends[i]])
      })
      SF <- sapply(out, function(x)
        x$cmpd$name)
      StrucPaths <- file.path(dirname(file),
                              sub("\\.fgt", "", basename(file)),
                              paste0(SF, ".sfd"))
      hasStruc <- file.exists(StrucPaths)
      if (any(hasStruc)) {
        Struc <- lapply(StrucPaths[hasStruc], readFGT)
        for (i in which(hasStruc)) {
          if (!is.null(Struc[[i]]))
            out[[i]] <- c(out[[i]], structure = Struc[[i]][order(-sapply(Struc[[i]], function(x)
              x[[1]]$totalscore))])
        }
      }
      return(out)
    } else {
      NULL
    }
  }


#' Extract essential information from full MSFINDER results
#'
#' @param x Full result set as returned by \code{readFGT} (list)
#' @param rank Which sum formula ranks should be returned (numeric)
#'
#' @return Data.frame
#' @export
#'
#' @examples
#' fgt <- list.files(system.file("extdata/MSFINDER", package = "msdialhelpers"), 
#'   "*.fgt", 
#'   full.names = TRUE)[1]
#' full <- readFGT(fgt)
#' summarizeMSFresults(full)
#'
#' fgt <- list.files(system.file("extdata/MSFINDER", package = "msdialhelpers"), 
#'   "*.fgt", 
#'   full.names = TRUE)
#' full <- lapply(fgt, readFGT)
#' names(full) <- 1:length(full) # make use of plyr::ldply .id column
#' plyr::ldply(full, summarizeMSFresults, rank = 1:2)
summarizeMSFresults <- function(x, rank = 1) {
  if (max(rank) > length(x)) {
    warning("Result set shorter than requested number of ranks")
    rank <- rank[-which(rank > length(x))]
  }
  plyr::ldply(rank, function(i) {
    tmp <- x[[i]]$cmpd # sum formula info
    useColnames <- c("name",
                     "exactmass",
                     "accuratemass",
                     "massdifference",
                     "totalscore")
    out1 <- if (!is.null(tmp))
      tmp[, useColnames]
    else
      data.frame(matrix(
        NA,
        nrow = 1,
        ncol = length(useColnames),
        dimnames = list(NULL, useColnames)
      ),
      stringsAsFactors = F)
    colnames(out1)[1] <- "formula" # avoid duplicate name
    colnames(out1)[5] <- "totalscore_formula"
    tmp <- x[[i]]$structure1 # structure info
    useColnames <- c("name",
                     "id",
                     "inchikey",
                     "smiles",
                     "ontology",
                     "totalscore")
    out2 <- if (!is.null(tmp))
      tmp$cmpd[, useColnames]
    else
      data.frame(matrix(
        NA,
        nrow = 1,
        ncol = length(useColnames),
        dimnames = list(NULL, useColnames)
      ),
      stringsAsFactors = F)
    cbind(out1, out2, rank = i, stringsAsFactors = F)
  }, .id = NULL)
}

#' Read back MSFINDER MAT files
#'
#' Read back MSFINDER query files.
#'
#' @param x path to .mat file (character)
#' @param digits number of decimal places to round mz and intensity to, resp. (numeric(2))
#'
#' @return data.frame
#' @export
#'
#' @examples
#' matfile <- list.files(system.file("extdata/MSFINDER", package = "msdialhelpers"), 
#'   "*.mat", 
#'   full.names = TRUE)[1]
#' readMAT(matfile)
#'
#' matfiles <- list.files(system.file("extdata/MSFINDER", package = "msdialhelpers"), 
#'   "*.mat", 
#'   full.names = TRUE)
#' names(matfiles) <- basename(matfiles)
#' plyr::ldply(matfiles, readMAT)
readMAT <- function(x, digits = c(4, 0)) {
  x <- readLines(x)
  meta <- 1:(grep("^MSTYPE", x)[1] - 1) # first process meta data
  tmp <- strsplit(x[meta], ": ")
  Names <- sapply(tmp, "[", 1)
  Vals <- sapply(tmp, "[", 2)
  out1 <- data.frame(t(matrix(Vals, dimnames = list(Names, NULL))), stringsAsFactors = F)
  isNumeric <- suppressWarnings(!is.na(apply(out1, 2, as.numeric))) # which columns can be converted numeric
  out1[, isNumeric] <- as.numeric(out1[, isNumeric])
  x <- x[-meta] # done with meta data, remove
  nspec <- length(grep("^MSTYPE", x))
  seclen <- diff(c(grep("^MSTYPE", x), length(x) + 1))
  splitvec <- rep(1:nspec, seclen)
  secs <- split(x, splitvec)
  process_sec <- function(x) {
    x <- x[grep("^[0-9\\.]+", x)] # discard non-peaks
    x <- strsplit(x, " ")
    fmtfun <- function(x)
      paste(round(as.numeric(x[1]), digits[1]),
            round(as.numeric(x[2]), digits[2]),
            sep = ":")
    paste(sapply(x, fmtfun), collapse = " ")
  }
  out2 <- data.frame(t(matrix(sapply(secs, process_sec))), stringsAsFactors = F)
  colnames(out2) <- sapply(secs, function(x)
    sub("MSTYPE: ", "", x[grep("^MSTYPE", x)]))
  out <- cbind(out1, out2, stringsAsFactors = F)
  return(out)
}
