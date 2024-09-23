#' Annotate MS/MS fragments with MS-FINDER
#'
#' TBD
#' @param x data.frame
#' @param in_column name of column containing spectra to annotate
#' @param id_column name of column containing compound IDs
#' @param formula_column name of column containing compound sum formulas
#' @param smiles_column name of column containing compound SMILES codes
#' @param precursormz_column name of column containing precursor m/z's
#' @param precursortype_column name of column containing precursor types. Known precursor types include "\[M+H\]+", "\[M+Na\]+", "\[M-H\]-", "\[M\]+." etc.
#' @param ini list of MSFINDER parameters as returned by \code{getMSFINDERini()}
#' @param ionmode "Positive" or "Negative"
#' @param collision_energy in eV
#' @param matbase working directory
#' @param skipRun don't run MSFINDER but only read results
#' @param keepFolder delete 'matbase' folder after finishing?
#' @param verbose passed to \code{\link{triggerMSFconsole}()}
#' @return data.frame
#' @export
#' @importFrom rlang .data
#' @importFrom rlang :=
annotateFragmentsMSFINDER <- function(x,
                                      in_column = "s",
                                      id_column = "alignment_id",
                                      formula_column = "formula",
                                      smiles_column = "smiles",
                                      precursormz_column = "precursor_mz",
                                      precursortype_column = "precursor_type",
                                      ini = NULL,
                                      ionmode = c("Positive", "Negative")[1],
                                      collision_energy = 20,
                                      matbase = "./MSFINDER_tmp",
                                      skipRun = FALSE,
                                      keepFolder = FALSE,
                                      verbose = TRUE) {
  if (is.null(ini)) {
    init <- getMSFINDERini()
  }
  ## create and populate temporary folder for MSFINDER
  MSFin <-
    data.frame(
      id = x[[id_column]],
      formula = x[[formula_column]],
      smiles = x[[smiles_column]],
      precursormz = x[[precursormz_column]],
      precursortype = x[[precursortype_column]],
      stringsAsFactors = FALSE
    ) |>
    tibble::as_tibble() |>
    dplyr::group_by(.data$id) |>
    dplyr::mutate("cand" := 1:dplyr::n()) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(.data$formula) & !is.na(.data$smiles)) |>
    dplyr::mutate("matname" := sprintf("id%06d_cand%06d", .data$id, .data$cand)) |>
    dplyr::mutate("matfilename" := file.path(.data$matbase, paste0(.data$matname, ".mat"))) |>
    dplyr::mutate("ionmode" := ionmode) |>
    dplyr::mutate("collision_energy" := collision_energy)
  MSFin$s <- x[[in_column]][match(MSFin$id, x[[id_column]])]
  if (!skipRun) {
    if (dir.exists(matbase)) {
      unlink(matbase, recursive = TRUE, force = TRUE)
    } else {
      dir.create(matbase)
    }
    if (!keepFolder) {
      on.exit(unlink(matbase, recursive = TRUE, force = TRUE))
    }
    message(
      sprintf(
        "Creating %d .mat files (%d spectra each with up to %d annotation hypotheses)",
        nrow(MSFin),
        length(unique(MSFin$id)),
        max(MSFin$cand)
      )
    )
    pb <- utils::txtProgressBar(max = nrow(MSFin))
    for (i in 1:nrow(MSFin)) {
      utils::setTxtProgressBar(pb, i)
      msf <- MSFin[i, ]
      writeMAT(
        name = msf$matname,
        precursormz = msf$precursormz,
        precursortype = msf$precursortype,
        ionmode = msf$ionmode,
        formula = msf$formula,
        SMILES = msf$smiles,
        collision_energy = msf$collision_energy,
        ms1spec = msf$s[[1]][, 1:2],
        ms2spec = msf$s[[1]][, 1:2],
        outfile = msf$matfilename
      )
    }
    close(pb)
    ## trigger MSFINDER process
    message("Triggering MSFINDER process")
    triggerMSFconsole(
      utils::shortPathName(normalizePath(matbase)),
      analysisType = "annotate",
      MSFINDER_ini = ini,
      readResults = FALSE,
      verbose = verbose
    )
  }
  ## read MSFINDER results
  fgtfiles <- list.files(matbase, "\\.fgt$", full.names = TRUE)
  message(sprintf(
    "Reading back %d .fgt files and attaching results to output",
    length(fgtfiles)
  ))
  if (length(fgtfiles) != nrow(MSFin)) {
    warning("MSFINDER results do not match number of input spectra")
  }
  fgt0 <- lapply(fgtfiles, readFGT)
  fgt <- lapply(fgt0, function(x)
    x[[1]]$production[, c(3, 4, 1, 2, 5)])
  MSFres <- data.frame(fgtfile = fgtfiles) |>
    tibble::as_tibble() |>
    dplyr::mutate("id" := as.integer(sub(
      "id([0-9]+)_.*", "\\1", basename(.data$fgtfile)
    )), cand = as.integer(sub(
      ".*_cand([0-9]+)[^0-9]+", "\\1", basename(.data$fgtfile)
    ))) |>
    dplyr::mutate("fgt" := .data$fgt) |>
    dplyr::left_join(MSFin[, c("id", "cand", "s")], by = c("id", "cand")) |>
    dplyr::rename("s0" := .data$s) |>
    dplyr::mutate(s = purrr::pmap(list(.data$s0, .data$fgt), function(x, y) {
      aidx <- nummatch(x$mz, y$accuratemass, delta_x = .Machine$double.eps)
      x$formula[aidx[, 1]] <- y$formula[aidx[, 2]]
      x$exactmass[aidx[, 1]] <- y$exactmass[aidx[, 2]]
      x$error[aidx[, 1]] <- y$error[aidx[, 2]]
      x$label <- x$formula
      return(x)
    }))
  x$MSFres <- vector("list", length = nrow(x))
  for (id in unique(MSFres$id)) {
    x$MSFres[[which(x[[id_column]] == id)]] <- MSFres$fgt[MSFres$id == id]
    x[[in_column]][[which(x[[id_column]] == id)]] <- MSFres$s[[which(MSFres$id == id)[1]]]
  }
  return(x)
}
