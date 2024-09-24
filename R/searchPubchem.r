#' Add PubChem metadata to MS-DIAL alignment results
#'
#' Query Pubchem for compound metadata using the webchem package. A cache file
#' is used to store obtained results for subsequent queries. The function (1)
#' translates supplied InChIKeys to Pubchem CIDs and (2) queries PubChem for
#' compound details for these CIDs. Currently, PubChem 'Synonyms', 'Properties',
#' 'Uses' and 'Toxicity' sections are retrieved.
#' @param x data.frame
#' @param in_column name of column containing InChIKey's. Can be "NISTres" in
#'   which case query is performed on a previous NIST search result.
#' @param cache_file full path of results cache file or 'NULL' to disable
#'   caching
#' @return data.frame
#' @export
#' @importFrom rlang .data
#' @importFrom rlang :=
#' @examples
#' \dontrun{
#' ik <- data.frame(inchikey = "MXWJVTOOROXGIU-UHFFFAOYSA-N")
#' searchPubchem(ik)
#' }
searchPubchem <- function(x,
                          in_column = c("inchikey", "NISTres")[1],
                          cache_file = "~/pubchem_cache.rda") {
  ## extract unique inchikeys from NIST results
  stopifnot(in_column %in% colnames(x))
  if (in_column == "NISTres") {
    NISTres <- do.call("rbind", x$NISTres) |>
      tibble::as_tibble() |>
      dplyr::filter(!duplicated(.data$inchikey))
    ikeys <- NISTres$inchikey
  } else {
    ikeys <- x[[in_column]]
  }
  ## check cache for inchikeys
  if (!is.null(cache_file) && file.exists(cache_file)) {
    load(cache_file) # "pc_cache"
    pc_cache0 <- pc_cache
  }
  ikeys_found <- ikeys %in% pc_cache$inchikey
  ikeys_needed <- !ikeys_found
  if (sum(ikeys_found) > 0) {
    message(sprintf(
      "Found chached info for %d of %d unique InChiKey's",
      sum(ikeys_found),
      length(ikeys)
    ))
    pc_cache <- pc_cache[pc_cache$inchikey %in% ikeys[ikeys_found], ]
  }
  if (sum(ikeys_needed) > 0) {
    ikeys <- ikeys[ikeys_needed]
    message(sprintf("Querying Pubchem for %d unique InChiKey's", length(ikeys)))
    cid0 <- webchem::get_cid(
      ikeys,
      from = "inchikey",
      domain = "compound",
      match = "first",
      verbose = TRUE
    )
    cid <- cid0 |>
      dplyr::filter(!is.na(.data$cid)) |>
      dplyr::mutate("cid" := as.numeric(.data$cid)) |>
      dplyr::filter(!duplicated(.data$cid))
    if (nrow(cid) > 0) {
      ## get Pubchem synonym names for compounds
      message(sprintf("Querying details for %d unique CID's", nrow(cid)))
      syn0 <- webchem::pc_synonyms(cid$cid, from = "cid", verbose = T)
      syn <- plyr::ldply(syn0, paste, collapse = "|", .id = "cid") |>
        dplyr::as_tibble() |>
        dplyr::rename("synonyms" := .data$V1) |>
        dplyr::mutate("cid" := as.numeric(levels(.data$cid))[.data$cid])
      ## get Pubchem properties for compounds
      prop0 <- webchem::pc_prop(as.numeric(cid$cid), verbose = TRUE)
      prop <- prop0 |>
        tibble::as_tibble() |>
        dplyr::rename("cid" := .data$CID)
      ## get Pubchem 'Uses' and 'Toxicity' sections
      uses0 <- webchem::pc_sect(cid$cid, section = "Uses", verbose = TRUE)
      uses <- uses0 |>
        dplyr::filter(!is.na(.data$Name)) |>
        dplyr::filter(.data$Result != "cpdat") |>
        dplyr::rename("cid" := .data$CID) |>
        dplyr::group_by(.data$cid) |>
        dplyr::summarize("Result" := paste(.data$Result, collapse = "|")) |>
        dplyr::rename("Uses" := .data$Result) |>
        dplyr::mutate("cid" := as.numeric(.data$cid))
      tox0 <- webchem::pc_sect(cid$cid, section = "Toxicity Summary", verbose = TRUE)
      tox <- tox0 |>
        dplyr::filter(!is.na(.data$Name)) |>
        dplyr::rename("cid" := .data$CID) |>
        dplyr::group_by(.data$cid) |>
        dplyr::summarize("Result" := paste(.data$Result, collapse = "|")) |>
        dplyr::rename("Tox" := .data$Result) |>
        dplyr::mutate("cid" := as.numeric(.data$cid))
      ## combine results
      pcres <- NISTres |>
        dplyr::select("inchikey") |>
        dplyr::left_join(cid, by = c("inchikey" = "query")) |>
        dplyr::filter(!is.na(.data$cid)) |>
        dplyr::left_join(syn, by = "cid") |>
        dplyr::left_join(prop, by = "cid") |>
        dplyr::left_join(uses, by = "cid") |>
        dplyr::left_join(tox, by = "cid")
    }
  }
  if (!is.null(pc_cache) && nrow(pc_cache) > 0) {
    if (exists("pcres") && nrow(pcres) > 0) {
      pcres <- rbind(pcres, pc_cache)
      pcres <- pcres[!duplicated(pcres$inchikey), ]
      if (any(!pcres$inchikey %in% pc_cache0$inchikey)) {
        new_entry <- pcres[!pcres$inchikey %in% pc_cache0$inchikey, ]
        message(sprintf(
          "Updating cache file with %d new entries.",
          nrow(new_entry)
        ))
        pc_cache_out <- rbind(pc_cache0, new_entry)
        pc_cache <- pc_cache_out
        save(pc_cache, file = cache_file)
      } else {
        message("No new entries, leaving cache file unchanged.")
      }
    } else {
      pcres <- pc_cache
    }
  }
  ## attach results to respective object(s)
  out <- if (exists("pcres") && nrow(pcres) > 0) {
    if (in_column == "NISTres") {
      x |>
        dplyr::mutate("NISTres" := purrr::map(NISTres, function(y) {
          dplyr::left_join(y, pcres, by = "inchikey", suffix = c("", "$$")) |>
            dplyr::select(-tidyselect::ends_with("$$"))
        }))
    } else {
      x |>
        dplyr::left_join(pcres, by = "inchikey", suffix = c("", "$$")) |>
        dplyr::select(-tidyselect::ends_with("$$"))
    }
  } else {
    x
  }
  ## met_smiles <- sapply(out$NISTres, function(x) x$CanonicalSMILES[1])
  ## flt <- !is.na(met_smiles)
  ## out$smiles[flt] <- met_smiles[flt]
  return(out)
}
