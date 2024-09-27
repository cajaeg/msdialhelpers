#' Calculate m/z from chemical formula
#'
#' @param x chemical formula
#' @param adduct adduct ("\[M+H\]+", "\[M\]+" etc.) or NULL for neutral mass
#' @param round digits to round m/z to
#'
#' @return data.frame
#' @export
#'
#' @examples
#' fml2mz("C6H6O", "[M]+") # (radical) cation
#' fml2mz("C6H6O") # neutral mass
fml2mz <- function(x, adduct = NULL, round = 4) {
  isotopes <- NULL
  utils::data("isotopes", package = "enviPat", envir = environment()) # isotopes
  if (is.character(x)) {
    fml <- x
    m <- base::round(enviPat::check_chemform(isotopes, fml)[, "monoisotopic_mass"], round)
  } else if (is.numeric(x)) {
    fml <- NA_character_
    m <- x
  } else {
    return(NA)
  }
  if (is.null(adduct)) {
    r <- data.frame(
      name = "M",
      nmol = 1,
      charge = 0,
      massdiff = 0
    )
  } else {
    r <- ion2rule(adduct)
  }
  adductmz <- (m * r$nmol + r$massdiff) / pmax(1, abs(r$charge))
  n_adducts <- nrow(r)
  data.frame(
    formula = rep(fml, n_adducts),
    neutral_mass = rep(m, n_adducts),
    adduct = r$name,
    nmol = r$nmol,
    charge = r$charge,
    mz = base::round(adductmz, round)
  )
}
