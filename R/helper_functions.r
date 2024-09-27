#' check if arg is present in arg list, if not, set it to 'val' 
#'
#' (internal function)
#' @param args arg list
#' @param name arg
#' @param val default value
#'
#' @return list
#' @noRd
#'
checkArgs <- function(args, name, val) {
  if(!name %in% names(args))
    args[[name]] <- val
  return(args)
}


#' Check if a package is installed
#'
#' (internal function)
#' @param pkg package to check (character(1))
#' @param libPaths library paths to check, defaults to \code{.libPaths()}
#'
#' @return logical
#' @noRd
#' 
isInstalled <- function(pkg, libPaths = .libPaths()) {
  libPaths <- match.arg(libPaths)
  any(grepl(pkg, basename(
    list.dirs(libPaths, recursive = FALSE)
  )))
}


#' m/z or RT clustering based on 'hclust'
#'
#' @param x numeric vector
#' @param h tree height at which to cut, see \code{\link[stats]{cutree}()}
#' @param k desired number of groups, see \code{\link[stats]{cutree}()}
#' @param method cluster method, see \code{\link[stats]{hclust}()}
#' @param useFastCluster use \code{hclust.vector()} from \code{fastCluster}
#'   package instead of the standard \code{stats::hclust()}. Default is to use
#'   \code{fastCluster} if it is installed.
#'
#' @return numeric vector representing group structure
#' @export
#'
#' @examples
#' set.seed(123)
#' x <- abs(rnorm(100, 30, 15))
#' g <- hcgroup(x, h = 5, useFastCluster = FALSE)
#' g1 <- hcgroup(x, h = 5, useFastCluster = TRUE)
#' identical(g, g1)
#' plot(data.frame(x = x, y = 1), type = "h", col = g, ylim = c(0, 1))
hcgroup <- function(x,
                    h = NULL,
                    k = NULL,
                    method = "centroid",
                    useFastCluster = NULL) {
  smallNumber <- .Machine$double.xmin
  if (is.null(useFastCluster) ||
      (is.logical(useFastCluster) && useFastCluster))
    useFastCluster <- isInstalled("fastcluster")
  if (length(x) >= 2) {
    cl <- if (useFastCluster) {
      fastcluster::hclust.vector(c(smallNumber, x), method = method)
    } else {
      hc <- stats::hclust(stats::dist(c(smallNumber, x)) ^ 2, method = method)
      hc$height <- sqrt(hc$height)
      hc
    }
    stats::cutree(cl, h = h, k = k)[-1] - 1
  }
  else
    1
}


#' replace zero values with NA
#'
#' (internal function)
#' @param x numeric vector
#'
#' @return numeric vector
#' @noRd
#'
zero2NA <- function(x) {
  x[x==0] <- NA
  return(x)
}


#' replace NA with zero
#'
#' (internal function)
#' @param x numeric vector
#'
#' @return numeric vector
#' @noRd
#'
NA2zero <- function(x) {
  x[is.na(x)] <- 0
  return(x)
}


#' replace Inf with zero
#'
#' (internal function)
#' @param x numeric vector
#'
#' @return numeric vector
#' @noRd
#'
Inf2zero <- function(x) {
  x[!is.finite(x)] <- 0
  return(x)
}


#' simple replacement for 'plyr::ldply()'
#'
#' (internal function)
#' @param l list
#'
#' @return data.frame
#' @noRd
#'
do.rbind <- function(l) {
  data.frame(do.call(rbind, l), id = rep(1:length(l), sapply(l, nrow)))
}


#' Skip plot in multi-figure layout by plotting "nothing"
#'
#' @param message optional text message to show, e.g. "no data"
#'
#' @return NULL
#' @export
#'
emptyplot <- function(message = NULL) {
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 2), ylim = c(0, 2))
  if(!is.null(message))
    graphics::text(1, 1, labels = message)
  invisible(NULL)
}


#' Calculate isotope pattern using 'enviPat'
#' 
#' (internal function)
#' @param fml chemical formula
#' @param resolution mass spectral resolution (m/delta_m)
#' @param ... passed to 'enviPat::isowrap()'
#'
#' @return list
#' @noRd
#'
#' @examples
#' \dontrun{
#' fml2iso("C6H6O")
#' }
fml2iso <- function(fml, resolution = 20000, ...) {
  isotopes <- NULL
  utils::data("isotopes", package = "enviPat", envir = environment()) # "isotopes"
  fixMS <- function(ms) {
    colnames(ms) <- c("mz", "i")
    ms[, 1] <- round(ms[, 1], 5)
    ms[, 2] <- round(ms[, 2], 1)
    ms
  }
  out <- suppressMessages(
    enviPat::isowrap(
      isotopes,
      enviPat::check_chemform(isotopes, fml),
      resmass = FALSE,
      resolution = resolution,
      ...
    )
  )
  return(lapply(out, fixMS))
}


#' ion2rule
#'
#' (internal function used by fml2mz())
#' @param ions character vector, e.g. "\[M+H\]+"
#'
#' @return data.frame
#' @noRd
#'
#' @examples
#' \dontrun{
#' ion2rule("[M+H]+")
#' }
ion2rule <- function(ions = "[M+H]+") {
  checkSymbol <- function(ion) {
    regexpr("\\[[0-9]{0,2}M.*\\][0-9]{0,2}[\\+\\-]{1,2}", ion) != -1
  }
  shortCuts <- cbind(
    c("M+H", "M+Na", "M+K", "M+NH4", "M+", "M", "M-H", "M+Cl-", "M-"),
    c(
      "[M+H]+",
      "[M+Na]+",
      "[M+K]+",
      "[M+NH4]+",
      "[M]+",
      "[M]+",
      "[M-H]-",
      "[M+Cl]-",
      "[M]-"
    )
  )
  emass <- NULL
  chemical_elements <- NULL
  utils::data(emass, envir = environment(), package = "msdialhelpers")
  utils::data(chemical_elements, envir = environment(), package = "msdialhelpers")
  out <- lapply(ions, function(ion) {
    if (ion %in% shortCuts[, 1])
      ion <- shortCuts[, 2][which(shortCuts[, 1] == ion)]
    if (!checkSymbol(ion))
      stop("invalid ion")
    nmol <- sub(".*[^0-9M]([0-9]?M).*", "\\1", ion)
    nmol <- sub("M", "", nmol)
    nmol <- as.numeric(ifelse(nmol == "", 1, nmol))
    ch <- sub(".*[^0-9]([0-9]{0,2}[\\+\\-])$", "\\1", ion)
    sgn <- sub("[^\\+\\-]", "", ch)
    sgn <- ifelse(sgn == "+", 1, -1)
    ch <- sub("[\\+\\-]", "", ch)
    ch[ch == ""] <- "1"
    ch <- as.numeric(ch)
    ch <- ch * sgn
    x <- ion
    x <- sub("^.*\\[", "", x)
    x <- sub("\\].*", "", x)
    x <- sub("[0-9]?M", "", x)
    starts <- gregexpr("[\\+\\-]", x)[[1]]
    ends <- c(starts[-1] - 1, nchar(x))
    n <- length(starts)
    spl <- lapply(1:n, function(i)
      substr(x, starts[i], ends[i]))
    massdiff <- lapply(spl, function(y) {
      sgn <- sub("^([\\+\\-]).*", "\\1", y)
      sgn <- ifelse(sgn == "+", 1, -1)
      el <- sub("^[\\+\\-]", "", y)
      if (regexpr("^[0-9]+[A-Za-z]+", el) != -1)
        el <- gsub("([0-9]+)([A-Za-z]+)", "\\2\\1", el)
      el <- fml2tbl(el)
      masses <- sapply(colnames(el), function(a) {
        chemical_elements[, 2][which(chemical_elements[, 1] == a)[1]]
      })
      return(sum(masses * el[1, ]) * sgn)
    })
    massdiff <- sum(unlist(massdiff), na.rm = TRUE) + ch *
      -emass
    return(
      data.frame(
        name = ion,
        nmol = nmol,
        charge = ch,
        massdiff = massdiff,
        stringsAsFactors = FALSE
      )
    )
  })
  return(do.call("rbind", out))
}


#' Tabulate elements in chemical formula
#'
#' (internal function used by fml2mz())
#' @param fml chemical formula
#' @param elements chemical elements to consider, or NULL for all elements
#'   occurring in 'fml'
#'
#' @return data.frame
#' @noRd
#' 
#' @examples
#' \dontrun{
#' fml2tbl("C6H6O")
#' fml2tbl("C6H6O", elements = c("C", "H", "O", "S"))
#' }
#' 
fml2tbl <- function (fml, elements = NULL) {
  getElements <- function(form) {
    starts <- gregexpr("[A-Z]", form)[[1]]
    stops <- c(starts[-1] - 1, nchar(form))
    elements <- sapply(1:length(starts), function(i)
      substr(form, starts[i], stops[i]))
    cnts <- sub("[A-Za-z]+", "", elements)
    cnts[cnts == ""] <- 1
    cnts <- as.numeric(cnts)
    elements <- sub("[0-9]+", "", elements)
    return(sort(rep(elements, cnts)))
  }
  checkFormula <- function(form) {
    el <- getElements(form)
    paste0(names(table(el)), table(el), collapse = "")
  }
  if (!is.character(fml))
    return(NULL)
  fml <- gsub("[^A-Za-z0-9]", "", fml)
  fml <- checkFormula(fml)
  if (is.null(elements)) {
    elements <- unique(unlist(lapply(fml, getElements)))
  }
  out <- matrix(ncol = length(elements), nrow = length(fml))
  for (i in seq(along = elements)) {
    el <- elements[i]
    rgx <- sprintf(".*%s([0-9]+).*", el)
    cnt <- sub(rgx, "\\1", fml)
    cnt[regexpr("[[:alpha:]]", cnt) == 1] <- "0"
    cnt[cnt == ""] <- "1"
    cnt <- as.numeric(cnt)
    out[, i] <- cnt
  }
  colnames(out) <- elements
  return(out[, order(colnames(out)), drop = FALSE])
}
