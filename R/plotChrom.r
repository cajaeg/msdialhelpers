#' Plot extracted ion chromatograms (EIC, BPC, TIC)
#'
#' @param chrom chromatogram as obtained by \code{\link{getChrom}()}
#' @param label_peaks label major peaks with retention times
#' @param fwhm full-width half maximum for peak detection (passed to
#'   \code{\link[xcms]{peaksWithMatchedFilter}()})
#' @param label.k number of chromatogram section to analyze: in each section the
#'   major peak will be labelled
#' @param legend display a legend of m/z values
#' @param max.legend how many m/z's to show in legend
#' @param set.mfrow set multi-figure layout according to number of chromatograms
#' @param expand.ylim include some extra space for peak labels
#' @param ... passed to underlying function \code{graphics::matplot()}
#'
#' @return (invisibly) original 'chrom' object with attribute 'peaks' attached if
#'   peak detection was performed
#' @export
#'
#' @examples
#' # see getChrom()
plotChrom <-
  function(chrom,
           label_peaks = TRUE,
           fwhm = 5,
           label.k = 5,
           legend = NULL,
           max.legend = 15,
           set.mfrow = TRUE,
           expand.ylim = 1.05,
           ...) {
    if (inherits(chrom, "list")) {
      if (set.mfrow) {
        opar <- graphics::par(mfrow = grDevices::n2mfrow(length(chrom)))
        on.exit(graphics::par(opar))
      }
      args <- list(...)
      argNames <- names(args)
      flt <- sapply(chrom, inherits, "matrix")
      gylim <- c(0, max(0, unlist(lapply(chrom[flt], as.vector)), na.rm = TRUE) * expand.ylim)
      args <- checkArgs(args, "ylim", gylim)
      args <- checkArgs(args, "xlab", "")
      args <- checkArgs(args, "ylab", "")
      args[["label_peaks"]] <- label_peaks
      args[["fwhm"]] <- fwhm
      args[["label.k"]] <- label.k
      args[["legend"]] <- legend
      args[["max.legend"]] <- max.legend
      args[["expand.ylim"]] <- expand.ylim
      for (i in 1:length(chrom)) {
        if (flt[i]) {
          args[["chrom"]] <- chrom[[i]]
          do.call(plotChrom, args)
          if (!is.null(attr(chrom[[i]], "from_file")))
            graphics::mtext(
              sub("(.+)(\\..+)$", "\\1", basename(attr(
                chrom[[i]], "from_file"
              ))),
              side = 3,
              adj = 0.01,
              cex = 0.66
            )
        } else {
          emptyplot("no data")
        }
      }
      return(invisible(NULL))
    }
    args <- list(...)
    argNames <- names(args)
    args <- checkArgs(args, "type", "l")
    args <- checkArgs(args, "lty", 1:5)
    args <- checkArgs(args, "col", 1:6)
    if (!is.null(.tmp <- attr(chrom, "rt"))) {
      chrom_rt <- .tmp / 60
      xlab <- "RT (min)"
    } else if ("x" %in% argNames) {
      chrom_rt <- args[["x"]] / 60
      xlab <- "RT (min)"
    } else {
      warning("no RT provided, using index")
      chrom_rt <- seq_len(nrow(chrom))
      fwhm <- 0.1
      xlab <- "Index"
    }
    args[["x"]] <- chrom_rt
    args <- checkArgs(args, "xlab", xlab)
    args[["y"]] <- chrom
    args <- checkArgs(args, "ylab", "Intensity (a.u.)")
    args <- checkArgs(args, "ylim", c(0, max(1, chrom, na.rm = TRUE) * expand.ylim))
    if (is.null(args[["ylim"]]))
      args[["ylim"]] <- c(0, max(1, chrom, na.rm = TRUE) * expand.ylim)
    if (nrow(chrom) > 0) {
      do.call(graphics::matplot, args = args)
      if (!any(is.finite(chrom)))
        graphics::text(stats::median(graphics::par("usr")[1:2]), 
                       stats::median(graphics::par("usr")[3:4]), 
                       labels = "no data")
      if (label_peaks) {
        rt <- chrom_rt * 60
        n_ions <- ncol(chrom)
        int <- if (n_ions == 1)
          as.vector(chrom[, 1])
        else
          apply(chrom, 1, function(x)
            max(c(0, x), na.rm = TRUE))
        pks <- data.frame(suppressWarnings(
          xcms::peaksWithMatchedFilter(
            int = int,
            rt = rt,
            fwhm = fwhm,
            max = 50,
            snthresh = 2
          )
        ))
        if (nrow(pks) > 0) {
          pks$rt <- pks$rt / 60
          pks <- pks[order(pks$rt), ]
          pks$grp <- hcgroup(pks$rt, k = min(nrow(pks), label.k))
          flt <- unlist(tapply(-pks$maxo, pks$grp, function(x)
            rank(x, ties.method = "first") == 1, simplify = FALSE))
          pks <- pks[flt, ]
          graphics::text(
            pks$rt,
            pks$maxo,
            labels = round(pks$rt, 2),
            pos = 3,
            offset = 0.1
          )
        }
      }
      show.legend <- (is.logical(legend) && legend) ||
        is.character(legend) || is.numeric(legend) ||
        (is.null(legend) && ncol(chrom) > 1)
      show.legend <- show.legend && !is.null(attr(chrom, "mz"))
      if (show.legend) {
        legtext <- sprintf("%.3f", round(attr(chrom, "mz"), 3))
        if (is.character(legend) || is.numeric(legend)) {
          legtext <- paste(legtext, legend)
        }
        if (!is.null(legtext)) {
          maxleg <- max.legend
          legtext <- legtext[1:min(maxleg, length(legtext))]
          graphics::legend(
            "topright",
            bty = "n",
            legend = legtext,
            lty = args[["lty"]],
            col = args[["col"]]
          )
        }
      }
    } else {
      emptyplot("no data")
    }
    rv <- chrom
    attr(rv, "rt") <- chrom_rt * 60
    attr(rv, "peaks") <- if (label_peaks &&
                             exists("pks", envir = environment()))
      pks
    else
      NULL
    invisible(rv)
  }
