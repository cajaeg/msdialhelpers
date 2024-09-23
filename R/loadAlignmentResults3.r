#' #' Read alignment result file from older MSDIAL versions (<4)
#' #'
#' #' @param txtFile txt file resulting from 'Export Alignment Result'
#' #' @param postCurationFilt filter rows on Post curation result ("Highly correlated to ..." and similar). The default ("") removes most problematic compounds.
#' #' @param dotProductThr filter rows on dotProduct
#' #' @param fillThr filter rows on fill\%
#' #' @param medIntThr filter rows on median intensity (-1 to disable)
#' #'
#' #' @return data.frame
#' #' @export
#' #'
#' #' @examples
#' #' fp <- file.path(system.file(package="msdialhelpers"), "extdata/MSDIAL_Alignment_result.txt")
#' #' res <- loadAlignmentResults3(fp, dotProductThr=-1, fillThr=-1, medIntThr=-1) # for larger projects, 100/0.6/1000 seem reasonable settings
#' #' head(res)
#' loadAlignmentResults3 <- function(txtFile,
#'                                   postCurationFilt = NULL,
#'                                   dotProductThr = -1,
#'                                   fillThr = -1,
#'                                   medIntThr = -1) {
#'   stopifnot(file.exists(txtFile))
#'   useColumns = c(
#'     "Metabolite name",
#'     "Average Mz",
#'     "Average Rt(min)",
#'     "Adduct ion name",
#'     "Formula",
#'     "INCHIKEY",
#'     "SMILES",
#'     "Ontology",
#'     "MS/MS included",
#'     "Post curation result",
#'     "Fill %",
#'     "Dot product",
#'     "Reverse dot product",
#'     "MS1 isotopic spectrum",
#'     "MS/MS spectrum"
#'   )
#'   newColumnNames = c(
#'     "name",
#'     "mz",
#'     "rt",
#'     "adduct",
#'     "baseform",
#'     "inchikey",
#'     "smiles",
#'     "ontology",
#'     "hasMS2",
#'     "postcuration",
#'     "fillperc",
#'     "dotprod",
#'     "revdotprod",
#'     "ms1spec",
#'     "ms2spec"
#'   )
#'   tmp <- read.table(
#'     txtFile,
#'     sep = "\t",
#'     header = T,
#'     skip = 3,
#'     as.is = T,
#'     check.names = F,
#'     comment.char = "",
#'     na.strings = c("NA", "No record")
#'   )
#'   stopifnot("MS/MS spectrum" %in% colnames(tmp)) # something is wrong with input file
#'   samplecols <- seq(grep("MS/MS spectrum", colnames(tmp)) + 1, ncol(tmp))
#'   flt1 <- if (!is.null(postCurationFilt))
#'     tmp$`Post curation result` == postCurationFilt
#'   else
#'     rep(TRUE, nrow(tmp))
#'   flt2 <-
#'     tmp$`Dot product` >= dotProductThr &
#'     tmp$`Fill %` > fillThr &
#'     apply(tmp[, samplecols], 1, median, na.rm = T) > medIntThr # intensity filter
#'   tmp <- tmp[flt1 & flt2, , drop = F]
#'   out_meta <- tmp[, useColumns]
#'   names(out_meta) <- newColumnNames
#'   out_meta <- transform(out_meta, rt_min = rt)
#'   out_meta <- transform(out_meta, rt = round(rt_min * 60, 2))
#'   # out_meta <- transform(out_meta, name=sub(";.*$", "", name)) # clean-up name
#'   out_meta <- transform(out_meta, hasMS2 = as.logical(hasMS2))
#'   out <- tmp[, samplecols]
#'   attr(out, "meta") <- out_meta
#'   return(out)
#' }
