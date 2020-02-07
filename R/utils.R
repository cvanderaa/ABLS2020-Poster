
####---- Utility functions for parsing SCP data ---####

library(Features)
library(MsCoreUtils)
library(SingleCellExperiment)

##' @title SingleCellExperiment from tabular data
##'
##' @description
##'
##' Convert tabular data from a spreadsheet or a `data.frame` into a
##' `SingleCellExperiment` object.
##'
##' @param table File or object holding the quantitative data. Can be
##'     either a `character(1)` with the path to a text-based
##'     spreadsheet (comma-separated values by default, but see `...`)
##'     or an object that can be coerced to a `data.frame`. It is
##'     advised not to encode characters as factors.
##'
##' @param ecol A `numeric` indicating the indices of the columns to
##'     be used as assay values. Can also be a `character`
##'     indicating the names of the columns. Caution must be taken if
##'     the column names are composed of special characters like `(`
##'     or `-` that will be converted to a `.` by the `read.csv`
##'     function. If `ecol` does not match, the error message will
##'     dislpay the column names as seen by the `read.csv` function.
##'
##' @param fnames An optional `character(1)` or `numeric(1)`
##'     indicating the column to be used as row names.
##'
##' @param ... Further arguments that can be passed on to [read.csv]
##'     except `stringsAsFactors`, which is always `FALSE`.
##'
##' @param name An `character(1)` to name assay. If not set,
##'     `features` is used.
##'
##' @return An instance of class [SingleCellExperiment].
##'
##' @author Laurent Gatto, Christophe Vanderaa
##' 
##' @seealso The code relies on [Features::readSummarizedExperiment].
##'
##' @importFrom utils read.csv
##' @import SingleCellExperiment
##' @importFrom Features readSummarizedExperiment
##'
##' @md
##' @export
##'
readSingleCellExperiment <- function(table, 
                                     ecol, 
                                     fnames,
                                     ...){
  ## Read data as SummarizedExperiment
  sce <- readSummarizedExperiment(table, ecol, fnames, ...)
  sce <- as(sce, "SingleCellExperiment")
  return(sce)
}



#' Split SingleCellExperiment into an ExperimentList
#'
#' The fonction creates an [ExperimentList] containing [SingleCellExperiment] 
#' objects from a [SingleCellExperiment] object. `f` is used to split `x`` along 
#' the rows (`f`` was a feature variable name) or samples/columns (f was a phenotypic variable name). If f is passed as a factor, its length will be matched to nrow(x) or ncol(x) (in that order) to determine if x will be split along the features (rows) or sample (columns). Hence, the length of f must match exactly to either dimension.
#' @param x a single [SingleCellExperiment] object 
#' @param f a factor or a character of length 1. In the latter case, `f`` will 
#'     be matched to the row and column data variable names (in that order). If 
#'     a match is found, the respective variable is extracted, converted to a 
#'     factor if needed
#'
#' @return
#' @export
#'
#' @examples
splitSCE <- function(x, f){
  ## TODO convert function to a "split" method
  ## TODO split this function in two method with signatures
  ## c("SingleCellExperiment", "factor"), c("SingleCellExperiment", "character")
  
  ## Check that f is a factor
  if (is.character(f)) {
    if (length(f) != 1) 
      stop("Character must be of lenght one.")
    if (f %in% colnames(rowData(x))) {
      f <- rowData(x)[, f]
    }
    else if (f %in% colnames(colData(x))) {
      f <- colData(x)[, f]
    }
    else {
      stop(f, " not found in any feature/phenodata variables.")
    }
    if (!is.factor(f)) 
      f <- factor(f)
  }
  ## Check that the factor matches one of the dimensions
  if (!length(f) %in% dim(x)) 
    stop("length(f) not compatible with dim(x).")
  if (length(f) == nrow(x)) { ## Split along rows
    xl <- lapply(split(rownames(x), f = f), function(i) x[i, ])
  } else { ## Split along columns
    xl <- lapply(split(colnames(x), f = f), function(i) x[, i])
  }
  ## Convert list to an ExperimentList
  do.call(ExperimentList, xl)
}


nrowAssays <- function(x){
  out <- numeric(length(x))
  for(i in 1:length(x)){
    out[i] <- nrow(x[[i]])
  }
  return(out)
}

aggregateFeaturesSCE <- function(object, 
                                 i, 
                                 fcol, 
                                 name = "newAssay", 
                                 fun = MsCoreUtils::robustSummary, 
                                 ...){
  if (isEmpty(object))
    return(object)
  if (name %in% names(object))
    stop("There's already an assay named '", name, "'.")
  if (missing(fcol))
    stop("'fcol' is required.")    
  if (missing(i))
    i <- main_assay(object)
  assay_i <- assay(object, i)
  rowdata_i <- rowData(object[[i]])
  if (!fcol %in% names(rowdata_i))
    stop("'fcol' not found in the assay's rowData.")
  groupBy <- rowdata_i[[fcol]]
  
  if (anyNA(assay_i)) {
    msg <- paste("Your data contains missing values.",
                 "Please read the relevant section in the",
                 "aggregateFeatures manual page for details the",
                 "effects of missing values on data aggregation.")
    message(paste(strwrap(msg), collapse = "\n"))
  }
  aggregated_assay <- aggregate_by_vector(assay_i, groupBy, fun, ...)
  aggregated_rowdata <- Features::reduceDataFrame(rowdata_i, rowdata_i[[fcol]],
                                                  simplify = TRUE, drop = TRUE,
                                                  count = TRUE)
  
  se <- SingleCellExperiment(aggregated_assay,
                             rowData = aggregated_rowdata[rownames(aggregated_assay), ])
  hits <- findMatches(rownames(aggregated_assay), groupBy)
  elementMetadata(hits)$names_to <- rowdata_i[[fcol]][hits@to]
  elementMetadata(hits)$names_from <- rownames(assay_i)[hits@to]
  
  assayLinks <- AssayLink(name = name,
                          from = ifelse(is.character(i), i, names(object)[i]),
                          fcol = fcol,
                          hits = hits)
  addAssay(object,
           se,
           name = name,
           assayLinks = assayLinks)
  
}


#' Combine two SingleCellExperiment object
#' 
#' The function takes two [SingleCellExperiment] objects and combines them in 
#' a single object. The assays are combined by 
#'
#' @param x A [SingleCellExperiment] object
#' @param y A [SingleCellExperiment] object
#' @param i A `numeric` or `character` of length 1 indicating which assay
#'     to take from object x Default takes the first assay. 
#' @param j  A `numeric` or `character` of length 1 indicating which assay
#'     to take from object x Default takes the first assay. 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
combineSCE <- function(x, 
                       y, 
                       i = 1, 
                       j = 1,
                       ...) {
  stop("Not finished")
  if (class(x) != class(y)) 
    stop(paste0("objects must be the same class, but are ", 
                class(x), ", ", class(y)))
  assays(x)[i] <- combine(assayData(x)[i], assayData(y)[j])
  phenoData(x) <- combine(phenoData(x), phenoData(y))
  featureData(x) <- combine(featureData(x), featureData(y))
  experimentData(x) <- combine(experimentData(x), experimentData(y))
  protocolData(x) <- combine(protocolData(x), protocolData(y))
  if (validObject(x)) 
    return(x)
}


print.progress <- function(current, total, ...){
  if(length(current)!=1 && length(total) != 1) stop("Arguments \'current\' and \'total\' must have length 1.")
  perc <- format(current/total*100, digits = 0,...)
  cat(paste0("\rProgress: ", perc, " %   "))
  if(current == total) cat ("\n")
  flush.console()
}

