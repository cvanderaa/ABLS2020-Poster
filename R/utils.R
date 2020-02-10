
####---- Utility functions for parsing SCP data ---####

library(Features)
library(MsCoreUtils)
library(SingleCellExperiment)
library(sva)

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

scp_aggregateFeatures <- function(object, 
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
                             rowData = aggregated_rowdata[rownames(aggregated_assay), ],
                             colData = colData(object[[i]]))
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
                       ...) {
  stop("Not finished")
  ## TODO remove this check if function become method
  if (!all(c(class(x), class(y) == "SingleCellExperiment"))) 
    stop("'x' and 'y' must be objects of class 'SingleCellExperiment'")
  
  
  
  t1 <- assays(x)[[1]]
  t2 <- assays(y)[[1]]
  tcomb <- combine(t1, t2)
  
  dim(t1)
  dim(t2)
  dim(tcomb)
  
  assays(x)[i] <- combine(assays(x)[1], assays(y)[1])
  phenoData(x) <- combine(phenoData(x), phenoData(y))
  featureData(x) <- combine(featureData(x), featureData(y))
  experimentData(x) <- combine(experimentData(x), experimentData(y))
  protocolData(x) <- combine(protocolData(x), protocolData(y))
  if (validObject(x)) 
    return(x)
}


combineFeatures <- function(x, 
                            assays, 
                            name = "combined"){
  ## Check that all objects to combine are of same class
  cl <- unique(vapply(experiments(x[, , assays]), class, character(1L)))
  if (length(cl) != 1) 
    stop("Only experiments of same class can be combined. ", 
         paste0("'", cl, "'" , collapse = ", "), " classes were found.")
  ## Create the combined assay and add it to the Features object
  addAssay(object = x, 
           x = do.call(combineSCE, experiments(x[, , assays])), 
           assayLinks = AssayLinks(names = name), ## TODO think about creating assay links 
           name = name)
}



####---- Potential new Features functions ----#### 


getAssayLink <- function(x, ## SummarizedExperiment object to link from
                         y, ## SummarizedExperiment object to link to
                         from, ## name of object to link from
                         name, ## name of object to link to
                         fcol) {  ## the name of the shared feature variable in the respective rowData
  ## Get mathc between x and y rows
  fcolx <- rowData(x)[, fcol]
  fcoly <- rowData(y)[, fcol]
  hits <- findMatches(fcolx, fcoly)
  ## Add the row values used for matching
  elementMetadata(hits)$names_from <- rownames(x)[hits@to]
  elementMetadata(hits)$names_to <- fcolx[hits@to]
  AssayLink(name = name,
            from = from,
            fcol = fcol,
            hits = hits)
}

####---- SCP specific functions ----####



scp_filterSCR <- function(x, 
                          samples, 
                          carrier, 
                          i = 1,
                          thresh = 0.1){
  # Check arguments
  if(!inherits(x, "SingleCellExperiment")) 
    stop("'x' must be an 'SingleCellExperiment' object")
  if(length(carrier) != 1) stop("'carrier' must have length 1.")
  # Compute ratios
  ratios <- apply(assays(x)[[i]][, samples, drop = FALSE], 2, 
                  function(val) val/assays(x)[[i]][, carrier])
  # Filter data
  x[which(rowMeans(ratios, na.rm = TRUE) <= thresh), ]
}




scp_normalizeFeatures <- function(x, 
                                  i,
                                  fcol,
                                  FUN = "mean", 
                                  na.rm = TRUE, 
                                  name = NULL){
  if (FUN == "mean") 
    FUN <- function(val) val - mean(val, na.rm = na.rm)
  if(is.numeric(i))
    i <- names(x)[i]
  if (is.null(name))
    name <- paste0(i, "_fnorm")
  el <- experiments(x)[[i]]
  assays(el)[[1]] <- t(apply(assays(el)[[1]], 1, FUN))
  return(addAssay(x, el, name = name, 
                  assayLinks = getAssayLink(x = experiments(x)[[i]],
                                            y = el, 
                                            from = i,
                                            name = name, 
                                            fcol = fcol)))
}


scp_normalizeSamples <- function(x, 
                                 i, 
                                 fcol,
                                 FUN = "mean", 
                                 na.rm = TRUE, 
                                 name = NULL){
  if(is.numeric(i))
    i <- names(x)[i]
  if (is.null(name))
    name <- paste0(i, "_fnorm")
  if (FUN == "mean") FUN <- function(val) val - mean(val, na.rm = na.rm)
  if (FUN == "median") FUN <- function(val) val - median(val, na.rm = na.rm)
  
  el <- experiments(x)[[i]]
  assays(el)[[1]] <- apply(assays(el)[[1]], 2, FUN)
  return(addAssay(x, el, name = name,
                  assayLinks = getAssayLink(x = experiments(x)[[i]],
                                            y = el, 
                                            from = i,
                                            name = name, 
                                            fcol = fcol)))
}


scp_impute <- function(x, 
                       i, 
                       name, 
                       fcol, 
                       method, ...){
  if(is.numeric(i))
    i <- names(x)[i]
  if (is.null(name))
    name <- paste0(i, "_imputed")
  el <- experiments(x)[[i]]
  assays(el)[[1]] <- method(assays(el)[[1]], ...)
  return(addAssay(x, el, name = name,
                  assayLinks = getAssayLink(x = experiments(x)[[i]],
                                            y = el, 
                                            from = i,
                                            name = name, 
                                            fcol = fcol)))
}
KNNimputation <- function(dat, k = 3){
  # Create a copy of the data, NA values to be filled in later
  dat.imp<-dat
  
  # Calculate similarity metrics for all column pairs (default is Euclidean distance)
  dist.mat<-as.matrix( dist(t(dat)) )
  #dist.mat<-as.matrix(as.dist( dist.cosine(t(dat)) ))
  
  # Column names of the similarity matrix, same as data matrix
  cnames<-colnames(dist.mat)
  
  # For each column in the data... 
  for(X in cnames){
    
    # Find the distances of all other columns to that column 
    distances<-dist.mat[, X]
    
    # Reorder the distances, smallest to largest (this will reorder the column names as well)
    distances.ordered<-distances[order(distances, decreasing = F)]
    
    # Reorder the data matrix columns, smallest distance to largest from the column of interest
    # Obviously, first column will be the column of interest, column X
    dat.reordered<-dat[ , names(distances.ordered ) ]
    
    # Take the values in the column of interest
    vec<-dat[, X]
    
    # Which entries are missing and need to be imputed...
    na.index<-which( is.na(vec) )
    
    # For each of the missing entries (rows) in column X...
    for(i in na.index){
      
      # Find the most similar columns that have a non-NA value in this row
      closest.columns<-names( which( !is.na(dat.reordered[i, ])  ) )
      
      # If there are more than k such columns, take the first k most similar
      if( length(closest.columns)>k ){
        # Replace NA in column X with the mean the same row in k of the most similar columns
        vec[i]<-mean( dat[ i, closest.columns[1:k] ] )
      }
      
      # If there are less that or equal to k columns, take all the columns
      if( length(closest.columns)<=k ){
        # Replace NA in column X with the mean the same row in all of the most similar columns
        vec[i]<-mean( dat[ i, closest.columns ])
      }
    }
    # Populate a the matrix with the new, imputed values
    dat.imp[,X]<-vec
  }
  return(dat.imp)
}

scp_batchCorrect <- function(x, 
                             i, 
                             fcol,
                             batch, 
                             target, 
                             name = NULL){
  if(is.numeric(i))
    i <- names(x)[i]
  if (is.null(name))
    name <- paste0(i, "_batchCorrected")
  el <- experiments(x)[[i]]
  if(length(batch) == 1)
    batch <- colData(el)[, batch]
  batch <- as.factor(batch)
  
  if(length(target) == 1){
    target <- model.matrix(~ as.factor(colData(el)[, target]))
  }
  assays(el)[[1]] <- ComBat(dat = assays(el)[[1]], 
                            batch = batch, 
                            mod = target, 
                            par.prior = T)
  return(addAssay(x, el, name = name,
                  assayLinks = getAssayLink(x = experiments(x)[[i]],
                                            y = el, 
                                            from = i,
                                            name = name, 
                                            fcol = fcol)))
}
