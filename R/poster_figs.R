
library(scpdata)
library(SingleCellExperiment)
library(Features)

peps <- SingleCellExperiment(assays = list(peptide = exprs(specht2019_peptide)), 
                             rowData = fData(specht2019_peptide),
                             colData = pData(specht2019_peptide))

fts <- Features(experiments = list(peptide = peps), 
                colData = colData(peps))
fts <- scp_normalizeFeatures(fts, "peptide", "mean", name = "peptide_norm")
fts <- aggregateFeatures(fts, i = "peptide_norm", fcol = "protein", 
                         name = "protein", fun = colMeans, na.rm = TRUE)
fts <- scp_normalizeFeatures(fts, "protein", "mean", name = "protein_fnorm")
fts <- scp_normalizeSamples(fts, "protein_fnorm", "median", name = "protein_snorm")
scBc <- batchCorrect(scNorm, batch = "raw.file", target = "celltype")



### Function 
scp_normalizeFeatures <- function(x, i, FUN = "mean", na.rm = TRUE, name =NULL){
  if (FUN == "mean") FUN <- function(val) val - mean(val, na.rm = na.rm)
  
  el <- experiments(x)[[i]]
  assays(el)[[1]] <- t(apply(assays(el)[[1]], 1, FUN))
  
  if (is.null(name)) {
    experiments(x)[[i]] <- el
    return(el)
  } else {
    return(addAssay(x, el, name = name))
  }
  
}


scp_normalizeSamples <- function(x, i, FUN = "mean", na.rm = TRUE, name =NULL){
  if (FUN == "mean") FUN <- function(val) val - mean(val, na.rm = na.rm)
  if (FUN == "median") FUN <- function(val) val - median(val, na.rm = na.rm)
  
  el <- experiments(x)[[i]]
  assays(el)[[1]] <- apply(assays(el)[[1]], 2, FUN)
  
  if (is.null(name)) {
    experiments(x)[[i]] <- el
    return(el)
  } else {
    return(addAssay(x, el, name = name))
  }
  
}


batchCorrect <- function(x, batch, target){
  if(is.character(batch)){
    batch <- pData(obj)[, batch]
  } else if (!is.factor(batch)){
    stop("'batch' should be either a column name (character) in pData(obj) or a factor")
  }
  if(is.character(target)){
    target <- model.matrix(~ as.factor(pData(obj)[, target]))
  } else if (!is.matrix(target)){
    stop("'target' should be either a column name (character) in pData(obj) or a design matrix")
  }
  exprs(obj) <- ComBat(dat = exprs(obj), batch = batch, mod = target, par.prior = T)
  return(obj)
}

