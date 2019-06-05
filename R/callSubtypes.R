

#' callOneSubtype
#' Make subtype calls for one sample
#' @export
#' @param mods xgboost model list
#' @param X gene expression matrix, genes in rows, samples in columns
#' @param ci cluster label, and index into mods
#' @return preds of one cluster model.
#' @examples
#' calli <- callOneSubtype(mods, X, 4)
#'
callOneSubtype <- function(mods, X, ci) {

  # Xbin needs to have the same columns as the training matrix...
  print(paste0('calling subtype ', ci))
  mi <- mods[[ci]]
  Xbin <- dataProc(X, mods, ci)
  pred <- predict(mi$bst, Xbin)
  return(pred)
}


#' callSubtypes
#' Make subtype calls for each sample
#' @export
#' @param mods xgboost model list
#' @param X gene expression matrix, genes in rows, samples in columns
#' @return table, column 1 is best call, remaining columns are subtype prediction scores.
#' @examples
#' calls <- callSubtypes(mods, X)
#'
callSubtypes <- function(mods, X) {

  pList <- lapply(1:6, function(mi) callOneSubtype(mods, X, mi))  # was lapply(names(mods), ... )
  pMat  <- do.call('cbind', pList)
  colnames(pMat) <- 1:6 # names(mods)
  bestCall <- apply(pMat, 1, function(pi) colnames(pMat)[which(pi == max(pi)[1])])

  return(cbind(data.frame(BestCall=bestCall), pMat, stringsAsFactors=F))
}


#' callEnsemble
#' Make subtype calls for each sample
#' @export
#' @param ens list, result of fitEnsembleModel
#' @param X gene expression matrix, genes in rows, samples in columns
#' @return table, column 1 is best call, remaining columns are subtype prediction scores.
#' @examples
#' calls <- callEnsemble(mods, X, Y)
#'
callEnsemble <- function(ens, X) {

  eList <- lapply(ens, function(ei) callSubtypes(ei, X))
  eRes <- eRes[,-1] # remove best calls
  eRes <- Reduce('+', eList) / length(eList)
  colnames(eRes) <- 1:6 # names(mods)
  bestCall <- apply(eRes, 1, function(pi) colnames(eRes)[which(pi == max(pi)[1])])

  return(cbind(data.frame(BestCall=bestCall), eRes))
}
