

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
callSubtypes <- function(mods, X, Y) {

  pList <- lapply(names(mods), function(mi) callOneSubtype(mods, X, mi))
  pMat  <- do.call('cbind', pList)
  colnames(pMat) <- 1:6 # names(mods)
  bestCall <- apply(pMat, 1, function(pi) colnames(pMat)[which(pi == max(pi)[1])])

  return(cbind(data.frame(Y=Y, BestCall=bestCall), pMat))
}


