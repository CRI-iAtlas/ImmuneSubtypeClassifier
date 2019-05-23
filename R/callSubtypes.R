

#' callSubtypes
#' Make subtype calls for each sample
#' @export
#' @param mods xgboost model list
#' @param X gene expression matrix, genes in rows, samples in columns
#' @return table, column 1 is best call, remaining columns are subtype prediction scores.
#' @examples
#' mods <- fitSubtypeModel(Xs, Ys, params)
#'
callSubtypes <- function(mods, X) {
  Xbin  <- t(dataProc(X, mods))
  pList <- lapply(mods, function(mi) predict(mi$bst, Xbin))
  pMat  <- do.call('cbind', pList)
  bestCall <- apply(pMat, 1, function(pi) colnames(pMat)[which(pi == max(pi)[1])])
  return(list(Calls=cbind(data.frame(BestCall=bestCall), pMat), Xbin=Xbin))
}
