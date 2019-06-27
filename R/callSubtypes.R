

#' geneMatch
#' Match the incoming data to what was used in training
#' @export
#' @param X gene expression matrix, genes in rows, samples in columns
#' @return X, but with genes matching the ebpp set, missing genes are NAs
#' @examples
#' Xprime <- geneMatch(X)
#'
geneMatch <- function(X, geneid='pairs') {

  data('ebpp_gene')

  if (geneid == 'symbol') {
    idx <- match(table = rownames(X), x = ebpp_genes_sig$Symbol)
  } else if (geneid == 'entrez') {
    idx <- match(table = rownames(X), x = ebpp_genes_sig$Entrez)
  } else if (geneid == 'ensembl') {
    ensemble <- str_split(rownames(X), pattern = '\\.')
    ensemble <- unlist(lapply(ensemble, function(a) a[1]))
    idx <- match(table = ensemble, x = ebpp_genes_sig$Ensembl)
  } else if (geneid == 'pairs') {
    return(X)
  } else {
    print("For geneids, please use:  symbol, entrez, ensembl")
    return(NA)
  }

  X2 <- X[idx,]  ### Adds NA rows in missing genes
  rownames(X2) <- ebpp_genes_sig$Symbol
  return(X2)
}


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

  pList <- lapply(1:6, function(mi) callOneSubtype(mods, X, mi))
  pMat  <- do.call('cbind', pList)
  colnames(pMat) <- 1:6 # names(mods)
  bestCall <- apply(pMat, 1, function(pi) colnames(pMat)[which(pi == max(pi)[1])])

  return(cbind(data.frame(BestCall=bestCall), pMat, stringsAsFactors=F))
}


#' callEnsemble
#' Make subtype calls for each sample
#' @export
#' @param X gene expression matrix, genes in row.names, samples in column.names
#' @param path the path to the ensemble model, stored as RData, and named 'ens'
#' @param geneids either hgnc for gene symbols or entrez ids. will be matched to the EB++ matrix
#' @return table, column 1 is best call, remaining columns are subtype prediction scores.
#' @examples
#' calls <- callEnsemble(mods, X, Y)
#'
callEnsemble <- function(X, path='data', geneids='symbol') {

  if (path == 'data') {
    data("ensemble_model")
  } else {
    load(path)
  }

  X <- geneMatch(X, geneids)

  eList <- lapply(ens, function(ei) callSubtypes(mods=ei, X=X))
  eRes <- Reduce('+', eList) / length(eList)
  eRes <- eRes[,-1] # remove best calls
  colnames(eRes) <- 1:6 # names(mods)
  bestCall <- apply(eRes, 1, function(pi) colnames(eRes)[which(pi == max(pi)[1])])

  return(cbind(data.frame(BestCall=bestCall), eRes))
}
