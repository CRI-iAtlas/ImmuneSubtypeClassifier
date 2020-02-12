

#' geneMatch
#' Match the incoming data to what was used in training
#' @export
#' @param X gene expression matrix, genes in rows, samples in columns
#' @return X, but with genes matching the ebpp set, missing genes are NAs
#' @examples
#' Xprime <- geneMatch(X)
#'
geneMatch <- function(X, geneid='pairs') {  ## add datasource param  (RNAseq or io360)

  data(ebpp_gene)

  if (geneid == 'symbol') {
    idx <- match(table = rownames(X), x = ebpp_genes_sig$Symbol)  ### this is just for the EBPP genes ###

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

  # idx will be 485 elements long... non matched ebpp_sig_genes
  # will show as NAs in the list.
  
  # SO... we calculate sum of NAs over size of ebpp_genes_sig

  matchError <- sum(is.na(idx)) / nrow(ebpp_genes_sig)

  # NAs in idx will enter NA rows in X2 
  
  X2 <- X[idx,]  ### Adds NA rows in missing genes
  rownames(X2) <- ebpp_genes_sig$Symbol

  return(list(Subset=X2, matchError=matchError))
}


#' geneMatchErrorReport
#' Check whether the incoming data matches the 485 model gene IDs
#' @export
#' @param X gene expression matrix, genes in rows, samples in columns
#' @return list with percent missing genes and a vector of missing genes
#' @examples
#' missingGenes <- geneMatchErrorReport(X)
#'
geneMatchErrorReport <- function(X, geneid='pairs') {
  data(ebpp_gene)
  
  if (geneid == 'symbol') {
    idx <- match(table = rownames(X), x = ebpp_genes_sig$Symbol)  ### this is just for the EBPP genes ###
    
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
  
  # idx will be 485 elements long... non matched ebpp_sig_genes
  # will show as NAs in the list.
  
  # SO... we calculate sum of NAs over size of ebpp_genes_sig
  
  matchError <- sum(is.na(idx)) / nrow(ebpp_genes_sig)
  
  # NAs in idx will enter NA rows in X2 
  
  g <- ebpp_genes_sig[is.na(idx),]  ### Adds NA rows in missing genes

  return(list(matchError=matchError, missingGenes=g))
}

reportError <- function(err) {
  print("**************************************")
  print("    Gene Match Error Report           ")
  print("                                      ")
  print(paste0("  percent missing genes: ",err*100,"           "))
  print("                                      ")
  print("see ?geneMatchErrorReport for details ")
  print("                                      ")
  print("**************************************")
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

  return(data.frame(SampleID=colnames(X), BestCall=bestCall, pMat, stringsAsFactors=F))
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
callEnsemble <- function(X, path='data', geneids='symbol') {  ## add new parameter, RNA-seq or io360

  ## if datasource == 'RNA-seq'
  
     data('subtype_caller_model')  ## This is only for a EBPP classifier ##

  ## else if the datasource == 'io360 ###

     ## data('subtype_caller_io360_model')
  
  if (path == 'data') {  ## and datasource == 'RNAseq'
    data("ensemble_model")
    
  ## else if path == 'data' and datasource == 'io360'
    ##data(io360_model)
    
  } else {
    load(path)
  }
  
  
  res0 <- geneMatch(X, geneids) ## datasource param here.
  
  X <- res0$Subset
  matchError <- res0$matchError
  reportError(matchError)
  
  eList <- lapply(ens, function(ei) callSubtypes(mods=ei, X=X))
  ePart <- lapply(eList, function(a) a[,3:8])
  eStack <- array( unlist(ePart) , c(ncol(X), 6, length(ens)) )
  eMeds  <- apply( eStack , 1:2 , median )
  eMeds <- as.data.frame(eMeds)
  colnames(eMeds) <- 1:6 # names(mods)

  #bestCall <- apply(eMeds, 1, function(pi) colnames(eMeds)[which(pi == max(pi)[1])])
  predCall <- predict(scaller, as.matrix(eMeds)) + 1

  sampleIDs <- eList[[1]][,1]

  res0 <- data.frame(SampleIDs=sampleIDs, BestCall=predCall, eMeds)
  colnames(res0)[3:8] <- 1:6
  return(res0)
}



#' parCallEnsemble
#' Parallel version: Make subtype calls for each sample
#' @export
#' @param X gene expression matrix, genes in row.names, samples in column.names
#' @param path the path to the ensemble model, stored as RData, and named 'ens'
#' @param geneids either hgnc for gene symbols or entrez ids. will be matched to the EB++ matrix
#' @param numCores number of cores to use with parallel lib
#' @return table, column 1 is best call, remaining columns are subtype prediction scores.
#' @examples
#' calls <- callEnsemble(mods, X, Y)
#'
parCallEnsemble <- function(X, path='data', geneids='symbol', numCores=2) {

  data('subtype_caller_model')

  if (path == 'data') {
    data("ensemble_model")
  } else {
    load(path)
  }

  X <- geneMatch(X, geneids)
  matchError <- res0$matchError
  reportError(matchError)

  cl <- makeForkCluster(numCores)

  #eList <- lapply(ens, function(ei) callSubtypes(mods=ei, X=X))
  eList <- parLapply(cl=cl, X=1:length(ens), fun=function(ei) callSubtypes(mods=ens[[ei]], X=X))

  stopCluster(cl)

  ePart <- lapply(eList, function(a) a[,3:8])
  eStack <- array( unlist(ePart) , c(ncol(X), 6, length(ens)) )
  eMeds  <- apply( eStack , 1:2 , median )
  eMeds <- as.data.frame(eMeds)
  colnames(eMeds) <- 1:6 # names(mods)

  #bestCall <- apply(eMeds, 1, function(pi) colnames(eMeds)[which(pi == max(pi)[1])])
  predCall <- predict(scaller, as.matrix(eMeds)) + 1

  sampleIDs <- eList[[1]][,1]

  res0 <- data.frame(SampleIDs=sampleIDs, BestCall=predCall, eMeds)
  colnames(res0)[3:8] <- 1:6
  return(res0)
}
