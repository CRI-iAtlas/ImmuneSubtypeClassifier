
#' testBinFun
#' Get difference in sums
#'
#' @param G Binary values around a previously selected pivot point
#' @param Ybin Binary phenotype vector.
#' @return test result, numeric value, the rank sum test.
#' @examples
#' res1 <- testBinFun(G, Ybin)
#' @export
testBinFun <- function(G, Ybin) {
  testres <- (sum(G[Ybin == 1])) - (sum(G[Ybin == 0]))
  return(testres)
}

#' testFun
#' Get difference in mean rank sums for a single gene
#'
#' @param G Gene expression profile for single sample
#' @param Ybin Binary phenotype vector.
#' @return test result, numeric value, the rank sum test.
#' @examples
#' res1 <- testFun(G, Ybin)
#' @export
testFun <- function(G, Ybin) {
  rankg <- rank(G)
  testres <- (sum(rankg[Ybin == 0]) / sum(Ybin == 0)) - (sum(rankg[Ybin == 1]) / sum(Ybin == 1))
  return(testres)
}

#' featureSelection
#' Subset the genes, given a matrix
#' @export
#' @param Xmat Matrix of gene expression data.
#' @param Ybin Binary phenotype vector.
#' @return Xsub, subset of Xmat by genes
#' @examples
#' Xsub <- featureSelection(Xmat, Ybin, 0.1)
#'
featureSelection <- function(Xmat, Ybin, testRes, ptail=0.05) {
  idx <- which( (testRes < quantile(testRes, ptail, na.rm = T)) |
                  (testRes > quantile(testRes, 1.0-ptail, na.rm = T)) )
  Xsub <- Xmat[idx,]
  Xsub[is.na(Xsub)] <- 0
  Xgenes <- rownames(Xmat)[idx]
  return(list(Xsub=Xsub, Genes=Xgenes))
}

#' featureSelection
#' Subset the genes, given a matrix
#' @param x One gene expression sample profile.
#' @param breakVec a vector of probabilities
#' @return xbin, binned expression profile
#' @examples
#' Xsub <- featureSelection(Xmat, Ybin, 0.1)
#'
breakBin <- function(x, breakVec){
  brks <- quantile(as.numeric(x), probs=breakVec, na.rm = T)
  xbin <- .bincode(x = x, breaks = brks, include.lowest = T)
  xbin <- as.numeric(xbin)
  xbin
}



binaryGene <- function(values) {
  # gidx
  pivotvalue <- values[1]
  values <- values[-1]
  res0 <- sapply(values, function(b) as.numeric(b >= pivotvalue))
  res0[is.na(res0)] <- rbinom(n = sum(is.na(res0)), prob = 0.5, 1)  ## replace NAs with random values
  return(res0)
}


#' createPairsFeatures
#' given a matrix and gene pairs, create binary features
#' @param X One gene expression matrix
#' @param genePairs a vector of gene pairs from the ensemble list
#' @return xbin, binned expression profile
#' @examples
#' Xsub <- createPairsFeatures(Xmat, genepairs)
#'
createPairsFeatures <- function(X, genes) {

  genePairs <- strsplit(genes, ':')
  gidx <- rownames(X)

  # first convert the gene pairs to a named list
  # where each entry of the list is the genes for a given pivot-gene
  pairList <- list()
  for(gi in genePairs) {
    pairList[[gi[1]]] <- c(pairList[[gi[1]]], gi[2])
  }

  resList <- list() # then for each pivot gene
  for (gi in names(pairList)) {
    # assuming it's in the data ... really should be!
    if (gi %in% rownames(X)) {
      gs <- unique(c(gi,pairList[[gi]]))  ## can end up with the pivot gene in the genes... guarentee
      if (gs[1] != gi) {
        print('ERROR: first gene does not match pivot gene')
        return(NA)
      }
      idx <- match(table=rownames(X), x=gs) ## get index to genes for this pivot
      Xsub <- X[idx,]                                         ## subset the matrix, NAs for missing genes, pivot gene on top
      rownames(Xsub) <- gs                                    ## give gene IDs
      resList[[gi]] <- apply(Xsub, 2, function(a) binaryGene(a))  ## create binary values
    } else { # else we need to include some dummy rows
      randMat <- matrix(data=rbinom(n = length(pairList[[gi]]) * ncol(X), prob = 0.5, 1), ncol=ncol(X))
      colnames(randMat) <- colnames(X)
      resList[[gi]] <- randMat
    }
  }
  newMat <- do.call('rbind', resList)
  rownames(newMat) <- genes
  newMat <- t(newMat)
  return(newMat)
}


#' trainDataProc
#' Data preprocessing
#' @export
#' @param Xmat Matrix of gene expression, samples in columns, genes in rows
#' @param Yvec Vector of phenotype, strings, 1, 2, etc
#' @param testRes List of previously calculated test results, from testFun(.)
#' @param cluster cluster name, string '1', as in Yvec
#' @param ptail Binary phenotype vector.
#' @param breakVec vector of break points, used to bin expression data
#' @return List of Xbin and Ybin, the binned, subset, and binarized values.
#' @examples
#' mod1 <- trainDataProc(Xmat, Yvec, ptail, cluster, breakVec==c(0, 0.25, 0.5, 0.75, 0.85, 1.0))
#'
trainDataProc <- function(Xmat, Yvec, cluster=1, dtype='continuous', ptail=0.05, breakVec=c(0, 0.25, 0.5, 0.75, 1.0)) {

  Ybin <- ifelse(Yvec == cluster, yes = 1, no=0)

  if (dtype =='continuous') {
    testRes <- apply(Xmat, 1, FUN=function(a) testFun(a,Ybin))
    Xscl <- scale(Xmat) # scale each sample, in columns
    Xbinned <- apply(Xscl, 2, breakBin, breakVec) # bin each column
    rownames(Xbinned) <- rownames(Xmat)
    Xfeat <- featureSelection(Xbinned, Ybin, testRes, ptail)  # subset genes
    Xbin <- t(Xfeat$Xsub)
    genes <- Xfeat$Genes
    return(list(dat=list(Xbin=Xbin,Ybin=Ybin,Genes=genes), testRes=testRes, breakVec=breakVec))
  }
  else if (dtype =='binary') {
    testRes <- apply(Xmat, 1, FUN=function(a) testBinFun(a,Ybin))
    Xfeat <- featureSelection(Xmat, Ybin, testRes, ptail)  # subset genes
    Xbin <- t(Xfeat$Xsub)
    genes <- Xfeat$Genes
    return(list(dat=list(Xbin=Xbin,Ybin=Ybin,Genes=genes), testRes=testRes, breakVec=breakVec))
  }
  else {
    print('ERROR')
    return(NA)
  }

}


#' dataProc
#' Data preprocessing
#' @export
#' @param Xmat Matrix of gene expression, genes in columns, samples in rows
#' @param mods a model or list of models, containing breakpoints, used to bin expression data
#' @param ci the cluster label and index into list of models
#' @param dtype data type, continuous or binary
#' @param mtype model type, binned or pairs
#' @return Xbin, the binned, subset, and binarized values.
#' @examples
#' mod1 <- dataProc(X, mods)
#'
dataProc <- function(X, mods=NULL, ci=NA, dtype = 'continous', mtype = 'pairs') {

  Xmat <- as.matrix(X)

  if (length(mods) > 3) {
    breakVec <- mods[[ci]]$breakVec
    genes    <- mods[[ci]]$genes
  } else {
    breakVec <- mods$breakVec
    genes    <- mods$genes
  }

  if (dtype == 'continuous' & mtype == 'binned') {
      # here we have continuous expression values and we're doing the binned model
      # so we need to bin the data
      Xscl <- scale(Xmat) # scale each sample, in columns
      Xbin <- apply(Xscl, 2, breakBin, breakVec)
      rownames(Xbin) <- rownames(X)
      idx <-  match(table = rownames(X), x = genes)
      Xbin <- t(Xbin[idx,])
      colnames(Xbin) <- genes
  } else if (dtype == 'continuous' & mtype == 'pairs') {
      # here we have expression data, and we're using the pairs model
      # so we need to make pairs features.
      Xbin <- createPairsFeatures(X, genes)
  } else if (dtype == 'binary') {
    # here we already have pairs features.
    idx <-  match(table = rownames(X), x = genes)
    Xbin <- t(X[idx,])
    colnames(Xbin) <- genes
  } else {
    print('Error: dataProc()')
    Xbin <- NA
  }

  return(Xbin)
}

