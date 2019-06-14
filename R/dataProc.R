

#' testFun
#' Get difference in mean rank sums for a single gene
#'
#' @param rankg Gene expression profile for single sample
#' @param Ybin Binary phenotype vector.
#' @return test result, numeric value, the rank sum test.
#' @examples
#' res1 <- testFun(G, Ybin)
#' @export
testFun <- function(rankg, Ybin) {
  res0 <- (sum(rankg[Ybin == 0], na.rm = T) / sum(Ybin == 0, na.rm = T)) - (sum(rankg[Ybin == 1], na.rm = T) / sum(Ybin == 1, na.rm = T))
  return(res0)
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



binaryGene <- function(pivotvalue, values) {
  res0 <- sapply(values, function(b) as.numeric(b <= pivotvalue))
  res0[is.na(res0)] <- rbinom(n = sum(is.na(res0)), prob = 0.5, 1)  ## replace NAs with random values
  return(res0)
}


makeGenePairs <- function(genes, Xsub) {

  # for each gene
  resList <- list()
  for (gi in genes) {
    # do pairs
    gval <- as.numeric(Xsub[gi,])
    res0 <- lapply(1:ncol(Xsub), function(i) binaryGene(gval[i], Xsub[,i]))
    # make matrix of features.
    resMat <- do.call('cbind', res0)
    rownames(resMat) <- sapply(genes, function(gj) paste0(gi,':',gj))
    resList[[gi]] <- resMat
  }
  Xbin <- do.call('rbind', resList)
  colnames(Xbin) <- colnames(Xsub)
  return(Xbin)
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

  # first convert the gene pairs to a named list
  # where each entry of the list is the genes for a given pivot-gene
  genePairs <- strsplit(genes, ':')
  pairList <- list()
  for(gi in genePairs) {
    pairList[[gi[1]]] <- c(pairList[[gi[1]]], gi[2])
  }

  resList <- list() # then for each pivot gene
  for (gi in names(pairList)) {
    # assuming it's in the data ... really should be!
    if (gi %in% rownames(X)) {
      gs <- pairList[[gi]]                  ## can end up with the pivot gene in the genes...
      pval <- as.numeric(X[gi,])            ## pivot values across samples
      idx <- match(table=rownames(X), x=gs) ## get index to genes for this pivot
      Xsub <- X[idx,]                       ## subset the matrix, NAs for missing genes, pivot gene on top
      if (class(Xsub) == 'numeric' & length(gs) == 1) {
        Xsub <- matrix(data=Xsub, ncol=ncol(X), nrow=1)
        colnames(Xsub) <- colnames(X)
      }
      rownames(Xsub) <- gs                  ## give gene IDs
      res0 <- lapply(1:ncol(Xsub), function(a) binaryGene(pval[a], Xsub[,a]))  ## create binary values
      resList[[gi]] <- do.call('cbind', res0)
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
trainDataProc <- function(Xmat, Yvec, cluster=1, dtype='continuous', ptail=0.01, breakVec=c(0, 0.25, 0.5, 0.75, 1.0), mtype='pairs') {

  Ybin <- ifelse(Yvec == cluster, yes = 1, no=0)

  if (dtype =='continuous' & mtype == 'binned') {
    Xrank <- apply(Xmat, 2, rank)
    testRes <- apply(Xrank, 1, FUN=function(a) testFun(a,Ybin))
    Xscl <- scale(Xmat) # scale each sample, in columns
    Xbinned <- apply(Xscl, 2, breakBin, breakVec) # bin each column
    rownames(Xbinned) <- rownames(Xmat)
    Xfeat <- featureSelection(Xbinned, Ybin, testRes, ptail)  # subset genes
    Xbin <- t(Xfeat$Xsub)
    genes <- Xfeat$Genes
    return(list(dat=list(Xbin=Xbin,Ybin=Ybin,Genes=genes), testRes=testRes, breakVec=breakVec))
  } else if (dtype =='continuous' & mtype == 'pairs') {
    Xrank <- apply(Xmat, 2, rank) # rank the expression within sample
    testRes <- sapply(1:nrow(Xrank), function(gi) testFun(as.numeric(Xrank[gi,]), Ybin))  # get genes with big rank diffs.
    Xfeat <- featureSelection(Xmat, Ybin, testRes, ptail)  # subset genes
    Xbin <- makeGenePairs(Xfeat$Genes, Xfeat$Xsub)
    Xbin <- t(Xbin)
    genes <- colnames(Xbin)
    return(list(dat=list(Xbin=Xbin,Ybin=Ybin,Genes=genes), testRes=testRes, breakVec='pairs'))
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
dataProc <- function(X, mods=NULL, ci=NA, dtype = 'continuous', mtype = 'pairs') {

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

