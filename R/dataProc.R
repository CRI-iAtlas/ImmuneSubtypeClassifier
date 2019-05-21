
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
featureSelection <- function(Xmat, Ybin, testRes, tail=0.05) {
  idx <- which( (testRes < quantile(testRes, tail)) | (testRes > quantile(testRes, 1.0-tail)) )
  Xsub <- Xmat[idx,]
  Xsub[is.na(Xsub)] <- 0
  return(Xsub)
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
  brks <- quantile(as.numeric(x), probs=breakVec)
  xbin <- .bincode(x = x, breaks = brks)
  xbin <- as.numeric(xbin)
  xbin
}

#' trainDataProc
#' Data preprocessing
#' @export
#' @param Xmat Matrix of gene expression, samples in columns, genes in rows
#' @param Yvec Vector of phenotype, strings, 1, 2, etc
#' @param testRes List of previously calculated test results, from testFun(.)
#' @param cluster cluster name, string '1', as in Yvec
#' @param tail Binary phenotype vector.
#' @param breakVec vector of break points, used to bin expression data
#' @return List of Xbin and Ybin, the binned, subset, and binarized values.
#' @examples
#' mod1 <- trainDataProc(Xmat, Yvec, tail, cluster, breakVec==c(0, 0.25, 0.5, 0.75, 0.85, 1.0))
#'
trainDataProc <- function(Xmat, Yvec, testRes=NULL, cluster='1', tail=0.05, breakVec=c(0, 0.25, 0.5, 0.75, 0.85, 1.0)) {
  # check that X is a matrix
  Ybin <- ifelse(Yvec == cluster, yes = 1, no=0)
  if (is.null(testRes)) {
    testRes <- apply(Xmat, 1, FUN=function(a) testFun(a,Ybin))
  }
  Xscl <- scale(Xsub) # scale each sample, in columns
  Xbinned <- apply(Xscl, 2, breakBin, breakVec)
  # changed from:  feature select, then get bin breaks over all matrix, then bin.
  Xbin <- featureSelection(Xbinned, Ybin, testRes, 0.05)
  return(list(dat=list(Xbin=Xbin,Ybin=Ybin), testRes=testRes))
}


#' testDataProc
#' Data preprocessing
#' @export
#' @param Xmat Matrix of gene expression, genes in columns, samples in rows
#' @param breakVec vector of break points, used to bin expression data
#' @return List of Xbin and Ybin, the binned, subset, and binarized values.
#' @examples
#' mod1 <- testDataProc(Xmat, breakVec==c(0, 0.25, 0.5, 0.75, 0.85, 1.0))
#'
testDataProc <- function(Xmat, breakVec=c(0, 0.25, 0.5, 0.75, 0.85, 1.0)) {
  # check that X is a matrix
  Xscl <- scale(Xmat) # scale each sample, in columns
  Xbin <- apply(Xscl, 2, breakBin, breakVec)
  # changed from:  feature select, then get bin breaks over all matrix, then bin.
  #Xbin <- featureSelection(Xbinned, Ybin, testRes, 0.05)
  return(list(dat=list(Xbin=Xbin,Ybin=Ybin), testRes=testRes))
}



