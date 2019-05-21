
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
  testres <- (sum(rankg[y == 0]) / sum(y == 0)) - (sum(rankg[y == 1]) / sum(y == 1))
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

#' dataProc
#' Data preprocessing
#' @export
#' @param Xmat Matrix of gene expression, genes in columns, samples in rows
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
  Xsub <- featureSelection(Xmat, Ybin, testRes, 0.05)
  Xscl <- t(scale(Xsub)) # scale each sample
  brks <- quantile(as.numeric(Xscl), probs=breakVec)
  Xbin <- apply(Xscl, 2, function(x) .bincode(x = x, breaks = brks))
  Xbin <- apply(Xbin, 2, as.numeric)
  return(list(dat=list(Xbin,Ybin), testRes=testRes))
}



