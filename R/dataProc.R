
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
trainDataProc <- function(Xmat, Yvec, testRes=NULL, cores=2, cluster=1, dtype='continuous', ptail=0.05, breakVec=c(0, 0.25, 0.5, 0.75, 1.0)) {

  Ybin <- ifelse(Yvec == cluster, yes = 1, no=0)


  if (dtype =='continuous' & is.null(testRes)) {
    testRes <- apply(Xmat, 1, FUN=function(a) testFun(a,Ybin))
    Xscl <- scale(Xmat) # scale each sample, in columns
    Xbinned <- apply(Xscl, 2, breakBin, breakVec) # bin each column
    rownames(Xbinned) <- rownames(Xmat)
    Xfeat <- featureSelection(Xbinned, Ybin, testRes, ptail)  # subset genes
    Xbin <- t(Xfeat$Xsub)
    genes <- Xfeat$Genes
    return(list(dat=list(Xbin=Xbin,Ybin=Ybin,Genes=genes), testRes=testRes, breakVec=breakVec))
  }
  else if (dtype =='binary' & is.null(testRes)) {
    testRes <- apply(Xmat, 1, FUN=function(a) testBinFun(a,Ybin))
    Xfeat <- featureSelection(Xmat, Ybin, testRes, ptail)  # subset genes
    Xbin <- t(Xfeat$Xsub)
    genes <- Xfeat$Genes
    return(list(dat=list(Xbin=Xbin,Ybin=Ybin,Genes=genes), testRes=testRes, breakVec=breakVec))
  }


}


#' dataProc
#' Data preprocessing
#' @export
#' @param Xmat Matrix of gene expression, genes in columns, samples in rows
#' @param mods a model or list of models, containing breakpoints, used to bin expression data
#' @param ci the cluster label and index into list of models
#' @return Xbin, the binned, subset, and binarized values.
#' @examples
#' mod1 <- dataProc(X, mods)
#'
dataProc <- function(X, mods=NULL, ci=NA) {

  Xmat <- as.matrix(X)

  if (length(mods) > 3) {
    breakVec <- mods[[ci]]$breakVec
    genes    <- mods[[ci]]$genes
  } else {
    breakVec <- mods$breakVec
    genes    <- mods$genes
  }

  Xscl <- scale(Xmat) # scale each sample, in columns
  Xbin <- apply(Xscl, 2, breakBin, breakVec)
  rownames(Xbin) <- rownames(X)
  idx <-  match(table = rownames(X), x = genes)
  Xbin <- t(Xbin[idx,])
  colnames(Xbin) <- genes
  return(Xbin)
}



