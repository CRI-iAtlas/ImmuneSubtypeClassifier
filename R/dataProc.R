

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
  # collect results here
  resList <- list()
  # all pairs of genes here
  all.combos <- t(combn(genes,2))

  for (i in 1:nrow(all.combos)) {
    gi <- all.combos[i,1]  # the pivot gene
    gval <- as.numeric(Xsub[gi,]) # and it's value
    # get the paired genes out.
    gcom <- all.combos[all.combos[,1] == gi,2]
    # pick sample j, get pivot value j which is for gene i, and pivot all the paired genes
    res0 <- lapply(1:ncol(Xsub), function(j) binaryGene(gval[j], Xsub[gcom,j]))
    # make matrix of features.
    resMat <- do.call('cbind', res0)
    rownames(resMat) <- sapply(gcom, function(gj) paste0(gi,':',gj))
    resList[[gi]] <- resMat
  }
  Xbin <- do.call('rbind', resList)
  colnames(Xbin) <- colnames(Xsub)
  return(Xbin)
}


#' createPairsFeatures
#' given a matrix and gene pairs, create binary features
#' @param X One gene expression matrix
#' @param genePairs a vector of gene pairs
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
  return(newMat)
}


makeSetData <- function(Xmat) {

  data("geneSetSymbols")
  resultList <- list()

  featureNames <- c()
  for (j1 in 1:4) {
    for (j2 in (j1+1):5) {
      featureNames <- c(featureNames, paste0('s',j1,'s',j2))
    }
  }

  # for each sample
  for (i in 1:ncol(Xmat)) {
    res0 <- numeric(length=10)
    idx <- 1
    for (j1 in 1:4) {
      for (j2 in (j1+1):5) {
        set1 <- genesetsymbols[[j1]]
        set2 <- genesetsymbols[[j2]]
        vals1 <- Xmat[rownames(Xmat) %in% set1,i]
        vals2 <- Xmat[rownames(Xmat) %in% set2,i]
        res1 <- sapply(vals1, function(v1) sum(v1 > vals2, na.rm=T))
        res0[idx] <- sum(res1, na.rm = T) / (length(vals1) * length(vals2))
        idx <- idx+1
      }
    }
    resultList[[i]] <- as.numeric(res0)
  }
  resMat <- do.call(cbind, resultList)
  colnames(resMat) <- colnames(Xmat)
  rownames(resMat) <- featureNames
  return(resMat)
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

  # working with matrices
  Xmat <- as.matrix(X)

  if (length(mods) > 3) {
    mods <- mods[[ci]]
  }

  # get out the relevant items
  breakVec <- mods$breakVec
  genes    <- mods$bst$feature_names
  singleGenes <- genes[!str_detect(genes, ':')]
  singleGenes <- singleGenes[!singleGenes %in% c("s1s2","s1s3","s1s4","s1s5","s2s3","s2s4","s2s5","s3s4","s3s5","s4s5")]
  pairedGenes <- genes[str_detect(genes, ':')]

  # bin the expression data
  Xbinned <- apply(Xmat, 2, breakBin, breakVec)
  rownames(Xbinned) <- rownames(Xmat)

  # and subset the genes to those not in pairs
  Xbinned <- Xbinned[rownames(Xbinned) %in% singleGenes,]

  # here we have expression data, and we're using the pairs model
  # so we need to make pairs features.
  Xpairs <- createPairsFeatures(Xmat, pairedGenes)
  colnames(Xpairs) <- colnames(Xmat)

  # gene set features.
  Xset <- makeSetData(Xmat)

  # join the data types and transpose
  Xbin <- t(rbind(Xbinned, Xpairs, Xset))
  return(Xbin)
}



#' trainDataProc
#' Data preprocessing
#' @export
#' @param Xmat Matrix of gene expression, samples in columns, genes in rows
#' @param Yvec Vector of phenotype, strings, 1, 2, etc
#' @param testRes List of previously calculated test results, from testFun(.)
#' @param subtype cluster name, string '1', as in Yvec
#' @param ptail Binary phenotype vector.
#' @param breakVec vector of break points, used to bin expression data
#' @return List of Xbin and Ybin, the binned, subset, and binarized values.
#' @examples
#' mod1 <- trainDataProc(Xmat, Yvec, ptail, cluster, breakVec==c(0, 0.25, 0.5, 0.75, 0.85, 1.0))
#'
trainDataProc <- function(Xmat, Yvec, subtype=1, ptail=0.01, breakVec=c(0, 0.25, 0.5, 0.75, 1.0)) {

  # Create the binary subtype identifier
  Ybin <- ifelse(Yvec == subtype, yes = 1, no=0)

  # bin the expression data
  Xbinned <- apply(Xmat, 2, breakBin, breakVec) # bin each column
  rownames(Xbinned) <- rownames(Xmat)

  # rank the data, and use it for feature selection
  Xrank <- apply(Xmat, 2, rank)
  testRes <- sapply(1:nrow(Xrank), function(gi) testFun(as.numeric(Xrank[gi,]), Ybin))  # get genes with big rank diffs.

  # subset the expression data for pairs
  Xfeat <- featureSelection(Xmat, Ybin, testRes, ptail)  # subset genes
  Xpairs <- makeGenePairs(Xfeat$Genes, Xfeat$Xsub)

  # subset the binned genes
  Xbinned <- Xbinned[Xfeat$Genes,]

  # gene set features.
  Xset <- makeSetData(Xmat)

  # join the data types and transpose
  Xbin <- t(rbind(Xbinned, Xpairs, Xset))
  genes <- colnames(Xbin)

  return(list(dat=list(Xbin=Xbin,Ybin=Ybin,Genes=genes), testRes=testRes, breakVec=breakVec))
}

