


#' modelPerf
#' Checking model performance.
#' @export
#' @param bst A trained xgboost
#' @param Xbin Binned gene expression data
#' @param Ybin Binary phenotype vector.
#' @return list of err, the error in prediction, and a rocPlot
#' @examples
#' mod1 <- fitAllModels(ebppGeneExpr, phenotype)
#'
modelPerf <- function(bst, Xbin, Ybin, title='perf1') {
  pred <- predict(bst, Xbin)
  err <- mean(as.numeric(pred > 0.5) != Ybin)
  print(table((pred > 0.5), Ybin))

  df <- data.frame(predictions=pred, labels=Ybin, stringsAsFactors = F)
  rocplot <- ggplot(df, aes(m = predictions, d = labels))+ geom_roc(n.cuts=20,labels=FALSE)
  rocplot <- rocplot + style_roc(theme = theme_grey) + geom_rocci(fill="pink") + ggtitle(title)
  return(list(modelError=err, plot=rocplot))
}


#' Checking model performance when Y is multiple valued
#' @export
#' @param calls Calls from callSubtypes(.)
#' @param Ytest Multi-class phenotype vector.
#' @param subtype The subtype label
#' @return list of err, the error in prediction, and a rocPlot
#' @examples
#' mod1 <- fitAllModels(ebppGeneExpr, phenotype)
#'
modelPerf2 <- function(calls, Ytest, subtype='1') {

  pred <- calls[, which(names(calls) == subtype)]
  Ybin <- sapply(Ytest, function(a) if (a == subtype){1} else {0})

  err <- mean(as.numeric(pred > 0.5) != Ybin)

  df <- data.frame(predictions=pred, labels=Ybin, stringsAsFactors = F)
  rocplot <- ggplot(df, aes(m = predictions, d = labels))+ geom_roc(n.cuts=20,labels=FALSE)
  rocplot <- rocplot + style_roc(theme = theme_grey) + geom_rocci(fill="pink") + ggtitle(paste0('Subtype: ', subtype))
  return(list(modelError=err, plot=rocplot, confusionMatrix=table((pred > 0.5), Ybin)))
}


#' Checking model performance.
#' @export
#' @param mods list of models fit to each subtype
#' @param Ytest Multi-class phenotype vector.
#' @return list of lists of err, the error in prediction, and a rocPlot
#' @examples
#' mod1 <- fitAllModels(ebppGeneExpr, phenotype)
#'
subtypePerf <- function(calls, Ytest) {

  plotList <- lapply(unique(Ytest), function(mi) {
    modelPerf2(calls, Ytest, mi)
  })

  return(plotList)
}


