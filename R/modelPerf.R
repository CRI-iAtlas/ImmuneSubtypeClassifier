


#' modelPerf
#' Checking model performance.
#' @export
#' @param bst A trained xgboost
#' @param Xbin Binned gene expression data
#' @param Ybin Binary phenotype vector.
#' @return err, the error in prediction
#' @examples
#' mod1 <- fitAllModels(ebppGeneExpr, phenotype)
#'
modelPerf <- function(bst, Xbin, Ybin, cluster) {
  pred <- predict(bst, Xbin)
  err <- mean(as.numeric(pred > 0.5) != Ybin)
  print(table((pred > 0.5), Ybin))

  df <- data.frame(predictions=pred, labels=Ybin)
  rocplot <- ggplot(df, aes(m = predictions, d = labels))+ geom_roc(n.cuts=20,labels=FALSE)
  rocplot <- rocplot + style_roc(theme = theme_grey) + geom_rocci(fill="pink") + ggtitle(paste0('C', cluster, ', ROC'))
  return(list(modelError=err, plot=rocplot))
}

