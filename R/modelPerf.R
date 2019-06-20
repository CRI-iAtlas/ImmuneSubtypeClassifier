


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
modelPerf2 <- function(calls, Ytest, subtype=NA) {

  pred <- calls[, which(names(calls) == subtype)]
  Ybin <- ifelse(Ytest == subtype, yes = 1, no=0)

  err <- mean(as.numeric(pred > 0.5) != Ybin)

  df <- data.frame(predictions=pred, labels=Ybin, stringsAsFactors = F)

  baseplot <- ggplot(df, aes(m = predictions, d = labels))+ geom_roc(n.cuts=20,labels=FALSE)
  rocplot <- baseplot + annotate("text", x = .75, y = .25, label = paste("AUC =", round(calc_auc(baseplot)$AUC, 2)))
  rocplot <- rocplot + style_roc(theme = theme_grey) + geom_rocci(fill="pink") + ggtitle(paste0('Subtype: ', subtype))
  return(list(modelAUC=round(calc_auc(baseplot)$AUC, 2), modelError=err, plot=rocplot, confusionMatrix=table((pred > 0.5))))
}


#' Checking model performance.
#' @export
#' @param mods list of models fit to each subtype
#' @param Ytest Multi-class vector of subtype labels.
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



#' Build the importance table.
#' @export
#' @param mods list of models fit to each subtype
#' @param type the type of models, single, subtypes(list of 6), ensemble
#' @return table of informative features
#' @examples
#' mod1 <- getImportance(mods, type='ensemble')
#'
importantFeatures <- function(mods, mtype='ensemble', ci=1) {

  if (mtype == 'ensemble') {
    el <- 1:length(mods) # size of ensemble
    si <- 1:6           # for each subtype
  } else if (mtype == 'subtypelist') {
    el <- 1
    si <- 1:6
  } else if (mtype == 'single') {
    el <- 1
    si <- ci
  } else {
    print('Please use mtype: single, list, or ensemble')
    return(NA)
  }

  res0 <- list()
  for (j in si) {  #-- length of number of subtypes
    res1 <- list()
    for (i in el) {   #-- size of ensemble
      print(paste0(j, '  ', i))
      ei <- mods[[i]]  # get the ensemble member out, list of 6
      m <- ei[[j]]$bst    # get the subtype classifier out
      g <- xgboost::xgb.importance(model=m)
      print(head(g$Feature))
      res1[[i]]  <- g
    }
    res0[[j]] <- res1
  }

  allgenes <- data.frame()

  for (j in si) {
    print(j)
    for (i in el) {
      x <- res0[[j]][[i]]
      for (xi in 1:nrow(x)) {
        allgenes <- rbind(allgenes, data.frame(Subtype1=j, EnsembleMember=i, GeneNum=xi, Gene=x$Feature[xi], Gain=x$Gain[xi]))
      }
    }
  }

  return(allgenes)
}





