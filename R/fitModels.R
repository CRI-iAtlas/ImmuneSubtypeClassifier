
#' fitAllModels
#' Training all subtype models.
#' @export
#' @param x Gene expression matrix.
#' @param y Phenotype vector.
#' @return A list containing xgboost classifiers.
#' @examples
#' mod1 <- fitAllModels(ebppGeneExpr, phenotype)
#'
fitAllModels <- function(x,y){
  return(TRUE)
}


#' fitOneModel
#' Train a single subtype model.
#' @export
#' @param Xbin Gene expression matrix.
#' @param Ybin Phenotype vector.
#' @param params Params for xgboost.
#' @return A single xgboost classifier.
#' @examples
#' modC1 <- fitOneModel(ebppGeneExpr, phenotype)
#'
fitOneModel <- function(Xbin, Ybin, params=list(max_depth = 2, eta = 0.5, nrounds = 33, nthread = 5)){
  bst <- xgboost(data = Xbin,
                 label = Ybin,
                 max_depth=params$max_depth,
                 eta=params$eta,
                 nrounds=params$nrounds,
                 nthread=params$nthread,
                 objective="binary:logistic")
  return(bst)
}


#' cvfitOneModel
#' Train a single subtype model using cross validation
#' @export
#' @param Xbin Gene expression matrix.
#' @param Ybin Phenotype vector.
#' @return A single xgboost classifier.
#' @examples
#' modC1 <- fitOneModel(ebppGeneExpr, phenotype)
#'
cvFitOneModel <- function(Xbin, Ybin,
                          params=list(max_depth = 2, eta = 0.5, nrounds = 100, nthread = 5, nfold=5)){
  dtrain <- xgb.DMatrix(Xbin, label = Ybin)
  cvRes <-xgb.cv(data = dtrain,
                 nrounds=params$nrounds,
                 nthread=params$nthread,
                 nfold=params$nfold,
                 max_depth=params$max_depth,
                 eta=params$eta,
                 early_stopping_rounds=2,
                 metrics = list("error", "rmse","auc"),
                 objective = "binary:logistic")

    bst <- xgboost(data = Xbin,
                 label = Ybin,
                 max_depth=params$max_depth,
                 eta=params$eta,
                 nrounds = cvRes$best_iteration,
                 nthread=params$nthread,
                 objective = "binary:logistic")
  return(bst)
}
