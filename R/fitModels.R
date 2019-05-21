
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
#' @return A single xgboost classifier.
#' @examples
#' modC1 <- fitOneModel(ebppGeneExpr, phenotype)
#'
fitOneModel <- function(Xbin, Ybin){
  bst <- xgboost(data = Xbin, label = Ybin, max_depth = 2, eta = 0.5, nrounds = 33, nthread = 5, objective = "binary:logistic")
  return(bst)
}
