#' @import data.table
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise arrange desc
#' @importFrom xgboost xgb.importance
NULL




#' Print Gene Match Error Report
#'
#' Display a formatted report of the gene matching error rate.
#'
#' @param err Numeric, the match error proportion (0-1).
#'
#' @return NULL (called for side effect of printing).
#'
#' @keywords internal
reportError <- function(err) {
  message("**************************************")
  message("    Gene Match Error Report           ")
  message("                                      ")
  message(sprintf("  percent missing genes: %.1f%%", err * 100))
  message("                                      ")
  message("see ?geneMatchErrorReport for details ")
  message("                                      ")
  message("**************************************")
}


#' Match Gene IDs to Model Features
#'
#' Match the incoming data gene identifiers to the features used in training.
#' Supports multiple gene ID formats.
#'
#' @param X Gene expression matrix with genes in rows and samples in columns.
#' @param geneid Character string specifying the gene ID type in the row names.
#'   One of "symbol" (HGNC symbols), "entrez" (Entrez IDs), "ensembl"
#'   (Ensembl IDs with optional version suffix), or "pairs" (already matched,
#'   returns input unchanged). Default is "pairs".
#' @param sampleid Character string specifying the sample ID column name.
#'   Default is "SampleBarcode".
#'
#' @return A list containing:
#'   \item{Subset}{Matrix with rows reordered/subset to match model features.
#'     Missing genes appear as NA rows.}
#'   \item{matchError}{Numeric, proportion of model genes not found in input
#'     (0 = perfect match, 1 = no genes matched).}
#'   \item{missingGenes}{Genes needed by the model, but not found in the gene_map or expression matrix.}
#'
#' @details
#' The function loads the \code{ebpp_genes_sig} reference data containing
#' the gene identifiers used during model training. Input genes are matched
#' to this reference, and the output matrix is reordered accordingly.
#' Genes present in the reference but missing from the input will have
#' NA values in the output.
#'
#' @examples
#' \dontrun{
#' # Match by gene symbol
#' result <- geneMatch(expr_matrix, geneid = "symbol")
#' matched_data <- result$Subset
#' pct_missing <- result$matchError
#' }
#'
#' @export
geneMatch <- function(X,
                      model,
                      model_path,
                      geneid = "symbol",
                      sampleid = "SampleBarcode",
                      labelid = "Label",
                      error_limit = 0.0) {

  # for convenience
  X <- as.data.frame(X)

  # get the required genes from the model
  if (!is.null(model)) {
    modelgenes <- model$pair_list
  } else if (!is.null(model_path)) {
    model <- readRDS(model_path)
  } else {
    print("Error: geneMatch ... Please include a model!")
    return(NULL)
  }

  # get the model genes out
  modelgenes <- unique(as.vector(unlist(model$pair_list)))

  # get the EBpp gene table
  data(ebpp_gene, envir = environment())

  # get the gene map for these genes
  gene_map <- ebpp_genes_full[ebpp_genes_full$Symbol %in% modelgenes, ]

  if (geneid == "symbol") {
    idx <- match(x = gene_map$Symbol, table = colnames(X))
  } else if (geneid == "entrez") {
    idx <- match(x = gene_map$Entrez, table = colnames(X))
  } else if (geneid == "ensembl") {
    input_ens <- stringr::str_split(colnames(X), pattern = "\\.")
    input_ens <- vapply(input_ens, function(a) a[1], character(1))
    idx <- match(x = gene_map$Ensembl, table = input_ens)
  } else if (geneid == "pairs") {
    return(list(Subset = X, matchError = 0))
  } else {
    stop("For geneid, please use: symbol, entrez, ensembl, or pairs")
  }

  # report the match error
  matchError <- sum(is.na(idx)) / (length(modelgenes))
  missingGenes <- c()
  if (matchError > error_limit) {
    # Genes required by model but not found in input data
    missingGenes <- modelgenes[is.na(idx)]
    reportError(matchError)
    cat("Missing genes:\n")
    print(missingGenes)
    return(list(Subset = NULL, matchError = matchError, missingGenes = missingGenes))
  }

  # Subset COLUMNS, but keep sample ID column
  gene_cols <- idx[!is.na(idx)]
  # Instead of X2 <- X[, gene_cols, drop = FALSE]
  X2 <- X[, gene_cols, drop = FALSE]

  colnames(X2) <- gene_map$Symbol[!is.na(idx)]

  # Preserve sample ID column
  if (sampleid %in% colnames(X)) {
    X2[[sampleid]] <- X[[sampleid]]
  }

  # Preserve Label column
  if (labelid %in% colnames(X)) {
    X2[[labelid]] <- X[[labelid]]
  }


  return(list(Subset = X2, matchError = matchError, missingGenes = missingGenes))
}


#' Call Immune Subtypes Using Robencla Model
#'
#' Make immune subtype predictions for samples using a trained robencla
#' ensemble classifier.
#'
#' @import data.table
#' @import robencla
#'
#' @param X Gene expression matrix with genes in rows and samples in columns.
#'   Row names should be gene identifiers matching the \code{geneid} parameter.
#'   Column names should be sample identifiers.
#' @param model A trained robencla model object. If NULL (default), loads
#'   the package's built-in model.
#' @param model_path Character string, path to a saved robencla model file
#'   (.rda format). Used only if \code{model} is NULL and a custom model
#'   path is desired. Default is NULL (uses built-in model).
#' @param geneid Character string specifying the gene ID type in row names.
#'   One of "symbol", "entrez", "ensembl", or "pairs". Default is "symbol".
#' @param labelid Character string specifying the sample labels.
#'   If present, then prediction metrics can be computed with model$pred
#'
#' @return A data frame with columns:
#'   \item{SampleIDs}{Sample identifiers from input column names.}
#'   \item{BestCall}{Predicted subtype (1-6).}
#'   \item{1-6}{Prediction scores for each subtype.}
#'
#' @details
#' This function wraps the robencla model's predict method to provide a
#' simple interface for immune subtype classification. It handles gene
#' matching, data transformation, and reformats the output to match the
#' legacy interface.
#'
#' The robencla model uses named feature pairs and an ensemble of XGBoost
#' classifiers to make predictions. Each sample receives a score for each
#' of the 6 immune subtypes, and the BestCall is determined by the model's
#' final prediction layer.
#'
#' @examples
#' \dontrun{
#' # Using default built-in model
#' results <- callSubtypes(expr_matrix, geneid = "symbol")
#'
#' # Using a custom model
#' my_model <- readRDS("path/to/model.rds")
#' results <- callSubtypes(expr_matrix, model = my_model, geneid = "symbol")
#'
#' # View results
#' table(results$BestCall)
#' }
#'
#' @export
callSubtypes <- function(X,
                         model = NULL,
                         model_path = NULL,
                         geneid = "symbol",
                         sampleid = 'Barcode',
                         labelid=NULL) {
  if (is.null(model)) {
    if (!is.null(model_path)) {
      model <- readRDS(model_path)
    } else {
      # Try to auto-load from models directory
      models_dir <- "model"
      if (dir.exists(models_dir)) {
        rds_files <- list.files(models_dir, pattern = "\\.rds$", full.names = TRUE)
        if (length(rds_files) == 1) {
          message("Auto-loading model from: ", rds_files[1])
          model <- readRDS(rds_files[1])
        } else if (length(rds_files) > 1) {
          message("Multiple .rds files found in ", models_dir, ". Using: ", rds_files[1])
          model <- readRDS(rds_files[1])
        } else {
          # Fallback to default
          data("robencla_model", envir = environment())
          model <- emod
        }
      } else {
        # Fallback to default
        data("robencla_model", envir = environment())
        model <- emod
      }
    }
  }

  # check that all needed genes are available
  res0 <- geneMatch(X, model, model_path,
                    geneid=geneid, sampleid=sampleid, labelid=labelid)

  X <- res0$Subset

  print("Starting prediction...")
  model$predict(
    data_frame = X,
    label_name = labelid,
    sample_id = sampleid
  )
  print("...finshed prediction")

  results <- model$results()

  output <- data.frame(
    SampleIDs = results$SampleID,
    BestCall = as.integer(gsub("C", "", results$BestCall)),
    Label = as.integer(gsub("C", "", results$Label)),
    stringsAsFactors = FALSE
  )

  score_cols <- grep("^C[1-6]$", colnames(results), value = TRUE)

  if (length(score_cols) > 0) {
    scores <- results[, score_cols, drop = FALSE]
    colnames(scores) <- gsub("C", "", colnames(scores))
    scores <- scores[, as.character(1:6)]
    output <- cbind(output, scores)
  }

  return(output)
}

#' Get Prediction Scores Only
#'
#' Extract just the prediction score matrix from subtype calls.
#'
#' @param results Data frame returned by \code{\link{callSubtypes}} or
#'   \code{\link{callEnsemble}}.
#'
#' @return A numeric matrix with samples in rows and subtypes (1-6) in columns.
#'
#' @examples
#' \dontrun{
#' results <- callSubtypes(expr_matrix)
#' scores <- getScores(results)
#' heatmap(as.matrix(scores))
#' }
#'
#' @export
getScores <- function(results) {
  score_cols <- intersect(colnames(results), as.character(1:6))
  as.matrix(results[, score_cols, drop = FALSE])
}


#' Get Best Calls Only
#'
#' Extract just the best subtype calls from prediction results.
#'
#' @param results Data frame returned by \code{\link{callSubtypes}} or
#'   \code{\link{callEnsemble}}.
#'
#' @return A character vector of subtype calls (1-6).
#'
#' @examples
#' \dontrun{
#' results <- callSubtypes(expr_matrix)
#' calls <- getBestCalls(results)
#' table(calls)
#' }
#'
#' @export
getBestCalls <- function(results) {
  results$BestCall
}
