#' @import data.table
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise arrange desc
#' @importFrom xgboost xgb.importance
NULL



#' Get a list of genes that are used in the model.
#' @param model A trained robencla model object. If NULL (default), loads
#'   the package's built-in model.
#' @param model_path Character string, path to a saved robencla model file
#'   (.rda format). Used only if \code{model} is NULL and a custom model
#'   path is desired. Default is NULL (uses built-in model).
#' @param geneid Character string specifying the gene ID type in row names.
#'   One of "symbol", "entrez", "ensembl", or "pairs". Default is "symbol".
#'
#' @return A list containing:
#'   \item{model_genes}{List of the genes from the model object.}
#'   \item{gene_mape}{data.frame, TCGA PanCancer EB++ genes that match the model_genes.}
#'
#' @details
#' The function loads the \code{ebpp_genes_sig} reference data containing
#' the gene identifiers used during model training. Model genes are matched
#' to this reference, and gene table is returned..
#'
#' @examples
#' \dontrun{
#' # Match by gene symbol
#' result <- modelGenes(model=NULL, model_path='immune_optimized_99_pairs.rds')
#' gene_list <- result$model_genes
#' matching_gene_table <- result$gene_map
#' }
#'
#' @export
modelGenes <- function(model=NULL, model_path=NULL, geneid='symbol') {

  # get the required genes from the model
  if (!is.null(model)) {
    model_genes <- model$pair_list
  } else if (!is.null(model_path)) {
    model <- readRDS(model_path)
    # get the model genes out
    model_genes <- model$pair_list
  } else {
    print("Error: geneMatch ... Please include a model!")
    return(NULL)
  }

  # unlist and unique
  model_genes <- unique(as.vector(unlist(model_genes)))

  # get the EBpp gene table
  data(ebpp_gene, envir = environment())

  # get the gene map for these genes
  # get the gene map for these genes
  if(geneid == 'symbol') {
    gene_map <- ebpp_genes_full[ebpp_genes_full$Symbol %in% model_genes, ]
  } else if (geneid == 'entrez') {
    gene_map <- ebpp_genes_full[ebpp_genes_full$Entrez %in% model_genes, ]
  } else if (geneid == 'ensembl') {
    gene_map <- ebpp_genes_full[ebpp_genes_full$Ensembl %in% model_genes, ]
  } else {
    stop("For geneid, please use: symbol, entrez, ensembl, or pairs")
  }

  if(nrow(gene_map) < length(model_genes)) {
    print("****** WARNING ******")
    print("Number of genes in model doesn't match number in map.")
  }

  return(list(model_genes=model_genes, gene_map=gene_map))
}


#' Match Gene IDs to Model Features
#'
#' Match incoming data gene identifiers to features used in training.
#'
#' @param X data.frame with samples in rows, genes in columns.
#' @param model A trained robencla model object.
#' @param model_path Path to a saved robencla model .rds file.
#' @param geneid One of "symbol", "entrez", "ensembl", or "pairs".
#' @param sampleid Column name for sample IDs. Default "SampleBarcode".
#' @param labelid Column name for labels. Default NULL.
#' @param error_limit Proportion of missing genes tolerated (0.0-1.0). Default 0.0.
#'
#' @return List with Subset (matched data.frame), matchError (proportion missing),
#'   and missingGenes (character vector of unmatched genes).
#'
#' @export
geneMatch <- function(X,
                      model = NULL,
                      model_path = NULL,
                      geneid = "symbol",
                      sampleid = "SampleBarcode",
                      labelid = "Label",
                      error_limit = 0.0) {

  X <- as.data.frame(X)

  # Load model if needed
  if (is.null(model)) {
    if (!is.null(model_path)) {
      model <- readRDS(model_path)
    } else {
      stop("geneMatch: please provide either model or model_path")
    }
  }

  # pairs mode - no matching needed
  if (geneid == "pairs") {
    return(list(Subset = X, matchError = 0, missingGenes = character(0)))
  }

  # Get unique model genes
  modelgenes <- unique(as.vector(unlist(model$pair_list)))

  # Load reference gene table
  data("ebpp_gene", package = "ImmuneSubtypeClassifier", envir = environment())

  # Build gene map based on ID type
  gene_map <- switch(geneid,
                     symbol  = ebpp_genes_full[ebpp_genes_full$Symbol %in% modelgenes, ],
                     entrez  = ebpp_genes_full[ebpp_genes_full$Entrez  %in% modelgenes, ],
                     ensembl = ebpp_genes_full[ebpp_genes_full$Ensembl %in% modelgenes, ],
                     stop("geneid must be one of: symbol, entrez, ensembl, pairs")
  )

  # Match gene map to input columns
  if (geneid == "ensembl") {
    # Strip version suffixes from input (e.g. ENSG00000001.12 -> ENSG00000001)
    input_genes <- vapply(
      strsplit(colnames(X), "\\."),
      function(a) a[1],
      character(1)
    )
  } else {
    input_genes <- colnames(X)
  }

  id_col <- switch(geneid,
                   symbol  = gene_map$Symbol,
                   entrez  = gene_map$Entrez,
                   ensembl = gene_map$Ensembl
  )

  idx <- match(id_col, input_genes)

  # Report and handle missing genes
  matchError <- sum(is.na(idx)) / length(modelgenes)
  missingGenes <- modelgenes[!(modelgenes %in% id_col[!is.na(idx)])]

  if (matchError > 0.0) {
    message("**************************************")
    message("    Gene Match Error Report           ")
    message(sprintf("  percent missing genes: %.1f%%", matchError * 100))
    message("  Missing genes:")
    message(paste(" ", missingGenes, collapse = "\n"))
    message("  Filling missing genes with NA.")
    message("**************************************")
  }

  # Hard stop for catastrophic missingness
  if (matchError > error_limit) {
    stop(sprintf(
      "geneMatch: %.1f%% of model genes missing — too many to impute. Check geneid type and data orientation.\nMissing genes: %s",
      matchError * 100,
      paste(missingGenes, collapse=", ")
    ))
  }

  # Subset to matched gene columns
  matched_idx <- idx[!is.na(idx)]
  X2 <- X[, matched_idx, drop = FALSE]
  colnames(X2) <- gene_map$Symbol[!is.na(idx)]

  # Fill missing genes with NA
  if (length(missingGenes) > 0) {
    for (g in missingGenes) {
      X2[[g]] <- NA_real_
    }
  }
  # Preserve sampleid column if present
  if (!is.null(sampleid) && sampleid %in% colnames(X)) {
    X2[[sampleid]] <- X[[sampleid]]
  }

  # Preserve labelid column if present
  if (!is.null(labelid) && labelid %in% colnames(X)) {
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
#' @param X_or_path Gene expression matrix (or path) with genes in rows and samples in columns.
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
callSubtypes <- function(X_or_path=NULL,
                         model = NULL,
                         model_path = NULL,
                         geneid = "symbol",
                         sampleid = 'Barcode',
                         labelid=NULL,
                         error_limit = 0.2) {
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

  if(typeof(X_or_path) == 'character') {
    X <- readr::read_csv(X_or_path)
  } else {
    X <- X_or_path
  }

  # check that all needed genes are available
  res0 <- geneMatch(X, model, model_path,
                   geneid=geneid, sampleid=sampleid, labelid=labelid,
                   error_limit=error_limit)

  X <- res0$Subset

  if (is.null(X)) {
  stop(sprintf(
    "callSubtypes: gene matching failed (%.1f%% missing). See above for details.",
    res0$matchError * 100
  ))
}

  print("Starting prediction...")
  model$predict(
    data_frame = X,
    label_name = labelid,
    sample_id = sampleid
  )
  print("...finshed prediction")

  results <- model$results()

  # Defensive check on results structure
  if (is.null(results) || nrow(results) == 0) {
    stop("callSubtypes: model$results() returned empty or NULL — prediction may have failed")
  }

  # Find the sample ID column flexibly
  possible_id_cols <- c("SampleID", "SampleIDs", "Barcode", "SampleBarcode", sampleid)
  id_col <- intersect(possible_id_cols, colnames(results))[1]
  if (is.na(id_col)) {
    stop(sprintf(
      "callSubtypes: could not find sample ID column in results. Available columns: %s",
      paste(colnames(results), collapse=", ")
    ))
  }

  output <- data.frame(
    SampleIDs = results[[id_col]],
    BestCall  = as.integer(gsub("C", "", results$BestCall)),
    Label     = if (!is.null(labelid)) as.integer(gsub("C", "", results$Label)) else NA,
    stringsAsFactors = FALSE
  )

  score_cols <- grep("^C[1-6]$", colnames(results), value = TRUE)

  if (length(score_cols) > 0) {
    scores <- results[, score_cols, drop = FALSE]
    colnames(scores) <- gsub("C", "", colnames(scores))
    scores <- scores[, as.character(1:6)]
    output <- cbind(output, scores)
  }

  return(list(Model=model, Pred=output))
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
