

#' Train a Robencla Classifier for Immune Subtype Classification
#'
#' Trains a robencla ensemble classifier using named feature pairs and
#' XGBoost parameters tuned for immune subtype prediction.
#'
#' @param data A data.frame containing the training data with features,
#'   labels, and sample identifiers.
#' @param label_name Character string specifying the column name containing
#'   class labels. Default is "ClusterLabel".
#' @param sample_id Character string specifying the column name containing
#'   sample identifiers. Default is "SampleBarcode".
#' @param pair_list A named list where each element corresponds to a class
#'   and contains a character vector of feature names. Features are paired
#'   sequentially (1-2, 3-4, 5-6, ...) for named pair transformation.
#' @param max_depth Integer, maximum tree depth. Higher values increase
#'   model complexity. Default is 12.
#' @param eta Numeric, learning rate. Smaller values are more conservative.
#'   Default is 0.3.
#' @param nrounds Integer, maximum number of boosting rounds. Default is 64.
#' @param early_stopping_rounds Integer, training stops if performance
#'   doesn't improve for this many rounds. Default is 2.
#' @param nthreads Integer, number of parallel threads. Default is 4.
#' @param gamma Numeric, minimum loss reduction required to make a split.
#'   Higher values are more conservative. Default is 0.2.
#' @param lambda Numeric, L2 regularization term on weights. Default is 1.2.
#' @param alpha Numeric, L1 regularization term on weights. Default is 0.2.
#' @param ensemble_size Integer, number of models in each class ensemble.
#'   Default is 11.
#' @param sample_prop Numeric, proportion of samples used to train each
#'   ensemble member. Default is 0.8.
#' @param feature_prop Numeric, proportion of features used to train each
#'   ensemble member. Default is 0.8.
#' @param subsample Numeric, XGBoost subsample ratio of training instances.
#'   Default is 0.8.
#' @param combine_function Character, method for combining ensemble
#'   predictions. Currently only "median" is supported. Default is "median".
#' @param verbose Integer, verbosity level. 0 for silent. Default is 0.
#' @param trim_model Logical, whether to call trim() on the model after
#'   training to reduce object size. Default is TRUE.
#'
#' @return A trained Robencla model object.
#'
#' @details
#' The classifier uses named feature pairs (namedpairs mode) where features
#' are paired sequentially within each class's pair_list. XGBoost parameters
#' are tuned for immune subtype classification tasks.
#'
#' For more information on XGBoost parameters, see:
#' \url{https://xgboost.readthedocs.io/en/latest/parameter.html}
#'
#' @examples
#' \dontrun{
#' # Define feature pairs for each class
#' pair_list <- list(
#'   C1 = c("GeneA", "GeneB", "GeneC", "GeneD"),
#'   C2 = c("GeneE", "GeneF", "GeneG", "GeneH")
#' )
#'
#' # Train the classifier
#' model <- trainRobencla(
#'   data = training_data,
#'   label_name = "Label",
#'   sample_id = "Barcode",
#'   pair_list = pair_list
#' )
#'
#' # Make predictions
#' model$predict(data_frame = test_data,
#'               label_name = "Label",
#'               sample_id = "SampleID")
#' results <- model$results()
#' }
#'
#' @seealso \code{\link[robencla]{Robencla}}
#'
#' @export
trainRobencla <- function(data,
                          label_name = "ClusterLabel",
                          sample_id = "SampleBarcode",
                          pair_list,
                          max_depth = 12,
                          eta = 0.3,
                          nrounds = 64,
                          early_stopping_rounds = 2,
                          nthreads = 4,
                          gamma = 0.2,
                          lambda = 1.2,
                          alpha = 0.2,
                          ensemble_size = 11,
                          sample_prop = 0.8,
                          feature_prop = 0.8,
                          subsample = 0.8,
                          combine_function = "median",
                          verbose = 0,
                          trim_model = TRUE) {

  obj_name <- as.character(Sys.time())

  mod <- robencla::Robencla$new(obj_name)

  print(paste0("Robencla version: ", mod$version()))

  params <- list(
    max_depth = max_depth,
    eta = eta,
    nrounds = nrounds,
    early_stopping_rounds = early_stopping_rounds,
    nthreads = nthreads,
    gamma = gamma,
    lambda = lambda,
    alpha = alpha,
    size = ensemble_size,
    sample_prop = sample_prop,
    feature_prop = feature_prop,
    subsample = subsample,
    combine_function = combine_function,
    verbose = verbose
  )

  mod$train(
    data_frame = data,
    label_name = label_name,
    sample_id = sample_id,
    data_mode = c("namedpairs"),
    signatures = NULL,
    pair_list = pair_list,
    params = params
  )

  if (trim_model) {
    mod$trim()
  }

  return(mod)
}


#' Default Feature Pair List for Immune Subtype Classification
#'
#' Returns the default named pair list used for immune subtype classification.
#' Each class (C1-C6) has 10 feature pairs (20 genes) selected based on
#' discriminative power.
#'
#' @return A named list with elements C1 through C6, each containing a
#'   character vector of 20 gene names representing 10 feature pairs.
#'
#' @examples
#' pair_list <- getFeaturesPairList()
#' names(pair_list)
#'
#' @export
getFeaturesPairList <- function() {
  list(
    C1 = c("B2M", "COL3A1", "B2M", "COL1A2", "COL1A2", "HLA-B",
           "COL3A1", "HLA-B", "APOE", "COL6A3", "APOE", "SDC1",
           "APOE", "COL6A1", "APOE", "MMP14", "APOE", "MMP2",
           "COL1A2", "SPARC"),
    C2 = c("IFI27", "RHOB", "IFI6", "RHOB", "IFI27", "LRP1",
           "IFI27", "NPC2", "IFI6", "NPC2", "IFI27", "TAGLN",
           "RHOB", "STAT1", "IFI6", "LRP1", "IFI27", "IGFBP3",
           "IFI27", "IGFBP4"),
    C3 = c("COL3A1", "SPARC", "IGFBP4", "JUP", "CD59", "JUP",
           "NPC2", "SLC25A5", "IFI27", "NPC2", "IGFBP4", "PFN1",
           "HNRNPA2B1", "IGFBP4", "CCT5", "NPC2", "IFI27", "MET",
           "IGFBP4", "TPI1"),
    C4 = c("APOE", "COL3A1", "APOE", "COL1A2", "COL3A1", "SPARC",
           "APOC1", "DSP", "COL1A2", "SPARC", "APOC1", "JUP",
           "APOC1", "COL6A3", "APOC1", "LYZ", "APOC1", "SDC1",
           "DSP", "RHOB"),
    C5 = c("APOE", "FN1", "FN1", "SPARC", "APOE", "COL3A1",
           "APOE", "COL1A2", "B2M", "SPARC", "APOE", "B2M",
           "FN1", "SLC1A3", "COL3A1", "SLC1A3", "APOE", "HLA-B",
           "COL1A2", "SLC1A3"),
    C6 = c("B2M", "COL3A1", "B2M", "COL1A2", "COL6A3", "ENO1",
           "BSG", "COL6A3", "COL6A3", "HNRNPA2B1", "COL6A3", "MYL6",
           "BSG", "COL6A1", "COL6A3", "DSP", "BSG", "MMP2",
           "JUP", "MMP2")
  )
}



#' Default Gene List for Immune Subtype Classification
#'
#' Returns the default named genes used for immune subtype classification.
#'
#' @return A named list with the 36 used genes.
#'
#' @examples
#' gene_list <- getFeaturesGeneTable()
#'
#' @export
getFeaturesGeneTable <- function() {
  gs <- unique(unlist(getFeaturesPairList()))
  data(ebpp_gene, envir = environment())
  df <- ebpp_genes_full[ebpp_genes_full$Symbol %in% gs,]
  # 2. Convert factors to characters (Modifying 'df')
  cols_to_fix <- c("Symbol", "Entrez", "Ensembl")
  df[cols_to_fix] <- lapply(df[cols_to_fix], as.character)
  return(
    df
  )
}



#' Build and Save Robencla Classifier
#'
#' Convenience function to build a robencla classifier from a data file,
#' evaluate performance, and save the trained model.
#'
#' @import data.table
#' @param data_path Character string, path to the training data CSV file.
#' @param output_path Character string, path where the trained model will
#'   be saved as an .rda file.
#' @param label_name Character string, column name containing class labels.
#'   Default is "Label".
#' @param sample_id Character string, column name containing sample IDs.
#'   Default is "Barcode".
#' @param pair_list Named list of feature pairs. If NULL, uses
#'   \code{getFeaturesPairList()}. Default is NULL.
#' @param label_prefix Character string to prepend to numeric labels.
#'   Set to NULL to skip label transformation. Default is "C".
#' @param evaluate Logical, whether to run self-evaluation on training data
#'   and print metrics. Default is TRUE.
#' @param ... Additional arguments passed to \code{trainRobencla()}.
#'
#' @return The trained and trimmed Robencla model object (invisibly).
#'
#' @examples
#' \dontrun
#' model <- build_robencla_classifier(
#'   data_path = "data/EBpp_pancancer.csv.gz",
#'   output_path = "models/robencla_trained_model.rda",
#'   label_name = "Label",
#'   sample_id = "Barcode"
#' )
#'
#' @export
build_robencla_classifier <- function(data_path,
                                      output_path,
                                      label_name = "Label",
                                      sample_id = "Barcode",
                                      pair_list = NULL,
                                      label_prefix = "C",
                                      evaluate = TRUE,
                                      ...) {

  if (is.null(pair_list)) {
    pair_list <- getFeaturesPairList()
  }

  message("Reading training data from: ", data_path)
  data <- readr::read_csv(data_path, show_col_types = FALSE)

  if (!is.null(label_prefix)) {
    data[[label_name]] <- paste0(label_prefix, data[[label_name]])
  }

  message("Training robencla classifier...")
  model <- trainRobencla(
    data = data,
    label_name = label_name,
    sample_id = sample_id,
    pair_list = pair_list
  )

  if (evaluate) {
    message("Evaluating on training data...")
    model$predict(
      data_frame = data,
      label_name = label_name,
      sample_id = sample_id
    )

    results <- model$results()
    cat("\nConfusion Matrix:\n")
    print(table(Predicted = results$BestCalls, Actual = model$test_label))

    cat("\nClassification Metrics:\n")
    print(model$classification_metrics())
  }

  message("Saving model to: ", output_path)
  saveRDS(model, file = output_path)

  invisible(model)
}


# C1   C2   C3   C4   C5   C6
# C1 2360   28   51   17    0   11
# C2   20 2555    3   12    0    0
# C3   28    4 2321   26    4   10
# C4    4    3   19 1079   23    0
# C5    0    0    1   25  358    0
# C6    4    2    2    0    0  159

# Label  Accuracy Sensitivity Specificity Precision        F1
# C1           C1 0.9674663   0.9768212   0.9840608 0.9566275 0.9666189
# C2           C2 0.9674663   0.9857253   0.9946459 0.9864865 0.9861058
# C3           C3 0.9674663   0.9682937   0.9893048 0.9699122 0.9691023
# C4           C4 0.9674663   0.9309750   0.9938519 0.9565603 0.9435942
# C5           C5 0.9674663   0.9298701   0.9970265 0.9322917 0.9310793
# C6           C6 0.9674663   0.8833333   0.9991060 0.9520958 0.9164265
# Average Average 0.9674663   0.9458364   0.9929993 0.9589957 0.9521545
