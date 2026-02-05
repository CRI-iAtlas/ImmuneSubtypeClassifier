# renv::install('/users/dgibbs/Code/robencla/')
# renv::snapshot()
# renv::activate()

# renv::install('.')
#
# library(robencla)
# #
# # # 1. Load your training data
# # # (Replace with your actual data loading line)
# # # data <- read.csv("path/to/your/training_data.csv")
# #
# # # # 2. Re-train and re-save the model
# # # # This forces the model to be built using the NEW, FIXED class definition
library(ImmuneSubtypeClassifier)

build_robencla_classifier(
   data_path = '../data/PanCancer_TCGA/EBpp_pancancer.csv.gz', # Update this path
   output_path = "model/robencla_trained_model.rds",
   label_name = "Label",      # Verify this matches your data
   sample_id = "Barcode",      # Verify this matches your data
   pair_list = NULL                  # Uses default pairs
)



ebpp <- readr::read_csv('../data/PanCancer_TCGA/EBpp_pancancer.csv.gz')

# model_loaded <- readRDS('model/robencla_trained_model.rds')

# Force the input to a plain data frame and provide a dummy barcode
ebpp_df <- data.frame(ebpp[1:10,], check.names = FALSE)
ebpp_df$Barcode[1:10]
# has ID

ImmuneSubtypeClassifier::callSubtypes(ebpp_df, geneid = "symbol", sampleid = "Barcode")
