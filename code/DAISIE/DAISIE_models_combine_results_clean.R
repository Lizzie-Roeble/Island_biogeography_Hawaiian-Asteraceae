# Required Libraries
library(readr)      
library(dplyr)      
library(writexl)    

# Note: The standard combine_DAISIE_ML_results_raw.R simply cats all fines. But the model results are different dimensions (ie type1 vs type2, and number of starting parameters)
# This script attempts to assign column headers and correctly match the columns for both type1 and type2, and then combine all the starting params into one column 
# The script runs, but it removes the models that failed to run

# Define directories
base_dir <- "/scratch/p300503/DAISIE_Hawaii/outputs/Models_Hawaii_constant"
ind_model <- file.path(base_dir, "Individual_models")  # Directory for individual model results
output_file <- file.path(base_dir, "DAISIE_ML_results_cleaned.xlsx")  # Final output file

# Set working directory
setwd(ind_model)

# Define column names for each type
type1_cols <- c("index", "dataset", "model", "init", "island_age", "M", "type", 
                "p_type2", "res", "ode", "replicate", "clado", "mu", "K", 
                "gamma", "ana", "loglik", "n_pars", "cond")

type2_cols <- c("index", "dataset", "model", "init", "island_age", "M", "type", 
                "p_type2", "res", "ode", "replicate", "clado", "mu", "K", 
                "gamma", "ana", "clado_2", "mu_2", "K_2", "gamma_2", "ana_2", 
                "type2_prop", "loglik", "n_pars", "cond")

# List result files
files <- list.files(pattern = "results_R[12]_.*\\.txt")

# Function to process files
process_file <- function(file) {
  # Read the file, handling errors
  data <- tryCatch(read.table(file, header = FALSE, fill = TRUE, stringsAsFactors = FALSE), 
                   error = function(e) return(NULL))
  
  if (is.null(data) || nrow(data) == 0) return(NULL)
  
  # Remove failed models - signaled by rows where NAs appear right after R1 or R2
  keep_rows <- apply(data, 1, function(row) {
    r_pos <- which(row %in% c("R1", "R2"))
    if (length(r_pos) == 0) return(FALSE)
    return(!is.na(row[r_pos + 1]))
  })
  
  data <- data[keep_rows, ]
  if (nrow(data) == 0) return(NULL)

  # Check model type
  model_type <- as.character(data[[7]])[1]  # Extract type from column 7
  
  # Process Type 1 models
  if (model_type == "type1") {
    if (ncol(data) < length(type1_cols)) return(NULL)  # Skip malformed files
    result <- data[, 1:length(type1_cols), drop = FALSE]
    colnames(result) <- type1_cols
    
    # Handle starting parameters
    if (ncol(data) > length(type1_cols)) {
      result$start_parameters <- apply(data[, (length(type1_cols) + 1):ncol(data), drop = FALSE], 1, paste, collapse = ",")
    } else {
      result$start_parameters <- NA
    }
    
    # Add empty columns for Type 2 parameters
    result <- result %>%
      mutate(clado_2 = NA, mu_2 = NA, K_2 = NA, gamma_2 = NA, ana_2 = NA, type2_prop = NA)
    
  } else {
    # Process Type 2 models
    if (ncol(data) < length(type2_cols)) return(NULL)  # Skip malformed files
    result <- data[, 1:length(type2_cols), drop = FALSE]
    colnames(result) <- type2_cols
    
    # Handle starting parameters
    if (ncol(data) > length(type2_cols)) {
      result$start_parameters <- apply(data[, (length(type2_cols) + 1):ncol(data), drop = FALSE], 1, paste, collapse = ",")
    } else {
      result$start_parameters <- NA
    }
  }
  
  return(result)
}

# Process all files
results <- lapply(files, process_file)
results <- Filter(Negate(is.null), results)  # Remove NULLs

# Combine data into a final dataframe
final_data <- bind_rows(results)

# Write output to an Excel file
write_xlsx(final_data, output_file)

