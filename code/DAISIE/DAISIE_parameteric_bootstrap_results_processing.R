###############################################
# DAISIE 
# Parametric bootstrap: results - processing                           
# Updated with DAISIE5
###############################################

# Goal: Process results from the bootstrap runs for maxage and precise models (5.14myr and 6.15myr)
# Raw results only contain replicates with >= 1 nonendemic; keep the first 1000 (by replicate number) per model

# Load libraries
library(tidyverse)

# Column names
bootstrap_colnames <- c(
  "Replicate", "Island_age", "Mainland_pool", "Type", "Res", "Model", 
  "Solver", "Cladogenesis", "Extinction", "Carrying_Capacity", 
  "Colonization", "Anagenesis", "Likelihood", "Number_of_parameters", 
  "Number_parameters_type2", "Cladogenesis_MLE", "Extinction_MLE", 
  "Carrying_Capacity_MLE", "Colonization_MLE", "Anagenesis_MLE"
)

# Function to read bootstrap results, keeping only complete lines (exactly 20 fields)
read_bootstrap <- function(file) {
  lines <- readLines(file)
  n_fields <- sapply(strsplit(trimws(lines), "[[:space:]]+"), length)
  cat(file, "- removed", sum(n_fields != 20), "malformed line(s)\n")
  dat <- read.table(text = lines[n_fields == 20], header = FALSE)
  colnames(dat) <- bootstrap_colnames
  dat
}

# Read in the data

# MAXAGE 5.14myr
bootstrap_maxage_5.14myr <- read_bootstrap("../data/DAISIE5_sims/DAISIE_Hawaii_sim_maxage_5.14myr.txt")
# MAXAGE 6.15myr
bootstrap_maxage_6.15myr <- read_bootstrap("../data/DAISIE5_sims/DAISIE_Hawaii_sim_maxage_6.15myr.txt")

# PRECISE 5.14myr
bootstrap_precise_5.14myr <- read_bootstrap("../data/DAISIE5_sims/DAISIE_Hawaii_sim_precise_5.14myr.txt")
# PRECISE 6.15myr
bootstrap_precise_6.15myr <- read_bootstrap("../data/DAISIE5_sims/DAISIE_Hawaii_sim_precise_6.15myr.txt")

# === PROCESS MAXAGE 5.14myr ===

bootstrap_maxage_5.14myr_filtered <- bootstrap_maxage_5.14myr %>%
  arrange(Replicate) %>% # order by replicate number
  slice_head(n = 1000) # keep first 1000 replicates with an ML result

write.table(bootstrap_maxage_5.14myr_filtered,
            file = "../data/DAISIE5_sims/DAISIE5_parametric_bootstrap_results_maxage_5.14myr.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# === PROCESS MAXAGE 6.15myr ===

bootstrap_maxage_6.15myr_filtered <- bootstrap_maxage_6.15myr %>%
  arrange(Replicate) %>% # order by replicate number
  slice_head(n = 1000) # keep first 1000 replicates with an ML result

write.table(bootstrap_maxage_6.15myr_filtered,
            file = "../data/DAISIE5_sims/DAISIE5_parametric_bootstrap_results_maxage_6.15myr.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# === PROCESS PRECISE 5.14myr ===

bootstrap_precise_5.14myr_filtered <- bootstrap_precise_5.14myr %>%
  arrange(Replicate) %>% # order by replicate number
  slice_head(n = 1000) # keep first 1000 replicates with an ML result

write.table(bootstrap_precise_5.14myr_filtered,
            file = "../data/DAISIE5_sims/DAISIE5_parametric_bootstrap_results_precise_5.14myr.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# === PROCESS PRECISE 6.15myr ===

bootstrap_precise_6.15myr_filtered <- bootstrap_precise_6.15myr %>%
  arrange(Replicate) %>% # order by replicate number
  slice_head(n = 1000) # keep first 1000 replicates with an ML result

write.table(bootstrap_precise_6.15myr_filtered,
            file = "../data/DAISIE5_sims/DAISIE5_parametric_bootstrap_results_precise_6.15myr.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# Sanity checks

## how many rows did we get? (should be 1000 each)
nrow(bootstrap_maxage_5.14myr_filtered)
nrow(bootstrap_maxage_6.15myr_filtered)
nrow(bootstrap_precise_5.14myr_filtered)
nrow(bootstrap_precise_6.15myr_filtered)

# check for any accidental duplicate Replicate IDs (should be none)
any(duplicated(bootstrap_maxage_5.14myr_filtered$Replicate))
any(duplicated(bootstrap_maxage_6.15myr_filtered$Replicate))
any(duplicated(bootstrap_precise_5.14myr_filtered$Replicate))
any(duplicated(bootstrap_precise_6.15myr_filtered$Replicate))
