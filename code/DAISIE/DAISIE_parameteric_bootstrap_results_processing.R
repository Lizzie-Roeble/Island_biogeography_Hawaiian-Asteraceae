###############################################
# DAISIE 
# Parametric bootstrap: results - processing                           
###############################################

# Goal: Combine results from two boostrap runs for maxage and precise models
# run1 = 1000 simulations with no conditional filtering on >= 1 nonendemic
# run2 = filtered replicate (on >= 1 nonendemic) run to complete total 1000

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

# Read in the data

# Run1 MAXAGE
bootstrap_run1_maxage <- read.table("../outputs/DAISIE_parametric_bootstrap_results_maxage_m1_run1.txt", header = FALSE)
colnames(bootstrap_run1_maxage) <- bootstrap_colnames
# Run2 MAXAGE
bootstrap_run2_maxage <- read.table("../outputs/DAISIE_parametric_bootstrap_results_maxage_m1_run2.txt", header = FALSE)
colnames(bootstrap_run2_maxage) <- bootstrap_colnames

# Run1 PRECISE
bootstrap_run1_precise <- read.table("../outputs/DAISIE_parametric_bootstrap_results_precise_m1_run1.txt", header = FALSE)
colnames(bootstrap_run1_precise) <- bootstrap_colnames
# Run2 PRECISE
bootstrap_run2_precise <- read.table("../outputs/DAISIE_parametric_bootstrap_results_precise_m1_run2.txt", header = FALSE)
colnames(bootstrap_run2_precise) <- bootstrap_colnames

# STAC4 lists
stac4_replicates_maxage <- readLines("../outputs/DAISIE_maxage_m1_sims_filt_list.txt")
stac4_replicates_precise <- readLines("../outputs/DAISIE_precise_m1_sims_filt_list.txt")

# === PROCESS MAXAGE ===

# Filter Run1
bootstrap_run1_maxage_filtered <- bootstrap_run1_maxage %>%
  filter(Solver == "odeint::runge_kutta_cash_karp54") %>% # filter run 1 to karp54 solver
  group_by(Replicate) %>%
  slice_min(order_by = Likelihood, n = 1, with_ties = FALSE) %>% # for duplicate replicates, take one with best loglike
  ungroup() %>%
  filter(Replicate %in% stac4_replicates_maxage) # filter to only reps that meet condition (>= 1 non-endemic)

# Combine with extra sims
bootstrap_maxage_combined <- bind_rows(bootstrap_run1_maxage_filtered, bootstrap_run2_maxage)

write.table(bootstrap_maxage_combined,
            file = "../outputs/DAISIE_parametric_bootstrap_results_maxage_m1.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# === PROCESS PRECISE ===

# Filter Run1
bootstrap_run1_precise_filtered <- bootstrap_run1_precise %>%
  filter(Solver == "odeint::runge_kutta_cash_karp54") %>% # filter run 1 to karp54 solver
  group_by(Replicate) %>%
  slice_min(order_by = Likelihood, n = 1, with_ties = FALSE) %>% # for duplicate replicates, take one with best loglike
  ungroup() %>%
  filter(Replicate %in% stac4_replicates_precise) # filter to only reps that meet condition (>= 1 non-endemic)

# Combine with extra sims
bootstrap_precise_combined <- bind_rows(bootstrap_run1_precise_filtered, bootstrap_run2_precise)

write.table(bootstrap_precise_combined,
            file = "../outputs/DAISIE_parametric_bootstrap_results_precise_m1.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# Sanity checks

## how many rows did we get?
nrow(bootstrap_maxage_combined)
nrow(bootstrap_precise_combined)

# check for any accidental duplicate Replicate IDs (should be none)
any(duplicated(bootstrap_maxage_combined$Replicate))
any(duplicated(bootstrap_precise_combined$Replicate))
