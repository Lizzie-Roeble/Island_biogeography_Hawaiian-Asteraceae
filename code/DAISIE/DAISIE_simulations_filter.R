#########################################
# DAISIE 
# Model Simulations: conditional filter 
# Updated with DAISIE5
#########################################

# Aim: Load 10,000 simulations (maxage and precise, 5.14myr and 5.15myr) and filter simulations on condition of >= 1 nonendemic

# Load required libraries
# devtools::install_github("rsetienne/DAISIE", ref = "develop")
library(readxl)
library(tidyverse)
library(DAISIE)

# Load the 10,000 simulations (run by RE)

# Max age - 5.14myr
load("../data/DAISIE5_sims/DAISIE_Hawaii_sim_maxage_5.14myr.Rdata")
island_replicates_maxage_5.14myr_10000 <- sim_results
rm(sim_results)

# Precise age - 5.14myr
load("../data/DAISIE5_sims/DAISIE_Hawaii_sim_precise_5.14myr.Rdata")
island_replicates_precise_5.14myr_10000 <- sim_results
rm(sim_results)

# Max age - 6.15myr
load("../data/DAISIE5_sims/DAISIE_Hawaii_sim_maxage_6.15myr.Rdata")
island_replicates_maxage_6.15myr_10000 <- sim_results
rm(sim_results)

# Precise age - 6.15myr
load("../data/DAISIE5_sims/DAISIE_Hawaii_sim_precise_6.15myr.Rdata")
island_replicates_precise_6.15myr_10000 <- sim_results
rm(sim_results)

# Function to check if a replicate has any STAC 4 species
has_stac4 <- function(replicate_data) {
  species_data <- replicate_data[-1]
  any(sapply(species_data, function(x) x$stac == 4))
}

# Function to count replicates by STAC value
count_stac <- function(dataset, label) {
  for (stac_val in 1:4) {
    count <- sum(sapply(dataset, function(rep) {
      species_data <- rep[-1]
      any(sapply(species_data, function(x) x$stac == stac_val))
    }))
    print(paste(label, "- STAC", stac_val, ":", count))
  }
}

# 1. Count replicates by STAC value (1, 2, 3, 4) for all datasets
cat("==== STAC counts for MAXAGE 5.14myr ====\n")
count_stac(island_replicates_maxage_5.14myr_10000, "MAXAGE")

cat("==== STAC counts for PRECISE 5.14myr ====\n")
count_stac(island_replicates_precise_5.14myr_10000, "PRECISE")

cat("==== STAC counts for MAXAGE 6.15myr ====\n")
count_stac(island_replicates_maxage_6.15myr_10000, "MAXAGE")

cat("==== STAC counts for PRECISE 6.15myr ====\n")
count_stac(island_replicates_precise_6.15myr_10000, "PRECISE")

# 2. Identify replicates with at least one non-endemic species (STAC 4)
precise_5.14myr_m1_sims_stac4_filter <- sapply(island_replicates_precise_5.14myr_10000, has_stac4)
maxage_5.14myr_m1_sims_stac4_filter  <- sapply(island_replicates_maxage_5.14myr_10000, has_stac4)
precise_6.15myr_m1_sims_stac4_filter <- sapply(island_replicates_precise_6.15myr_10000, has_stac4)
maxage_6.15myr_m1_sims_stac4_filter  <- sapply(island_replicates_maxage_6.15myr_10000, has_stac4)

# 3. Count replicates with STAC 4
precise_5.14myr_count_stac4 <- sum(precise_5.14myr_m1_sims_stac4_filter)
maxage_5.14myr_count_stac4  <- sum(maxage_5.14myr_m1_sims_stac4_filter)
precise_6.15myr_count_stac4 <- sum(precise_6.15myr_m1_sims_stac4_filter)
maxage_6.15myr_count_stac4  <- sum(maxage_6.15myr_m1_sims_stac4_filter)

# 4. Get a list of replicate numbers with STAC 4 for all datasets
precise_5.14myr_stac4_replicate_numbers <- which(precise_5.14myr_m1_sims_stac4_filter)
maxage_5.14myr_stac4_replicate_numbers  <- which(maxage_5.14myr_m1_sims_stac4_filter)
precise_6.15myr_stac4_replicate_numbers <- which(precise_6.15myr_m1_sims_stac4_filter)
maxage_6.15myr_stac4_replicate_numbers  <- which(maxage_6.15myr_m1_sims_stac4_filter)

# Save replicate lists (all replicates with >= 1 STAC 4)
write.table(precise_5.14myr_stac4_replicate_numbers,
            file = "../data/DAISIE5_sims/DAISIE_precise_5.14myr_m1_sims_filt_list.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

write.table(maxage_5.14myr_stac4_replicate_numbers,
            file = "../data/DAISIE5_sims/DAISIE_maxage_5.14myr_m1_sims_filt_list.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

write.table(precise_6.15myr_stac4_replicate_numbers,
            file = "../data/DAISIE5_sims/DAISIE_precise_6.15myr_m1_sims_filt_list.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

write.table(maxage_6.15myr_stac4_replicate_numbers,
            file = "../data/DAISIE5_sims/DAISIE_maxage_6.15myr_m1_sims_filt_list.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# Subset full dataset to only the replicates with at least one STAC 4
DAISIE_precise_5.14myr_m1_sims_filt <- island_replicates_precise_5.14myr_10000[precise_5.14myr_stac4_replicate_numbers]
DAISIE_maxage_5.14myr_m1_sims_filt  <- island_replicates_maxage_5.14myr_10000[maxage_5.14myr_stac4_replicate_numbers]
DAISIE_precise_6.15myr_m1_sims_filt <- island_replicates_precise_6.15myr_10000[precise_6.15myr_stac4_replicate_numbers]
DAISIE_maxage_6.15myr_m1_sims_filt  <- island_replicates_maxage_6.15myr_10000[maxage_6.15myr_stac4_replicate_numbers]

# Save as .Rdata files
save(DAISIE_precise_5.14myr_m1_sims_filt,
     file = "../data/DAISIE5_sims/DAISIE_precise_5.14myr_sims_filt.Rdata")
save(DAISIE_maxage_5.14myr_m1_sims_filt,
     file = "../data/DAISIE5_sims/DAISIE_maxage_5.14myr_sims_filt.Rdata")
save(DAISIE_precise_6.15myr_m1_sims_filt,
     file = "../data/DAISIE5_sims/DAISIE_precise_6.15myr_sims_filt.Rdata")
save(DAISIE_maxage_6.15myr_m1_sims_filt,
     file = "../data/DAISIE5_sims/DAISIE_maxage_6.15myr_sims_filt.Rdata")

# 5. Get the first 1000 replicates with STAC 4 for all datasets
precise_5.14myr_stac4_replicate_numbers_1000 <- precise_5.14myr_stac4_replicate_numbers[1:min(1000, length(precise_5.14myr_stac4_replicate_numbers))]
maxage_5.14myr_stac4_replicate_numbers_1000  <- maxage_5.14myr_stac4_replicate_numbers[1:min(1000, length(maxage_5.14myr_stac4_replicate_numbers))]
precise_6.15myr_stac4_replicate_numbers_1000 <- precise_6.15myr_stac4_replicate_numbers[1:min(1000, length(precise_6.15myr_stac4_replicate_numbers))]
maxage_6.15myr_stac4_replicate_numbers_1000  <- maxage_6.15myr_stac4_replicate_numbers[1:min(1000, length(maxage_6.15myr_stac4_replicate_numbers))]

# Save replicate lists (first 1000 replicates with >= 1 STAC 4)
write.table(precise_5.14myr_stac4_replicate_numbers_1000,
            file = "../data/DAISIE5_sims/DAISIE_precise_5.14myr_m1_sims_filt_list_1000.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

write.table(maxage_5.14myr_stac4_replicate_numbers_1000,
            file = "../data/DAISIE5_sims/DAISIE_maxage_5.14myr_m1_sims_filt_list_1000.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

write.table(precise_6.15myr_stac4_replicate_numbers_1000,
            file = "../data/DAISIE5_sims/DAISIE_precise_6.15myr_m1_sims_filt_list_1000.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

write.table(maxage_6.15myr_stac4_replicate_numbers_1000,
            file = "../data/DAISIE5_sims/DAISIE_maxage_6.15myr_m1_sims_filt_list_1000.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# Subset the original list to only the selected replicates
DAISIE_precise_5.14myr_m1_sims_filt_1000 <- island_replicates_precise_5.14myr_10000[precise_5.14myr_stac4_replicate_numbers_1000]
DAISIE_maxage_5.14myr_m1_sims_filt_1000  <- island_replicates_maxage_5.14myr_10000[maxage_5.14myr_stac4_replicate_numbers_1000]
DAISIE_precise_6.15myr_m1_sims_filt_1000 <- island_replicates_precise_6.15myr_10000[precise_6.15myr_stac4_replicate_numbers_1000]
DAISIE_maxage_6.15myr_m1_sims_filt_1000  <- island_replicates_maxage_6.15myr_10000[maxage_6.15myr_stac4_replicate_numbers_1000]

# Save as .Rdata files
save(DAISIE_precise_5.14myr_m1_sims_filt_1000,
     file = "../data/DAISIE5_sims/DAISIE_precise_5.14myr_sims_filt_1000.Rdata")
save(DAISIE_maxage_5.14myr_m1_sims_filt_1000,
     file = "../data/DAISIE5_sims/DAISIE_maxage_5.14myr_sims_filt_1000.Rdata")
save(DAISIE_precise_6.15myr_m1_sims_filt_1000,
     file = "../data/DAISIE5_sims/DAISIE_precise_6.15myr_sims_filt_1000.Rdata")
save(DAISIE_maxage_6.15myr_m1_sims_filt_1000,
     file = "../data/DAISIE5_sims/DAISIE_maxage_6.15myr_sims_filt_1000.Rdata")
