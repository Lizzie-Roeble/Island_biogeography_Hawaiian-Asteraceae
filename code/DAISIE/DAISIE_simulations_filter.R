#########################################
# DAISIE 
# Model Simulations: conditional filter                            
#########################################

# Aim: Load 10,000 simulations (maxage and precise) and filter simulations on condition of >= 1 nonendemic
# Input: 10,000 simulations for the maxage and precise age DAISIE models (m1: DD all parameters) (DIAISE_precise_m1_sims_10000; DIAISE_maxage_m1_sims_10000)
# Output1: Simulations filtered to condition (DIAISE_precise_m1_sims_filt; DIAISE_maxage_m1_sims_filt)
# Output2: List of simulations that meet the condition (DAISIE_precise_m1_sims_filt_list.txt; DAISIE_maxage_m1_sims_filt_list.txt)
# Output3: First 1,000 simulations filtered to condition (DIAISE_precise_m1_sims_filt_1000; DIAISE_maxage_m1_sims_filt_1000)

# Load required libraries
# devtools::install_github("rsetienne/DAISIE", ref = "develop")
library(readxl)
library(tidyverse)
library(DAISIE)

# Load the 10,000 simulations (run by RE)
load("../outputs/DAISIE_sims_maxage_m1_10000.Rdata")
island_replicates_maxage_10000 <- island_replicates
load("../outputs/DAISIE_sims_precise_m1_10000.Rdata")
island_replicates_precise_10000 <- sim_results

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

# 1. Count replicates by STAC value (1, 2, 3, 4) for both datasets  
cat("==== STAC counts for PRECISE ====\n")  
count_stac(island_replicates_precise_10000, "PRECISE")  

cat("==== STAC counts for MAXAGE ====\n")  
count_stac(island_replicates_maxage_10000, "MAXAGE")  

# 2. Identify replicates with at least one non-endemic species (STAC 4)  
precise_m1_sims_stac4_filter <- sapply(island_replicates_precise_10000, has_stac4)  
maxage_m1_sims_stac4_filter <- sapply(island_replicates_maxage_10000, has_stac4)  

# 3. Count replicates with STAC 4  
precise_count_stac4 <- sum(precise_m1_sims_stac4_filter)  
maxage_count_stac4 <- sum(maxage_m1_sims_stac4_filter) 

# 4. Get a list of replicate numbers with STAC 4 for both datasets  
precise_stac4_replicate_numbers <- which(precise_m1_sims_stac4_filter)  
maxage_stac4_replicate_numbers <- which(maxage_m1_sims_stac4_filter)

# Save replicate lists  
write.table(precise_stac4_replicate_numbers,  
            file = "../outputs/DAISIE_precise_m1_sims_filt_list.txt",  
            quote = FALSE,  
            row.names = FALSE,  
            col.names = FALSE)  
write.table(maxage_stac4_replicate_numbers,  
            file = "../outputs/DAISIE_maxage_m1_sims_filt_list.txt",  
            quote = FALSE,  
            row.names = FALSE,  
            col.names = FALSE)  

# Subset full dataset to only the replicates with at least one STAC 4
DAISIE_precise_m1_sims_filt <- island_replicates_precise_10000[precise_stac4_replicate_numbers]
DAISIE_maxage_m1_sims_filt <- island_replicates_maxage_10000[maxage_stac4_replicate_numbers]

# Save as .Rdata files
save(DAISIE_precise_m1_sims_filt,
     file = "../outputs/DAISIE_precise_m1_sims_filt.Rdata")
save(DAISIE_maxage_m1_sims_filt,
     file = "../outputs/DAISIE_maxage_m1_sims_filt.Rdata")

# 5. Get the first 1000 replicates with STAC 4 for both datasets  
precise_stac4_replicate_numbers_1000 <- precise_stac4_replicate_numbers[1:min(1000, length(precise_stac4_replicate_numbers))]  
maxage_stac4_replicate_numbers_1000 <- maxage_stac4_replicate_numbers[1:min(1000, length(maxage_stac4_replicate_numbers))]  

# Subset the original list to only the selected replicates  
DAISIE_precise_m1_sims_filt_1000 <- island_replicates_precise_10000[precise_stac4_replicate_numbers_1000]  
DAISIE_maxage_m1_sims_filt_1000 <- island_replicates_maxage_10000[maxage_stac4_replicate_numbers_1000]  

# Save as .Rdata files  
save(DAISIE_precise_m1_sims_filt_1000,  
     file = "../outputs/DAISIE_sims_precise_m1_filt_1000.Rdata")  
save(DAISIE_maxage_m1_sims_filt_1000,  
     file = "../outputs/DAISIE_sims_maxage_m1_filt_1000.Rdata") 


