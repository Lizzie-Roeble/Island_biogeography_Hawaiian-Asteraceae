# Goal: get a list of the simulations that have at least one non-endemic
# DAISIE non-endemic = stac 4

# Load the 1000 replicate dataset
load("../outputs/DAISIE_Hawaii_sims_precise_RE.Rdata")
island_replicates_1000 <- sim_results

# Load the 10000 replicate dataset
load("../outputs/DAISIE_Hawaii_sims_precise_m1_10000.Rdata")
island_replicates_10000 <- sim_results

# Check if first 1000 of new dataset match old dataset
identical(island_replicates_1000, island_replicates_10000[1:1000])

# Function to check if a replicate has any STAC 4 species
has_stac4 <- function(replicate_data) {
  species_data <- replicate_data[-1]
  any(sapply(species_data, function(x) x$stac == 4))
}

# Count replicates by STAC value (1, 2, 3, 4)
for (stac_val in 1:4) {
  count <- sum(sapply(island_replicates_10000, function(rep) {
    species_data <- rep[-1]
    any(sapply(species_data, function(x) x$stac == stac_val))
  }))
  print(paste("STAC", stac_val, ":", count))
}

# Identify replicates with at least one non-endemic species (STAC 4)
precise_m1_sims_filt <- sapply(island_replicates_10000, has_stac4)

# Count replicates with STAC 4
sum(precise_m1_sims_filt)   

# Get replicate numbers with STAC 4
precise_m1_sims_filt <- which(precise_m1_sims_filt)

# Save 
write.table(precise_m1_sims_filt,
            file = "../outputs/Hawaii_precise_m1_sims_filtered_stac4.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
