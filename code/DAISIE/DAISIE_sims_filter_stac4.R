# Goal: get a list of the simulations that have atleast one non-endemic

load("../outputs/Hawaii_maxage_m1_sims.RData")
load("../outputs/Hawaii_precise_m1_sims.RData")

# Function to check if a replicate has any STAC 4 species
has_stac4 <- function(replicate_data) {
  species_data <- replicate_data[-1]
  any(sapply(species_data, function(x) x$stac == 4))
}

# Check both datasets
precise_m1_sim_filt <- sapply(Hawaii_precise_m1_sims, has_stac4)
maxage_m1_sims_filt <- sapply(Hawaii_maxage_m1_sims, has_stac4)

# Count replicates with STAC 4
sum(precise_m1_sim_filt)  # Hawaii_precise_m1_sims
sum(maxage_m1_sims_filt)   # Hawaii_maxage_m1_sims

# Get replicate numbers with STAC 4
precise_m1_sim_filt <- which(precise_m1_sim_filt)
maxage_m1_sims_filt <- which(maxage_m1_sims_filt)

# Save 
write.table(maxage_m1_sims_filt,
            file = "../outputs/Hawaii_maxage_m1_sims_filtered_stac4.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
write.table(precise_m1_sim_filt,
            file = "../outputs/Hawaii_precise_m1_sims_filtered_stac4.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)