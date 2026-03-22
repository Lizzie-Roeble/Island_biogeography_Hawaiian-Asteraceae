###############################
# DAISIE 
# Diversity Through Time (DTT)                           
###############################

# Load required libraries
library(DAISIE)

# Load simulations: 1000 simulations filtered on condition of >= 1 nonendemic
load("../outputs/DAISIE_sims_precise_m1_filt_1000.Rdata")
load("../outputs/DAISIE_sims_maxage_m1_filt_1000.Rdata")
island_replicates_maxage_m1 <- DAISIE_maxage_m1_sims_filt_1000
island_replicates_precise_m1 <- DAISIE_precise_m1_sims_filt_1000

# Plot DTT for maxage
pdf(file = "../outputs/DAISIE_dtt_maxage_m1_raw.pdf", width = 8, height = 8.5)  
DAISIE_plot_sims(DAISIE_maxage_m1_sims_filt_1000)
dev.off()

# Plot DTT for precise age
pdf(file = "../outputs/DAISIE_dtt_precise_m1_raw.pdf", width = 8, height = 8.5)  
DAISIE_plot_sims(DAISIE_precise_m1_sims_filt_1000)
dev.off()

