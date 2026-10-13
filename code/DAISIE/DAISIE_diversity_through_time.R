###############################
# DAISIE 
# Diversity Through Time (DTT)
# DAISIE 5
###############################

# Load required libraries
library(DAISIE)

# Load simulations: 1000 simulations filtered on condition of >= 1 nonendemic
# maxage - 5.14myr
load("../data/DAISIE5_sims/DAISIE_maxage_5.14myr_sims_filt_1000.Rdata")
island_replicates_maxage_5.14myr <- DAISIE_maxage_5.14myr_m1_sims_filt_1000
# precise - 5.14myr
load("../data/DAISIE5_sims/DAISIE_precise_5.14myr_sims_filt_1000.Rdata")
island_replicates_precise_5.14myr <- DAISIE_precise_5.14myr_m1_sims_filt_1000
# maxage - 6.15myr
load("../data/DAISIE5_sims/DAISIE_maxage_6.15myr_sims_filt_1000.Rdata")
island_replicates_maxage_6.15myr <- DAISIE_maxage_6.15myr_m1_sims_filt_1000
# precise - 6.15myr
load("../data/DAISIE5_sims/DAISIE_precise_6.15myr_sims_filt_1000.Rdata")
island_replicates_precise_6.15myr <- DAISIE_precise_6.15myr_m1_sims_filt_1000


# Plot DTT 
# maxage - 5.14myr
pdf(file = "../outputs/DAISIE_dtt_maxage_5.14myr_DAISIE5_raw.pdf", width = 8, height = 8.5)  
DAISIE_plot_sims(island_replicates_maxage_5.14myr)
dev.off()
# precise - 5.14myr
pdf(file = "../outputs/DAISIE_dtt_precise_5.14myr_DAISIE5_raw.pdf", width = 8, height = 8.5)  
DAISIE_plot_sims(island_replicates_precise_5.14myr)
dev.off()
# maxage - 6.15myr
pdf(file = "../outputs/DAISIE_dtt_maxage_6.15myr_DAISIE5_raw.pdf", width = 8, height = 8.5)  
DAISIE_plot_sims(island_replicates_maxage_6.15myr)
dev.off()
# precise - 6.15myr
pdf(file = "../outputs/DAISIE_dtt_precise_6.15myr_DAISIE5_raw.pdf", width = 8, height = 8.5)  
DAISIE_plot_sims(island_replicates_precise_6.15myr)
dev.off()
