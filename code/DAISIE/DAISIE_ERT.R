################################### 
# DAISIE
# Evolutionary Return Time (ERT) 
# DAISIE 5
###################################

# Install development branch of DAISIE 
# devtools::install_github("rsetienne/DAISIE", ref = "develop")

# Load required libraries
library(DAISIE)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(patchwork)
library(writexl)

# Clear workspace
rm(list = ls())

#-----------------------------------
# 1. FUNCTIONS
#-----------------------------------

# Calculate specific quantiles from a probability distribution
#
# @param probdist A 2-row matrix with values in first row and probabilities in second
# @param probs Vector of probabilities for which to calculate quantiles
# @return Named vector of quantiles
calculate_quantiles <- function(probdist, probs) {
  result <- NULL
  cdf <- cumsum(probdist[2, ])
  
  for (i in 1:length(probs)) {
    n <- max(which(cdf <= probs[i]))
    x <- probdist[1, n]
    
    if (cdf[n] == probs[i]) {
      result[i] <- x
    } else if (n < length(cdf)) {
      result[i] <- ((x + 1) * (probs[i] - cdf[n]) + x * (cdf[n + 1] - probs[i])) / (cdf[n + 1] - cdf[n])
    } else {
      result[i] <- x
    }
  }
  
  names(result) <- probs
  return(result)
}

# Calculate the recovery time by finding when diversity reaches target level
#
# @param target_diversity The diversity level to recover to
# @param diversity_timeseries Vector of diversity values over time
# @param timestep Time increment between diversity values
# @return Time when diversity reaches target
calculate_recovery_time <- function(target_diversity, diversity_timeseries, timestep) {
  # Find first point where diversity exceeds target
  for (i in 2:(length(diversity_timeseries) / timestep)) {
    if (diversity_timeseries[i] > target_diversity) {
      # Using linear interpolation between timesteps: y = ax + b -> x = (y - b) / a
      y2 <- diversity_timeseries[i]
      y1 <- diversity_timeseries[i - 1]
      dy <- y2 - y1
      x2 <- (i - 1) * timestep
      x1 <- (i - 2) * timestep
      dx <- x2 - x1
      a <- dy / dx
      b <- y2 - a * x2
      x <- (target_diversity - b) / a
      return(x)
    }
  }
  return(NA)  # Return NA if target not reached
}

# Calculate quantiles for species distributions over time
#
# @param marginal_prob_dist Result from DAISIE_margprobdist
# @param species_index Index for species type (1=endemic, 3=total)
# @param max_species Maximum number of species to consider
# @param probability_levels Vector of probability levels for quantiles
# @return Matrix of quantiles with time as first column
calculate_species_quantiles <- function(marginal_prob_dist, species_index, max_species, probability_levels) {
  num_species <- 0:(max_species - 1)
  result_matrix <- matrix(ncol = length(probability_levels) + 1, nrow = nrow(marginal_prob_dist[[species_index]]))
  result_matrix[, 1] <- 0:(nrow(marginal_prob_dist[[species_index]]) - 1)
  
  for (i in 1:nrow(marginal_prob_dist[[species_index]])) {
    prob_dist <- rbind(num_species, marginal_prob_dist[[species_index]][i, ])
    result_matrix[i, 2:(length(probability_levels) + 1)] <- calculate_quantiles(prob_dist, probability_levels)
  }
  
  return(result_matrix)
}

#-----------------------------------
# 2. DEFINE MODEL PARAMETERS
#-----------------------------------

# Model parameters
# Maxage - 6.15myr
pars_maxage <- c(1.92637708766939, 0.966436759805301, 74.1983894314959, 0.000509148912158619, 9.55293550834906)
# Precise - 6.15myr
pars_precise <- c(1.578244508, 0.79944063, 100.5920454, 0.000498832, 1.468723953)


# Simulation constants
mainland_pool_size <- 7000
time_step <- 0.1 
simulation_time <- 7  # Maximum time in Myr
time_points <- seq(0.000001, simulation_time, time_step)
max_species <- 300

#-----------------------------------
# 3. DEFINE DIVERSITY SCENARIOS
#-----------------------------------

# Hawaii-specific diversity values
total_species_count <- 103  # 101 described species + 2 independent P. sandwicenisum taxa
extinct_species_count <- 12
non_endemic_species_count <- 1

# Generalized diversity values
prehuman_diversity <- total_species_count # pre-human = pre-european arrival (1778) 
current_diversity_total <- total_species_count - extinct_species_count
current_diversity_endemic <- total_species_count - extinct_species_count - non_endemic_species_count
current_diversity_non_endemic <- non_endemic_species_count

# Future IUCN diversity projections
iucn_threatened_count <- 53  # All endemic
future_iucn_total <- current_diversity_total - iucn_threatened_count
future_iucn_endemic <- current_diversity_endemic - iucn_threatened_count
future_iucn_non_endemic <- non_endemic_species_count  

# Lineage-specific data for Diverity-Dependent models (m1)
# Current diversity by lineage (endemic, non-endemic)
current_diversity_by_lineage <- rbind(
  c(31, 0),  # silverswords 34 - 3 = 31 (extinct)
  c(20, 0),  # bidens 20 - 0
  c(16, 0),  # lipochaeta 20 - 4
  c(8, 0),   # tetramolopium 12 - 4
  c(4, 0),   # hesperomannia 4 - 0
  c(3, 0),   # artemisia 3 - 0
  c(3, 0),   # keysseria 3 - 0
  c(2, 0),   # remya 3 - 1
  c(1, 0),   # pseudo sandwicensium a 1 - 0
  c(1, 0),   # pseudo sandwicensium b 1 - 0
  c(1, 0),   # pseudo sandwicensium c 1 - 0
  c(0, 1)    # adenostemma 1 - 0
)

# Future IUCN-threatened diversity by lineage (endemic, non-endemic)
future_diversity_by_lineage <- rbind(
  c(10, 0),  # silverswords 34 - 3 (extinct) - 21 (threatened)
  c(6, 0),   # bidens 20 - 0 - 14
  c(10, 0),  # lipochaeta 20 - 4 - 6
  c(4, 0),   # tetramolopium 12 - 4 - 4
  c(0, 0),   # hesperomannia 4 - 0 - 4
  c(2, 0),   # artemisia 3 - 0 - 1
  c(1, 0),   # keysseria 3 - 0 - 2
  c(1, 0),   # remya 3 - 1 - 1
  c(1, 0),   # pseudo sandwicensium a 1 - 0 - 0
  c(1, 0),   # pseudo sandwicensium b 1 - 0 - 0
  c(1, 0),   # pseudo sandwicensium c 1 - 0 - 0
  c(0, 1)    # adenostemma 1 - 0 - 0
)

#-----------------------------------
# 4. ERT ANALYSIS: Max Age 6.15myr
#-----------------------------------

#---
# 4.1 Current to Pre-human for Max Age M1
#---

# First, convert the lineage-specific data to total counts for DAISIE_margprobdist
current_to_prehuman_start_maxage <- c(sum(current_diversity_by_lineage[,1]), 
                                      sum(current_diversity_by_lineage[,2]))

# Calculate expected species numbers over time using ExpEIN2 for M1
current_to_prehuman_exp_maxage <- DAISIE_ExpEIN2(
  tvec = time_points,
  pars = pars_maxage,
  M = mainland_pool_size,
  initEI = current_diversity_by_lineage
)

# Calculate marginal probability distributions using DAISIE_margprobdist2
current_to_prehuman_prob_dist_maxage <- DAISIE:::DAISIE_margprobdist2(
  tvec = time_points,
  pars = pars_maxage,
  M = mainland_pool_size,
  initEI = current_diversity_by_lineage,
  res = max_species,
  ddmodel = 11,
  methode = "ode45",
  reltolint = 1e-16,
  abstolint = 1e-16
)

# Set probability levels for confidence intervals
probability_levels <- c(0.05, 0.1, 0.5, 0.9, 0.95)

# Calculate quantiles for total species
total_species_quantiles_c2p_maxage <- calculate_species_quantiles(
  current_to_prehuman_prob_dist_maxage,
  species_index = 3,  # Total species
  max_species = max_species,
  probability_levels = probability_levels
)

# Calculate recovery times
recovery_time_extinct_maxage <- calculate_recovery_time(
  prehuman_diversity,
  current_to_prehuman_exp_maxage[[3]],
  time_step
)

recovery_time_extinct_95_maxage <- calculate_recovery_time(
  prehuman_diversity * 0.95,
  current_to_prehuman_exp_maxage[[3]],
  time_step
)

#---
# 4.2 Future IUCN to Current for Max Age M1
#---

# Convert the lineage-specific data to total counts for DAISIE_margprobdist
future_to_current_start_maxage <- c(sum(future_diversity_by_lineage[,1]), 
                                    sum(future_diversity_by_lineage[,2]))

# Calculate expected species numbers over time using ExpEIN2 for M1
future_to_current_exp_maxage <- DAISIE_ExpEIN2(
  tvec = time_points,
  pars = pars_maxage,
  M = mainland_pool_size,
  initEI = future_diversity_by_lineage
)

future_to_current_prob_dist_maxage <- DAISIE:::DAISIE_margprobdist2(
  tvec = time_points,
  pars = pars_maxage,
  M = mainland_pool_size,
  initEI = future_diversity_by_lineage,
  res = max_species,
  ddmodel = 11,
  methode = "ode45",
  reltolint = 1e-16,
  abstolint = 1e-16
)

# Calculate quantiles for total species
total_species_quantiles_f2c_maxage <- calculate_species_quantiles(
  future_to_current_prob_dist_maxage,
  species_index = 3,  # Total species
  max_species = max_species,
  probability_levels = probability_levels
)

# Calculate recovery times
recovery_time_future_maxage <- calculate_recovery_time(
  current_diversity_total,
  future_to_current_exp_maxage[[3]],
  time_step
)

recovery_time_future_95_maxage <- calculate_recovery_time(
  current_diversity_total * 0.95,
  future_to_current_exp_maxage[[3]],
  time_step
)

recovery_time_future_prehuman_maxage <- calculate_recovery_time(
  prehuman_diversity,
  future_to_current_exp_maxage[[3]],
  time_step
)

recovery_time_future_prehuman_95_maxage <- calculate_recovery_time(
  prehuman_diversity * 0.95,
  future_to_current_exp_maxage[[3]],
  time_step
)

#-----------------------------------
# 5. ERT ANALYSIS: Precise 6.15myr (DD)
#-----------------------------------

#---
# 5.1 Current to Pre-human for Precise M1
#---

# First, convert the lineage-specific data to total counts for DAISIE_margprobdist
current_to_prehuman_start_precise <- c(sum(current_diversity_by_lineage[,1]), 
                                       sum(current_diversity_by_lineage[,2]))

# Calculate expected species numbers over time using ExpEIN2 for M1
current_to_prehuman_exp_precise <- DAISIE_ExpEIN2(
  tvec = time_points,
  pars = pars_precise,
  M = mainland_pool_size,
  initEI = current_diversity_by_lineage
)

# Calculate marginal probability distributions using DAISIE_margprobdist2
current_to_prehuman_prob_dist_precise <- DAISIE:::DAISIE_margprobdist2(
  tvec = time_points,
  pars = pars_precise,
  M = mainland_pool_size,
  initEI = current_diversity_by_lineage,
  res = max_species,
  ddmodel = 11,
  methode = "ode45",
  reltolint = 1e-16,
  abstolint = 1e-16
)

# Calculate quantiles for total species
total_species_quantiles_c2p_precise <- calculate_species_quantiles(
  current_to_prehuman_prob_dist_precise,
  species_index = 3,  # Total species
  max_species = max_species,
  probability_levels = probability_levels
)

# Calculate recovery times
recovery_time_extinct_precise <- calculate_recovery_time(
  prehuman_diversity,
  current_to_prehuman_exp_precise[[3]],
  time_step
)

recovery_time_extinct_95_precise <- calculate_recovery_time(
  prehuman_diversity * 0.95,
  current_to_prehuman_exp_precise[[3]],
  time_step
)

#---
# 5.2 Future IUCN to Current for Precise M1
#---

# Convert the lineage-specific data to total counts for DAISIE_margprobdist
future_to_current_start_precise <- c(sum(future_diversity_by_lineage[,1]), 
                                     sum(future_diversity_by_lineage[,2]))

# Calculate expected species numbers over time using ExpEIN2 for M1
future_to_current_exp_precise <- DAISIE_ExpEIN2(
  tvec = time_points,
  pars = pars_precise,
  M = mainland_pool_size,
  initEI = future_diversity_by_lineage
)

future_to_current_prob_dist_precise <- DAISIE:::DAISIE_margprobdist2(
  tvec = time_points,
  pars = pars_precise,
  M = mainland_pool_size,
  initEI = future_diversity_by_lineage,
  res = max_species,
  ddmodel = 11,
  methode = "ode45",
  reltolint = 1e-16,
  abstolint = 1e-16
)

# Calculate quantiles for total species
total_species_quantiles_f2c_precise <- calculate_species_quantiles(
  future_to_current_prob_dist_precise,
  species_index = 3,  # Total species
  max_species = max_species,
  probability_levels = probability_levels
)

# Calculate recovery times
recovery_time_future_precise <- calculate_recovery_time(
  current_diversity_total,
  future_to_current_exp_precise[[3]],
  time_step
)

recovery_time_future_95_precise <- calculate_recovery_time(
  current_diversity_total * 0.95,
  future_to_current_exp_precise[[3]],
  time_step
)

recovery_time_future_prehuman_precise <- calculate_recovery_time(
  prehuman_diversity,
  future_to_current_exp_precise[[3]],
  time_step
)

recovery_time_future_prehuman_95_precise <- calculate_recovery_time(
  prehuman_diversity * 0.95,
  future_to_current_exp_precise[[3]],
  time_step
)

#-----------------------------------
# 6. SAVE RESULTS
#-----------------------------------

# Create a data frame to store the ERT results for both models
ert_results <- data.frame(
  Scenario = rep(c(
    "Recovery from current to pre-human diversity",
    "Recovery from current to 95% of pre-human diversity",
    "Recovery from future IUCN to current diversity",
    "Recovery from future IUCN to 95% of current diversity",
    "Recovery from future IUCN to pre-human diversity",
    "Recovery from future IUCN to 95% of pre-human diversity"
  ), 2),
  Model = rep(c("Max Age M1 (DD)", "Precise M1 (DD)"), each = 6),
  Time_Myr = c(
    # Maxage results
    recovery_time_extinct_maxage,
    recovery_time_extinct_95_maxage,
    recovery_time_future_maxage,
    recovery_time_future_95_maxage,
    recovery_time_future_prehuman_maxage,
    recovery_time_future_prehuman_95_maxage,
    # Precise results
    recovery_time_extinct_precise,
    recovery_time_extinct_95_precise,
    recovery_time_future_precise,
    recovery_time_future_95_precise,
    recovery_time_future_prehuman_precise,
    recovery_time_future_prehuman_95_precise
  )
)

print(ert_results)

write_xlsx(ert_results, "../outputs/DAISIE5_ERT_results.xlsx")

#-----------------------------------
# 7. PLOT ERTS
#-----------------------------------

#---
# 7.1 Prepare Plot Data for Max Age
#---

current_to_prehuman_plot_data_maxage <- data.frame(
  Time = time_points,
  Mean = current_to_prehuman_exp_maxage[[3]],
  Lower_95 = total_species_quantiles_c2p_maxage[, 2],
  Lower_90 = total_species_quantiles_c2p_maxage[, 3],
  Median = total_species_quantiles_c2p_maxage[, 4],
  Upper_90 = total_species_quantiles_c2p_maxage[, 5],
  Upper_95 = total_species_quantiles_c2p_maxage[, 6]
)

# Create plot data frame for scenario 2 (future to current)
future_to_current_plot_data_maxage <- data.frame(
  Time = time_points,
  Mean = future_to_current_exp_maxage[[3]],
  Lower_95 = total_species_quantiles_f2c_maxage[, 2],
  Lower_90 = total_species_quantiles_f2c_maxage[, 3],
  Median = total_species_quantiles_f2c_maxage[, 4],
  Upper_90 = total_species_quantiles_f2c_maxage[, 5],
  Upper_95 = total_species_quantiles_f2c_maxage[, 6]
)

#---
# 7.2 Prepare Plot Data for Precise
#---

# Create plot data frame for scenario 1 (current to pre-human)
current_to_prehuman_plot_data_precise <- data.frame(
  Time = time_points,
  Mean = current_to_prehuman_exp_precise[[3]],
  Lower_95 = total_species_quantiles_c2p_precise[, 2],
  Lower_90 = total_species_quantiles_c2p_precise[, 3],
  Median = total_species_quantiles_c2p_precise[, 4],
  Upper_90 = total_species_quantiles_c2p_precise[, 5],
  Upper_95 = total_species_quantiles_c2p_precise[, 6]
)

# Create plot data frame for scenario 2 (future to current)
future_to_current_plot_data_precise <- data.frame(
  Time = time_points,
  Mean = future_to_current_exp_precise[[3]],
  Lower_95 = total_species_quantiles_f2c_precise[, 2],
  Lower_90 = total_species_quantiles_f2c_precise[, 3],
  Median = total_species_quantiles_f2c_precise[, 4],
  Upper_90 = total_species_quantiles_f2c_precise[, 5],
  Upper_95 = total_species_quantiles_f2c_precise[, 6]
)

#---
# 7.3 Create Plots for Max Age M1 (Diversity-Dependent)
#---

# Create Plot 1: Current to pre-human for Max Age
current_to_prehuman_plot_maxage <- ggplot(current_to_prehuman_plot_data_maxage, aes(x = Time)) +
  geom_ribbon(aes(ymin = Lower_95, ymax = Upper_95), fill = "lightsteelblue1", alpha = 0.5) +
  geom_ribbon(aes(ymin = Lower_90, ymax = Upper_90), fill = "lightsteelblue2", alpha = 0.5) +
  geom_line(aes(y = Mean), color = "black", size = 0.6) +
  geom_hline(yintercept = prehuman_diversity, color = "darkblue", linetype = "dashed") +
  geom_hline(yintercept = current_diversity_total, color = "#C9202E", linetype = "dashed") +
  geom_vline(xintercept = recovery_time_extinct_maxage, color = "darkgrey", linetype = "dotted") +
  #annotate("text", x = recovery_time_extinct_maxage + 0.05, y = 20,
  #label = paste0("ERT = ", round(recovery_time_extinct_maxage, 2), " Myr"),
  #hjust = 0, vjust = 0, size = 3.5) +
  labs(
    title = "Max Age M1: Recovery from Current to Pre-European",
    x = "Time (Myr)",
    y = "Total Species"
  ) +
  theme_bw() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 200)) +
  scale_x_continuous(expand = c(0, 0)) +  
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  )
current_to_prehuman_plot_maxage

# Create Plot 2: Future IUCN to current for Max Age
future_to_current_plot_maxage <- ggplot(future_to_current_plot_data_maxage, aes(x = Time)) +
  geom_ribbon(aes(ymin = Lower_95, ymax = Upper_95), fill = "lightsteelblue1", alpha = 0.5) +
  geom_ribbon(aes(ymin = Lower_90, ymax = Upper_90), fill = "lightsteelblue2", alpha = 0.5) +
  geom_line(aes(y = Mean), color = "black", size = 0.6) +
  geom_hline(yintercept = prehuman_diversity, color = "darkblue", linetype = "dashed") +
  geom_hline(yintercept = current_diversity_total, color = "#C9202E", linetype = "dashed") +
  geom_hline(yintercept = future_iucn_total, color = "#F7941D", linetype = "dashed") +
  geom_vline(xintercept = recovery_time_future_maxage, color = "grey", linetype = "dotted") +
  geom_vline(xintercept = recovery_time_future_prehuman_maxage, color = "darkgrey", linetype = "dotted") +
  #annotate("text", x = recovery_time_future_maxage + 0.05, y = 20,
  #label = paste0("ERT (to current) = ", round(recovery_time_future_maxage, 2), " Myr"),
  # hjust = 0, vjust = 0, size = 3.5) +
  #annotate("text", x = recovery_time_future_prehuman_maxage + 0.05, y = 30,
  #label = paste0("ERT (to pre-human) = ", round(recovery_time_future_prehuman_maxage, 2), " Myr"),
  #hjust = 0, vjust = 0, size = 3.5) +
  labs(
    title = "Max Age M1: Recovery from Future IUCN to Current",
    x = "Time (Myr)",
    y = "Total Species") +
  theme_bw() +
  coord_cartesian(xlim = c(0, 3), ylim = c(0, 200)) +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  )
future_to_current_plot_maxage

#---
# 7.4 Create Plots for Precise M1 (Diversity-Dependent)
#---

# The confidence do not plot as expected due to DAISIE_margprobdist not being compatible with DD models

# Create Plot 1: Current to pre-human for Precise
current_to_prehuman_plot_precise <- ggplot(current_to_prehuman_plot_data_precise, aes(x = Time)) +
  geom_ribbon(aes(ymin = Lower_95, ymax = Upper_95), fill = "thistle1", alpha = 0.5) +
  geom_ribbon(aes(ymin = Lower_90, ymax = Upper_90), fill = "thistle2", alpha = 0.5) +
  geom_line(aes(y = Mean), color = "black", size = 0.6) +
  geom_hline(yintercept = prehuman_diversity, color = "darkblue", linetype = "dashed") +
  geom_hline(yintercept = current_diversity_total, color = "#C9202E", linetype = "dashed") +
  geom_vline(xintercept = recovery_time_extinct_precise, color = "darkgrey", linetype = "dotted") +
  #annotate("text", x = recovery_time_extinct_precise + 0.05, y = 20,
  #label = paste0("ERT = ", round(recovery_time_extinct_precise, 2), " Myr"),
  #hjust = 0, vjust = 0, size = 3.5) +
  labs(
    title = "Precise M1: Recovery from Current to Pre-European",
    x = "Time (Myr)",
    y = "Total Species") +
  theme_bw() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 200)) +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  )
current_to_prehuman_plot_precise

# Create Plot 2: Future IUCN to current for Precise
future_to_current_plot_precise <- ggplot(future_to_current_plot_data_precise, aes(x = Time)) +
  geom_ribbon(aes(ymin = Lower_95, ymax = Upper_95), fill = "thistle1", alpha = 0.5) +
  geom_ribbon(aes(ymin = Lower_90, ymax = Upper_90), fill = "thistle2", alpha = 0.5) +
  geom_line(aes(y = Mean), color = "black", size = 0.6) +
  geom_hline(yintercept = prehuman_diversity, color = "darkblue", linetype = "dashed") +
  geom_hline(yintercept = current_diversity_total, color = "#C9202E", linetype = "dashed") +
  geom_hline(yintercept = future_iucn_total, color = "#F7941D", linetype = "dashed") +
  geom_vline(xintercept = recovery_time_future_precise, color = "grey", linetype = "dotted") +
  geom_vline(xintercept = recovery_time_future_prehuman_precise, color = "darkgrey", linetype = "dotted") +
  #annotate("text", x = recovery_time_future_precise + 0.05, y = 20,
  #label = paste0("ERT (to current) = ", round(recovery_time_future_precise, 2), " Myr"),
  #hjust = 0, vjust = 0, size = 3.5) +
  #annotate("text", x = recovery_time_future_prehuman_precise + 0.05, y = 30,
  #label = paste0("ERT (to pre-human) = ", round(recovery_time_future_prehuman_precise, 2), " Myr"),
  #hjust = 0, vjust = 0, size = 3.5) +
  labs(
    title = "Precise M1: Recovery from Future IUCN to Current",
    x = "Time (Myr)",
    y = "Total Species") +
  theme_bw() +
  coord_cartesian(xlim = c(0, 3), ylim = c(0, 200)) +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  )
future_to_current_plot_precise

#---
# 7.5 Compare Mean return times between Maxage and Precise
#---

# Combined plot for Current to Pre-human recovery (Max Age and Precise)
current_to_prehuman_combined <- ggplot() +
  # Add mean lines only
  geom_line(data = current_to_prehuman_plot_data_maxage, aes(x = Time, y = Mean, color = "MaxAge"), size = 0.8) +
  geom_line(data = current_to_prehuman_plot_data_precise, aes(x = Time, y = Mean, color = "Precise"), size = 0.8) +
  
  # Reference lines
  geom_hline(yintercept = prehuman_diversity, color = "darkblue", linetype = "dashed") +
  geom_hline(yintercept = current_diversity_total, color = "#C9202E", linetype = "dashed") +
  
  # Recovery time markers
  geom_vline(xintercept = recovery_time_extinct_maxage, color = "#385582", linetype = "dotted") +
  geom_vline(xintercept = recovery_time_extinct_precise, color = "#754C29", linetype = "dotted") +
  
  # Plot styling
  labs(
    title = "Comparison of Models Max Age and Precise",
    subtitle = "Recovery from Current to Pre-human Diversity",
    x = "Time (Myr)",
    y = "Total Species",
    color = "Model"
  ) +
  scale_color_manual(values = c("MaxAge" = "#385582", "Precise" = "#754C29")) +
  theme_bw() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 200)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "bottom"
  )
current_to_prehuman_combined

# Combined plot for Threatened to Current & Pre-human recovery (Max Age and Precise)
future_to_current_combined <- ggplot() +
  # Add mean lines only
  geom_line(data = future_to_current_plot_data_maxage, aes(x = Time, y = Mean, color = "MaxAge"), size = 0.8) +
  geom_line(data = future_to_current_plot_data_precise, aes(x = Time, y = Mean, color = "Precise"), size = 0.8) +
  
  # Reference lines
  geom_hline(yintercept = prehuman_diversity, color = "darkblue", linetype = "dashed") +
  geom_hline(yintercept = current_diversity_total, color = "#C9202E", linetype = "dashed") +
  geom_hline(yintercept = future_iucn_total, color = "#F7941D", linetype = "dashed") +
  
  # Recovery time markers
  geom_vline(xintercept = recovery_time_future_maxage, color = "#385582", linetype = "dotted") +
  geom_vline(xintercept = recovery_time_future_prehuman_maxage, color = "#385582", linetype = "dotted") +
  geom_vline(xintercept = recovery_time_future_precise, color = "#754C29", linetype = "dotted") +
  geom_vline(xintercept = recovery_time_future_prehuman_precise, color = "#754C29", linetype = "dotted") +
  
  # Plot styling
  labs(
    subtitle = "Recovery from Threatened to Current & Pre-human Diversity",
    x = "Time (Myr)",
    y = "Total Species",
    color = "Model"
  ) +
  scale_color_manual(values = c("MaxAge" = "#385582", "Precise" = "#754C29")) +
  theme_bw() +
  coord_cartesian(xlim = c(0, 3), ylim = c(0, 200)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "bottom"
  )
future_to_current_combined

#---
# 7.6 Combine and Save Plots
#---

ERT_plots_maxage <- current_to_prehuman_plot_maxage + future_to_current_plot_maxage
ERT_plots_maxage

ERT_plots_precise <- current_to_prehuman_plot_precise + future_to_current_plot_precise
ERT_plots_precise

ERT_plots_maxage_precise <- current_to_prehuman_combined + future_to_current_combined
ERT_plots_maxage_precise

ggsave("../outputs/DAISIE5_ERT_plots_maxage_6.15myr_raw.pdf", ERT_plots_maxage, width = 8, height = 4)
ggsave("../outputs/DAISIE5_ERT_plots_precise_6.15myr_raw.pdf", ERT_plots_precise, width = 8, height = 4)
ggsave("../outputs/DAISIE5_ERT_plots_maxage_precise_6.15myr_raw.pdf", ERT_plots_maxage_precise , width = 6.5, height = 4)
