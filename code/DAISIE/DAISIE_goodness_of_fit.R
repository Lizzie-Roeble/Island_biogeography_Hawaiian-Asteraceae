################################################
# DAISIE: Goodness-of-fit     
###############################################

# Load required libraries
library(DAISIE)
library(ggplot2)
library(patchwork)

# Load simulations: 1000 simulations filtered on condition of >= 1 nonendemic
load("../outputs/DAISIE_sims_precise_m1_filt_1000.Rdata")
load("../outputs/DAISIE_sims_maxage_m1_filt_1000.Rdata")
island_replicates_maxage_m1 <- DAISIE_maxage_m1_sims_filt_1000
island_replicates_precise_m1 <- DAISIE_precise_m1_sims_filt_1000

# Observed diversity values
number_species_data <- 101
colonisations_data <- 12
endemic_species_data <- 100
non_endemic_species_data <- 1
endemic_clades_data <- 8 # cladogenetic clades / radiations
endemic_singletons_data <- 3 # also known as anagenetic spp
size_largest_clade_data <- 34

# Function to process simulation results
process_simulation <- function(island_replicates) {
  replicates <- length(island_replicates)
  
  # Initialize vectors with proper length for better performance
  number_colonists <- numeric(replicates)
  number_spec <- numeric(replicates)
  number_stac2 <- numeric(replicates)
  number_stac3 <- numeric(replicates)
  number_stac4 <- numeric(replicates)
  number_colonisations_including_stac3s <- numeric(replicates)
  number_stac2_singletons <- numeric(replicates)
  size_largest_clade <- numeric(replicates)
  size_smallest_clade <- numeric(replicates)
  
  # Fill with NA where appropriate
  size_largest_clade[] <- NA
  size_smallest_clade[] <- NA
  
  for(i in 1:replicates) {
    the_island <- island_replicates[[i]]
    number_colonists[i] <- the_island[[1]]$stt_all[nrow(the_island[[1]]$stt_all), "present"]
    
    if(number_colonists[i] == 0) {
      # Keep defaults (0 for numeric, NA for sizes)
      next
    }
    
    # Unlist for faster access to elements
    u_island <- unlist(the_island)
    
    # Calculate stats for inhabited islands
    number_spec[i] <- sum(the_island[[1]]$stt_all[nrow(the_island[[1]]$stt_all), c(2, 3, 4)])
    number_stac2[i] <- sum(u_island[names(u_island) == "stac"] == 2)
    number_stac3[i] <- sum(u_island[names(u_island) == "stac"] == 3)
    number_stac4[i] <- the_island[[1]]$stt_all[nrow(the_island[[1]]$stt_all), 2]
    
    # Process clade information if stac2 exists
    if(number_stac2[i] > 0) {
      # Preallocate lists with estimated size
      store_stac2s_brts <- vector("list", number_stac2[i])
      store_all_brts <- vector("list", length(the_island) - 1)
      
      stac2_counter <- 0
      all_counter <- 0
      
      for (isl in 1:(length(the_island) - 1)) {
        isl_data <- the_island[[isl + 1]]
        
        # Store branching times
        all_counter <- all_counter + 1
        store_all_brts[[all_counter]] <- isl_data$branching_times
        
        # Store stac2 branching times
        if(isl_data$stac == 2) {
          stac2_counter <- stac2_counter + 1
          store_stac2s_brts[[stac2_counter]] <- isl_data$branching_times
        }
      }
      
      # Trim lists to actual size
      if(stac2_counter < length(store_stac2s_brts)) {
        store_stac2s_brts <- store_stac2s_brts[1:stac2_counter]
      }
      if(all_counter < length(store_all_brts)) {
        store_all_brts <- store_all_brts[1:all_counter]
      }
      
      # Calculate metrics
      stac2_lengths <- lengths(store_stac2s_brts)
      number_stac2_singletons[i] <- sum(stac2_lengths == 2)
      size_largest_clade[i] <- max(stac2_lengths) - 1
      size_smallest_clade[i] <- min(lengths(store_all_brts)) - 1
    }
    
    # Handle colonisations including stac3s
    number_colonisations_including_stac3s[i] <- number_colonists[i]
    
    if(number_stac3[i] > 0) {
      count_stac3_cols <- 0
      for (isl in 1:(length(the_island) - 1)) {
        if(the_island[[isl + 1]]$stac == 3) {
          count_stac3_cols <- count_stac3_cols + length(the_island[[isl + 1]]$all_colonisations) - 1
        }
      }
      number_colonisations_including_stac3s[i] <- number_colonisations_including_stac3s[i] + count_stac3_cols
    }
  }
  
  # Return results without type1/type2 distinction
  return(list(
    number_spec = number_spec,
    number_colonists = number_colonists,
    number_stac2 = number_stac2,
    number_stac3 = number_stac3,
    number_stac4 = number_stac4,
    number_stac2_singletons = number_stac2_singletons,
    number_colonisations_including_stac3s = number_colonisations_including_stac3s,
    number_endemic_species = number_spec - number_stac4,
    size_largest_clade = size_largest_clade,
    size_smallest_clade = size_smallest_clade
  ))
}

# Process both parameter sets
results_maxage_m1 <- process_simulation(island_replicates_maxage_m1)
results_precise_m1 <- process_simulation(island_replicates_precise_m1)

# Histogram of simulated vs observed

# Create individual plots for maxage_m1
# Number of species (maxage_m1)
maxage_m1_species <- ggplot(data.frame(value = results_maxage_m1$number_spec), aes(x = value)) +
  geom_histogram(fill = "grey", color = NA, alpha = 0.7, bins = 20) +
  geom_vline(xintercept = median(results_maxage_m1$number_spec, na.rm = TRUE), color = "black", linewidth = 1) +
  geom_vline(xintercept = number_species_data, 
             color = "red", linewidth = 1) +
  labs(x = "Number of species", y = "Frequency") +
  theme_classic() +
  theme(plot.margin = margin(5, 5, 5, 5))
maxage_m1_species

# Colonization events (maxage_m1)
maxage_m1_colonists <- ggplot(data.frame(value = results_maxage_m1$number_colonists), aes(x = value)) +
  geom_histogram(fill = "grey", color = NA, alpha = 0.7, bins = 10) +
  geom_vline(xintercept = median(results_maxage_m1$number_colonists, na.rm = TRUE), color = "black", linewidth = 1) +
  geom_vline(xintercept = colonisations_data, 
             color = "red", linewidth = 1) +
  labs(x = "Number of colonization events", y = "Frequency") +
  theme_classic() +
  theme(plot.margin = margin(5, 5, 5, 5))
maxage_m1_colonists 

# Size of largest clade (maxage_m1)
maxage_m1_largest_clade <- ggplot(data.frame(value = results_maxage_m1$size_largest_clade), aes(x = value)) +
  geom_histogram(fill = "grey", color = NA, alpha = 0.7, bins = 15) +
  geom_vline(xintercept = median(results_maxage_m1$size_largest_clade, na.rm = TRUE), 
             color = "black", linewidth = 1) +
  geom_vline(xintercept = size_largest_clade_data, 
             color = "red", linewidth = 1) +
  labs(x = "Size of largest clade", y = "Frequency") +
  theme_classic() +
  theme(plot.margin = margin(5, 5, 5, 5))
maxage_m1_largest_clade

# Endemic clades (maxage_m1)
maxage_m1_endemic_clades <- ggplot(data.frame(value = results_maxage_m1$number_stac2), aes(x = value)) +
  geom_histogram(fill = "grey", color = NA, alpha = 0.7, bins = 12) +
  geom_vline(xintercept = median(results_maxage_m1$number_stac2, na.rm = TRUE), color = "black", linewidth = 1) +
  geom_vline(xintercept = endemic_clades_data, 
             color = "red", linewidth = 1) +
  labs(x = "Number of endemic clades", y = "Frequency") +
  theme_classic() +
  theme(plot.margin = margin(5, 5, 5, 5))
maxage_m1_endemic_clades

# Endemic singletons (maxage_m1)
maxage_m1_endemic_singletons <- ggplot(data.frame(value = results_maxage_m1$number_stac2_singletons), aes(x = value)) +
  geom_histogram(fill = "grey", color = NA, alpha = 0.7, bins = 10) +
  geom_vline(xintercept = median(results_maxage_m1$number_stac2_singletons, na.rm = TRUE), color = "black", linewidth = 1) +
  geom_vline(xintercept = endemic_singletons_data, 
             color = "red", linewidth = 1) +
  labs(x = "Endemic singletons", y = "Frequency") +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  theme_classic() +
  theme(plot.margin = margin(5, 5, 5, 5))
maxage_m1_endemic_singletons

# Non-endemic species (maxage_m1)
maxage_m1_non_endemic <- ggplot(data.frame(value = results_maxage_m1$number_stac4), aes(x = value)) +
  geom_histogram(fill = "grey", color = NA, alpha = 0.7, bins = 3) +
  geom_vline(xintercept = median(results_maxage_m1$number_stac4, na.rm = TRUE), color = "black", linewidth = 1) +
  geom_vline(xintercept = non_endemic_species_data, 
             color = "red", linewidth = 1) +
  labs(x = "Number of non-endemic species", y = "Frequency") +
  theme_classic() +
  theme(plot.margin = margin(5, 5, 5, 5))
maxage_m1_non_endemic

# Create individual plots for precise_m1
# Number of species (precise_m1)
precise_m1_species <- ggplot(data.frame(value = results_precise_m1$number_spec), aes(x = value)) +
  geom_histogram(fill = "grey", color = NA, alpha = 0.7, bins = 20) +
  geom_vline(xintercept = median(results_precise_m1$number_spec, na.rm = TRUE), color = "black", linewidth = 1) +
  geom_vline(xintercept = number_species_data, 
             color = "red", linewidth = 1) +
  labs(x = "Number of species", y = "Frequency") +
  theme_classic() +
  theme(plot.margin = margin(5, 5, 5, 5))
precise_m1_species

# Colonization events (precise_m1)
precise_m1_colonists <- ggplot(data.frame(value = results_precise_m1$number_colonists), aes(x = value)) +
  geom_histogram(fill = "grey", color = NA, alpha = 0.7, bins = 10) +
  geom_vline(xintercept = median(results_precise_m1$number_colonists, na.rm = TRUE), color = "black", linewidth = 1) +
  geom_vline(xintercept = colonisations_data, 
             color = "red", linewidth = 1) +
  labs(x = "Number of colonization events", y = "Frequency") +
  theme_classic() +
  theme(plot.margin = margin(5, 5, 5, 5))
precise_m1_colonists 

# Size of largest clade (precise_m1)
precise_m1_largest_clade <- ggplot(data.frame(value = results_precise_m1$size_largest_clade), aes(x = value)) +
  geom_histogram(fill = "grey", color = NA, alpha = 0.7, bins = 15) +
  geom_vline(xintercept = median(results_precise_m1$size_largest_clade, na.rm = TRUE), 
             color = "black", linewidth = 1) +
  geom_vline(xintercept = size_largest_clade_data, 
             color = "red", linewidth = 1) +
  labs(x = "Size of largest clade", y = "Frequency") +
  theme_classic() +
  theme(plot.margin = margin(5, 5, 5, 5))
precise_m1_largest_clade

# Endemic clades (precise_m1)
precise_m1_endemic_clades <- ggplot(data.frame(value = results_precise_m1$number_stac2), aes(x = value)) +
  geom_histogram(fill = "grey", color = NA, alpha = 0.7, bins = 12) +
  geom_vline(xintercept = median(results_precise_m1$number_stac2, na.rm = TRUE), color = "black", linewidth = 1) +
  geom_vline(xintercept = endemic_clades_data, 
             color = "red", linewidth = 1) +
  labs(x = "Number of endemic clades", y = "Frequency") +
  theme_classic() +
  theme(plot.margin = margin(5, 5, 5, 5))
precise_m1_endemic_clades

# Endemic singletons (precise_m1)
precise_m1_endemic_singletons <- ggplot(data.frame(value = results_precise_m1$number_stac2_singletons), aes(x = value)) +
  geom_histogram(fill = "grey", color = NA, alpha = 0.7, bins = 9) +
  geom_vline(xintercept = median(results_precise_m1$number_stac2_singletons, na.rm = TRUE), color = "black", linewidth = 1) +
  geom_vline(xintercept = endemic_singletons_data, 
             color = "red", linewidth = 1) +
  labs(x = "Endemic singletons", y = "Frequency") +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) + 
  theme_classic() +
  theme(plot.margin = margin(5, 5, 5, 5))
precise_m1_endemic_singletons

# Non-endemic species (precise_m1)
precise_m1_non_endemic <- ggplot(data.frame(value = results_precise_m1$number_stac4), aes(x = value)) +
  geom_histogram(fill = "grey", color = NA, alpha = 0.7, bins = 5) +
  geom_vline(xintercept = median(results_precise_m1$number_stac4, na.rm = TRUE), color = "black", linewidth = 1) +
  geom_vline(xintercept = non_endemic_species_data, 
             color = "red", linewidth = 1) +
  labs(x = "Number of non-endemic species", y = "Frequency") +
  theme_classic() +
  theme(plot.margin = margin(5, 5, 5, 5))
precise_m1_non_endemic

# Combine
goodness_fit_maxage_m1 <- maxage_m1_species + maxage_m1_colonists + maxage_m1_largest_clade + maxage_m1_endemic_clades + maxage_m1_endemic_singletons + maxage_m1_non_endemic
goodness_fit_maxage_m1

goodness_fit_precise_m1 <- precise_m1_species + precise_m1_colonists + precise_m1_largest_clade + precise_m1_endemic_clades + precise_m1_endemic_singletons + precise_m1_non_endemic
goodness_fit_precise_m1

# Save the plots 
ggsave("../outputs/DAISIE_goodness_of_fit_maxage_m1_raw.pdf", goodness_fit_maxage_m1, 
       width = 11, height = 7, units = "in")
ggsave("../outputs/DAISIE_goodness_of_fit_precise_m1_raw.pdf", goodness_fit_precise_m1, 
       width = 11, height = 7, units = "in")

########### not updated below ############
# Additional exploration

# Create a summary dataframe for quick comparison

summary_stats <- data.frame(
  Parameter = c("Number of Species", "Colonizations", "Endemic Species", "Endemic Clades", 
                "Endemic Singletons", "Non-endemic Species", "Size Largest Clade"),
  Observed = c(number_species_data, colonisations_data, endemic_species_data, 
               endemic_clades_data, endemic_singletons_data, non_endemic_species_data, 
               size_largest_clade_data),
  Median_m1 = c(median(results_m1$number_spec), median(results_m1$number_colonists), 
                median(results_m1$number_endemic_species), median(results_m1$number_stac2), 
                median(results_m1$number_stac2_singletons), median(results_m1$number_stac4),
                median(results_m1$size_largest_clade, na.rm = TRUE)),
  Median_m4 = c(median(results_m4$number_spec), median(results_m4$number_colonists), 
                median(results_m4$number_endemic_species), median(results_m4$number_stac2), 
                median(results_m4$number_stac2_singletons), median(results_m4$number_stac4),
                median(results_m4$size_largest_clade, na.rm = TRUE))
)
# Add absolute differences from observed
summary_stats$Diff_m1 <- abs(summary_stats$Median_m1 - summary_stats$Observed)
summary_stats$Diff_m4 <- abs(summary_stats$Median_m4 - summary_stats$Observed)
# Better model (lower difference is better)
summary_stats$Better_Model <- ifelse(summary_stats$Diff_m1 < summary_stats$Diff_m4, "m1", 
                                     ifelse(summary_stats$Diff_m1 > summary_stats$Diff_m4, "m4", "Tie"))
# Print summary statistics
print(summary_stats)
# Write to csv
write.csv(summary_stats, "Goodness-of-fit_model_comparison.csv", row.names = FALSE)

# Create comparison bar plot

summary_long <- melt(summary_stats[, c("Parameter", "Observed", "Median_m1", "Median_m4")], 
                     id.vars = "Parameter")
# Create bar plot
comparison_plot <- ggplot(summary_long, aes(x = Parameter, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.8) +
  scale_fill_manual(values = c("Observed" = "darkgrey", "Median_m1" = "#385582", "Median_m4" = "#754C29"),
                    labels = c("Observed" = "Observed", "Median_m1" = "Model m1", "Median_m4" = "Model m4")) +
  labs(title = "Comparison of Models m1 and m4 with Observed Data",
       x = "", y = "Value", fill = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "top"
  )
comparison_plot
# Save the comparison plot
#ggsave("Hawaii_model_comparison.png", comparison_plot, width = 10, height = 6, dpi = 300)

# Create a heatmap of differences

diff_matrix <- as.matrix(summary_stats[, c("Diff_m1", "Diff_m4")])
rownames(diff_matrix) <- summary_stats$Parameter
# Melt the difference matrix for ggplot
diff_data <- melt(diff_matrix)
colnames(diff_data) <- c("Parameter", "Model", "Difference")
diff_data$Model <- gsub("Diff_", "", diff_data$Model)
# Create heatmap
heatmap_plot <- ggplot(diff_data, aes(x = Model, y = Parameter, fill = Difference)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(diff_data$Difference)) +
  geom_text(aes(label = round(Difference, 1)), color = "black", size = 3) +
  labs(title = "Absolute Differences from Observed Values",
       fill = "Difference") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_text(size = 14, face = "bold")
  )
# Save the heatmap
#ggsave("Hawaii_difference_heatmap.png", heatmap_plot, width = 8, height = 6, dpi = 300)

