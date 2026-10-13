###############################################
# DAISIE: Goodness-of-fit 
# DAISIE5 
###############################################

# Goal: Compare simulated diversity (1000 filtered sims per model) to observed Hawaii data
# Datasets: maxage and precise, 5.14myr and 6.15myr

# Load required libraries
library(DAISIE)
library(ggplot2)
library(patchwork)
library(reshape2)

# Plot font (Gill Sans; on Windows use "Gill Sans MT")
plot_font <- "Gill Sans"

# PDF device using macOS Quartz (supports system fonts; no XQuartz/cairo needed)
pdf_device <- function(filename, width, height, ...) {
  grDevices::quartz(file = filename, type = "pdf", width = width, height = height)
}

# Load simulations: 1000 simulations filtered on condition of >= 1 nonendemic
load("../data/DAISIE5_sims/DAISIE_maxage_5.14myr_sims_filt_1000.Rdata")
load("../data/DAISIE5_sims/DAISIE_precise_5.14myr_sims_filt_1000.Rdata")
load("../data/DAISIE5_sims/DAISIE_maxage_6.15myr_sims_filt_1000.Rdata")
load("../data/DAISIE5_sims/DAISIE_precise_6.15myr_sims_filt_1000.Rdata")

island_replicates_maxage_5.14myr <- DAISIE_maxage_5.14myr_m1_sims_filt_1000
island_replicates_precise_5.14myr <- DAISIE_precise_5.14myr_m1_sims_filt_1000
island_replicates_maxage_6.15myr <- DAISIE_maxage_6.15myr_m1_sims_filt_1000
island_replicates_precise_6.15myr <- DAISIE_precise_6.15myr_m1_sims_filt_1000

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

# Process all four datasets
results_maxage_5.14myr <- process_simulation(island_replicates_maxage_5.14myr)
results_precise_5.14myr <- process_simulation(island_replicates_precise_5.14myr)
results_maxage_6.15myr <- process_simulation(island_replicates_maxage_6.15myr)
results_precise_6.15myr <- process_simulation(island_replicates_precise_6.15myr)

# List of results, names used for labelling below (keeps datasets from getting mixed up)
results_list <- list(
  maxage_5.14myr = results_maxage_5.14myr,
  precise_5.14myr = results_precise_5.14myr,
  maxage_6.15myr = results_maxage_6.15myr,
  precise_6.15myr = results_precise_6.15myr
)

# Sanity checks

## one row per dataset: 1000 sims, every sim has >= 1 non-endemic, and medians differ between datasets
data.frame(
  Dataset = names(results_list),
  N_sims = sapply(results_list, function(x) length(x$number_spec)), # should be 1000
  Min_nonendemic = sapply(results_list, function(x) min(x$number_stac4)), # should be >= 1
  Median_species = sapply(results_list, function(x) median(x$number_spec)),
  Median_colonists = sapply(results_list, function(x) median(x$number_colonists))
)

## maxage and precise objects should not be identical (should be FALSE)
identical(island_replicates_maxage_5.14myr, island_replicates_precise_5.14myr)
identical(island_replicates_maxage_6.15myr, island_replicates_precise_6.15myr)
identical(island_replicates_maxage_5.14myr, island_replicates_maxage_6.15myr)
identical(island_replicates_precise_5.14myr, island_replicates_precise_6.15myr)

# Histogram of simulated vs observed

# Function to create the six goodness-of-fit plots for one dataset and combine them
# Small whole-number counts use bins of width 1 centred on integers
# binwidth_species and binwidth_clade default to a whole-number width giving ~20 and ~15 bars (set a number to override)
make_goodness_of_fit_plot <- function(results, binwidth_species = NULL, binwidth_clade = NULL) {
  
  if (is.null(binwidth_species)) {
    binwidth_species <- max(1, ceiling(diff(range(results$number_spec, na.rm = TRUE)) / 20))
  }
  if (is.null(binwidth_clade)) {
    binwidth_clade <- max(1, ceiling(diff(range(results$size_largest_clade, na.rm = TRUE)) / 15))
  }
  
  # Number of species
  species_plot <- ggplot(data.frame(value = results$number_spec), aes(x = value)) +
    geom_histogram(fill = "grey", color = NA, alpha = 0.7, binwidth = binwidth_species, boundary = 0) +
    geom_vline(xintercept = median(results$number_spec, na.rm = TRUE), color = "black", linewidth = 1) +
    geom_vline(xintercept = number_species_data, 
               color = "red", linewidth = 1) +
    labs(x = "Number of species", y = "Frequency") +
    theme_classic(base_family = plot_font) +
    theme(plot.margin = margin(5, 5, 5, 5))
  
  # Colonization events
  colonists_plot <- ggplot(data.frame(value = results$number_colonists), aes(x = value)) +
    geom_histogram(fill = "grey", color = NA, alpha = 0.7, binwidth = 1, center = 0) +
    geom_vline(xintercept = median(results$number_colonists, na.rm = TRUE), color = "black", linewidth = 1) +
    geom_vline(xintercept = colonisations_data, 
               color = "red", linewidth = 1) +
    labs(x = "Number of colonization events", y = "Frequency") +
    theme_classic(base_family = plot_font) +
    theme(plot.margin = margin(5, 5, 5, 5))
  
  # Size of largest clade
  largest_clade_plot <- ggplot(data.frame(value = results$size_largest_clade), aes(x = value)) +
    geom_histogram(fill = "grey", color = NA, alpha = 0.7, binwidth = binwidth_clade, boundary = 0) +
    geom_vline(xintercept = median(results$size_largest_clade, na.rm = TRUE), 
               color = "black", linewidth = 1) +
    geom_vline(xintercept = size_largest_clade_data, 
               color = "red", linewidth = 1) +
    labs(x = "Size of largest clade", y = "Frequency") +
    theme_classic(base_family = plot_font) +
    theme(plot.margin = margin(5, 5, 5, 5))
  
  # Endemic clades
  endemic_clades_plot <- ggplot(data.frame(value = results$number_stac2), aes(x = value)) +
    geom_histogram(fill = "grey", color = NA, alpha = 0.7, binwidth = 1, center = 0) +
    geom_vline(xintercept = median(results$number_stac2, na.rm = TRUE), color = "black", linewidth = 1) +
    geom_vline(xintercept = endemic_clades_data, 
               color = "red", linewidth = 1) +
    labs(x = "Number of endemic clades", y = "Frequency") +
    theme_classic(base_family = plot_font) +
    theme(plot.margin = margin(5, 5, 5, 5))
  
  # Endemic singletons
  endemic_singletons_plot <- ggplot(data.frame(value = results$number_stac2_singletons), aes(x = value)) +
    geom_histogram(fill = "grey", color = NA, alpha = 0.7, binwidth = 1, center = 0) +
    geom_vline(xintercept = median(results$number_stac2_singletons, na.rm = TRUE), color = "black", linewidth = 1) +
    geom_vline(xintercept = endemic_singletons_data, 
               color = "red", linewidth = 1) +
    labs(x = "Endemic singletons", y = "Frequency") +
    scale_x_continuous(breaks = seq(0, 10, by = 2)) +
    theme_classic(base_family = plot_font) +
    theme(plot.margin = margin(5, 5, 5, 5))
  
  # Non-endemic species
  non_endemic_plot <- ggplot(data.frame(value = results$number_stac4), aes(x = value)) +
    geom_histogram(fill = "grey", color = NA, alpha = 0.7, binwidth = 1, center = 0) +
    geom_vline(xintercept = median(results$number_stac4, na.rm = TRUE), color = "black", linewidth = 1) +
    geom_vline(xintercept = non_endemic_species_data, 
               color = "red", linewidth = 1) +
    labs(x = "Number of non-endemic species", y = "Frequency") +
    theme_classic(base_family = plot_font) +
    theme(plot.margin = margin(5, 5, 5, 5))
  
  # Combine
  species_plot + colonists_plot + largest_clade_plot + endemic_clades_plot + endemic_singletons_plot + non_endemic_plot
}

# Create combined plots for each dataset
goodness_fit_maxage_5.14myr <- make_goodness_of_fit_plot(results_maxage_5.14myr)
goodness_fit_maxage_5.14myr

goodness_fit_precise_5.14myr <- make_goodness_of_fit_plot(results_precise_5.14myr)
goodness_fit_precise_5.14myr

goodness_fit_maxage_6.15myr <- make_goodness_of_fit_plot(results_maxage_6.15myr)
goodness_fit_maxage_6.15myr

goodness_fit_precise_6.15myr <- make_goodness_of_fit_plot(results_precise_6.15myr)
goodness_fit_precise_6.15myr

# Save the plots 
ggsave("../outputs/DAISIE5_goodness_of_fit_maxage_5.14myr_raw.pdf", goodness_fit_maxage_5.14myr, 
       width = 11, height = 7, units = "in", device = pdf_device)
ggsave("../outputs/DAISIE5_goodness_of_fit_precise_5.14myr_raw.pdf", goodness_fit_precise_5.14myr, 
       width = 11, height = 7, units = "in", device = pdf_device)
ggsave("../outputs/DAISIE5_goodness_of_fit_maxage_6.15myr_raw.pdf", goodness_fit_maxage_6.15myr, 
       width = 11, height = 7, units = "in", device = pdf_device)
ggsave("../outputs/DAISIE5_goodness_of_fit_precise_6.15myr_raw.pdf", goodness_fit_precise_6.15myr, 
       width = 11, height = 7, units = "in", device = pdf_device)

########### Additional exploration ###########

# Create a summary dataframe for quick comparison
model_names <- names(results_list)

summary_stats <- data.frame(
  Parameter = c("Number of Species", "Colonizations", "Endemic Species", "Endemic Clades", 
                "Endemic Singletons", "Non-endemic Species", "Size Largest Clade"),
  Observed = c(number_species_data, colonisations_data, endemic_species_data, 
               endemic_clades_data, endemic_singletons_data, non_endemic_species_data,
               size_largest_clade_data)
)

# Add median of simulations for each dataset
for (model in model_names) {
  results <- results_list[[model]]
  summary_stats[[paste0("Median_", model)]] <- c(
    median(results$number_spec), median(results$number_colonists), 
    median(results$number_endemic_species), median(results$number_stac2), 
    median(results$number_stac2_singletons), median(results$number_stac4),
    median(results$size_largest_clade, na.rm = TRUE)
  )
}

# Add absolute differences from observed
for (model in model_names) {
  summary_stats[[paste0("Diff_", model)]] <- abs(summary_stats[[paste0("Median_", model)]] - summary_stats$Observed)
}

# Closest model per parameter (lower difference is better)
diff_cols <- paste0("Diff_", model_names)
summary_stats$Closest_Model <- apply(summary_stats[, diff_cols], 1, function(x) {
  best <- which(x == min(x))
  if (length(best) > 1) "Tie" else sub("Diff_", "", names(x)[best])
})

# Print summary statistics
print(summary_stats)

# Write to csv
write.csv(summary_stats, "../outputs/DAISIE5_goodness_of_fit_model_comparison.csv", row.names = FALSE)

# Create comparison bar plot
summary_long <- melt(summary_stats[, c("Parameter", "Observed", paste0("Median_", model_names))], 
                     id.vars = "Parameter")

# Create bar plot
comparison_plot <- ggplot(summary_long, aes(x = Parameter, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.8) +
  scale_fill_manual(values = c("Observed" = "darkgrey", 
                               "Median_maxage_5.14myr" = "#385582", 
                               "Median_maxage_6.15myr" = "#8FA9D1",
                               "Median_precise_5.14myr" = "#754C29",
                               "Median_precise_6.15myr" = "#C9A27E"),
                    labels = c("Observed" = "Observed", 
                               "Median_maxage_5.14myr" = "Maxage 5.14myr", 
                               "Median_maxage_6.15myr" = "Maxage 6.15myr",
                               "Median_precise_5.14myr" = "Precise 5.14myr",
                               "Median_precise_6.15myr" = "Precise 6.15myr")) +
  labs(title = "Comparison of DAISIE5 Models with Observed Data",
       x = "", y = "Value", fill = "") +
  theme_minimal(base_family = plot_font) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "top"
  )

comparison_plot

# Save the comparison plot
#ggsave("../outputs/DAISIE5_model_comparison.png", comparison_plot, width = 10, height = 6, dpi = 300)

# Create a heatmap of differences
diff_matrix <- as.matrix(summary_stats[, diff_cols])
rownames(diff_matrix) <- summary_stats$Parameter

# Melt the difference matrix for ggplot
diff_data <- melt(diff_matrix)
colnames(diff_data) <- c("Parameter", "Model", "Difference")
diff_data$Model <- gsub("Diff_", "", diff_data$Model)

# Create heatmap
heatmap_plot <- ggplot(diff_data, aes(x = Model, y = Parameter, fill = Difference)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(diff_data$Difference)) +
  geom_text(aes(label = round(Difference, 1)), color = "black", size = 3, family = plot_font) +
  labs(title = "Absolute Differences from Observed Values",
       fill = "Difference") +
  theme_minimal(base_family = plot_font) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14, face = "bold")
  )

heatmap_plot

# Save the heatmap
#ggsave("../outputs/DAISIE5_difference_heatmap.png", heatmap_plot, width = 8, height = 6, dpi = 300)
