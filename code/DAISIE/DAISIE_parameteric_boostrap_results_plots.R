###############################################
# DAISIE 
# Parametric bootstrap: results - plots and summary    
# Updated with DAISIE5
###############################################

# Goal 1: Plot the results of the parametric bootstrap analysis for all four datasets (maxage and precise, 5.14myr and 6.15myr)
# Goal 2: Summary table of ML parameter estimates, median value from simulations, and confidence intervals
# Note: plots from pre-filtered results (first 1000 replicates >= 1 nonendemic)

library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)
library(writexl)

# Plot font (Gill Sans; on Windows use "Gill Sans MT")
plot_font <- "Gill Sans"

# PDF device using macOS Quartz 
pdf_device <- function(filename, width, height, ...) {
  grDevices::quartz(file = filename, type = "pdf", width = width, height = height)
}

# Helper functions

# create pairwise scatter plot
create_scatter_plot <- function(data, x_var, y_var, x_true, y_true,
                                log_x = FALSE, log_y = FALSE) {
  
  # Lookup table for nice axis labels
  var_labels <- c(
    Cladogenesis = "Cladogenesis",
    Extinction = "Extinction",
    Carrying_Capacity = "Carrying Capacity",
    Colonization = "Colonization",
    Anagenesis = "Anagenesis"
  )
  
  # Basic plot with points and MLE point
  p <- ggplot(data, aes_string(x = x_var, y = y_var)) +
    geom_point(color = "lightblue", alpha = 0.5, size = 1.5) +
    geom_point(data = data.frame(x = x_true, y = y_true), aes(x = x, y = y),
               color = "red", size = 3) +
    labs(
      x = var_labels[[x_var]],
      y = var_labels[[y_var]]
    ) +
    theme_classic(base_family = plot_font)
  
  # Apply log scaling if requested
  if (log_x) {
    p <- p + scale_x_log10(labels = scales::label_number())
  }
  if (log_y) {
    p <- p + scale_y_log10(labels = scales::label_number())
  }
  
  # Improve Colonization axis formatting (avoid scientific notation, add accuracy)
  if (x_var == "Colonization") {
    p <- p + scale_x_continuous(labels = scales::label_number(accuracy = 0.0001))
  }
  if (y_var == "Colonization") {
    p <- p + scale_y_continuous(labels = scales::label_number(accuracy = 0.0001))
  }
  
  return(p)
}

# helper function to remove axis labels for plot matrix
remove_axis_labels <- function(plot, remove_x = TRUE, remove_y = TRUE) {
  plot +
    theme(
      axis.title.x = if (remove_x) element_blank() else element_text(),
      axis.title.y = if (remove_y) element_blank() else element_text()
    )
}

# Main bootstrap analysis function
run_daisie_bootstrap_analysis <- function(input_path, output_label) {
  
  # Load bootstrap results (already pre-filtered results)
  bootstrap_results <- read.table(input_path, header = FALSE)
  colnames(bootstrap_results) <- c(
    "Replicate", "Island_age", "Mainland_pool", "Type", "Res", "Model", 
    "Solver", "Cladogenesis", "Extinction", "Carrying_Capacity", 
    "Colonization", "Anagenesis", "Likelihood", "Number_of_parameters", 
    "Number_parameters_type2", "Cladogenesis_MLE", "Extinction_MLE", 
    "Carrying_Capacity_MLE", "Colonization_MLE", "Anagenesis_MLE"
  )
  
  # Report number of replicates used
  cat(output_label, "- number of replicates:", nrow(bootstrap_results), "\n")
  
  # Extract MLE values from first row
  original_mle <- c(
    Cladogenesis = bootstrap_results$Cladogenesis_MLE[1],
    Extinction = bootstrap_results$Extinction_MLE[1],
    Carrying_Capacity = bootstrap_results$Carrying_Capacity_MLE[1],
    Colonization = bootstrap_results$Colonization_MLE[1],
    Anagenesis = bootstrap_results$Anagenesis_MLE[1]
  )
  
  # Calculate medians
  parameter_stats <- bootstrap_results %>%
    summarize(
      Cladogenesis_median = median(Cladogenesis, na.rm = TRUE),
      Extinction_median = median(Extinction, na.rm = TRUE),
      K_median = median(Carrying_Capacity, na.rm = TRUE),
      Colonization_median = median(Colonization, na.rm = TRUE),
      Anagenesis_median = median(Anagenesis, na.rm = TRUE)
    )
  
  # Histograms
  clado_plot <- ggplot(data.frame(value = bootstrap_results$Cladogenesis), aes(x = value)) +
    geom_histogram(fill = "lightblue", color = NA, alpha = 0.7, bins = 20) +
    geom_vline(xintercept = parameter_stats$Cladogenesis_median, color = "black", linewidth = 1) +
    geom_vline(xintercept = original_mle["Cladogenesis"], color = "red", linewidth = 1) +
    labs(x = "Cladogenesis", y = "Cladogenesis") +
    theme_classic(base_family = plot_font) +
    theme(plot.margin = margin(5, 5, 5, 5))
  
  ext_plot <- ggplot(data.frame(value = bootstrap_results$Extinction), aes(x = value)) +
    geom_histogram(fill = "lightblue", color = NA, alpha = 0.7, bins = 20) +
    geom_vline(xintercept = parameter_stats$Extinction_median, color = "black", linewidth = 1) +
    geom_vline(xintercept = original_mle["Extinction"], color = "red", linewidth = 1) +
    labs(x = "Extinction", y = "Frequency") +
    theme_classic(base_family = plot_font) +
    theme(plot.margin = margin(5, 5, 5, 5))
  
  k_plot <- ggplot(data.frame(value = bootstrap_results$Carrying_Capacity), aes(x = value)) +
    geom_histogram(fill = "lightblue", color = NA, alpha = 0.7, bins = 25) +
    geom_vline(xintercept = parameter_stats$K_median, color = "black", linewidth = 1) +
    geom_vline(xintercept = original_mle["Carrying_Capacity"], color = "red", linewidth = 1) +
    scale_x_log10() +
    labs(x = "Carrying Capacity", y = "Frequency") +
    theme_classic(base_family = plot_font) +
    theme(plot.margin = margin(5, 5, 5, 5))
  
  col_plot <- ggplot(data.frame(value = bootstrap_results$Colonization), aes(x = value)) +
    geom_histogram(fill = "lightblue", color = NA, alpha = 0.7, bins = 20) +
    geom_vline(xintercept = parameter_stats$Colonization_median, color = "black", linewidth = 1) +
    geom_vline(xintercept = original_mle["Colonization"], color = "red", linewidth = 1) +
    labs(x = "Colonization", y = "Frequency") +
    theme_classic(base_family = plot_font) +
    theme(plot.margin = margin(5, 5, 5, 5))+
    scale_x_continuous(labels = scales::label_number(accuracy = 0.0001))
  
  ana_plot <- ggplot(data.frame(value = bootstrap_results$Anagenesis), aes(x = value)) +
    geom_histogram(fill = "lightblue", color = NA, alpha = 0.7, bins = 20) +
    geom_vline(xintercept = parameter_stats$Anagenesis_median, color = "black", linewidth = 1) +
    geom_vline(xintercept = original_mle["Anagenesis"], color = "red", linewidth = 1) +
    labs(x = "Anagenesis", y = "Frequency") +
    theme_classic(base_family = plot_font) +
    theme(plot.margin = margin(5, 5, 5, 5))
  
  hist_combined <- clado_plot + ext_plot + k_plot + col_plot + ana_plot
  
  ggsave(paste0("../outputs/DAISIE5_parametric_bootstrap_results_", output_label, "_histo_raw.pdf"), 
         hist_combined, width = 11, height = 7, units = "in", device = pdf_device)
  
  # Pairwise scatter matrix
  parameter_df <- bootstrap_results %>%
    select(Cladogenesis, Extinction, Carrying_Capacity, Colonization, Anagenesis) %>%
    filter(if_all(everything(), is.finite))
  
  clado_ext_plot <- create_scatter_plot(parameter_df, "Cladogenesis", "Extinction", original_mle["Cladogenesis"], original_mle["Extinction"])
  clado_k_plot <- create_scatter_plot(parameter_df, "Cladogenesis", "Carrying_Capacity", original_mle["Cladogenesis"], original_mle["Carrying_Capacity"], log_y = TRUE)
  clado_col_plot <- create_scatter_plot(parameter_df, "Cladogenesis", "Colonization", original_mle["Cladogenesis"], original_mle["Colonization"])
  clado_ana_plot <- create_scatter_plot(parameter_df, "Cladogenesis", "Anagenesis", original_mle["Cladogenesis"], original_mle["Anagenesis"])
  ext_k_plot <- create_scatter_plot(parameter_df, "Extinction", "Carrying_Capacity", original_mle["Extinction"], original_mle["Carrying_Capacity"], log_y = TRUE)
  ext_col_plot <- create_scatter_plot(parameter_df, "Extinction", "Colonization", original_mle["Extinction"], original_mle["Colonization"])
  ext_ana_plot <- create_scatter_plot(parameter_df, "Extinction", "Anagenesis", original_mle["Extinction"], original_mle["Anagenesis"])
  k_col_plot <- create_scatter_plot(parameter_df, "Carrying_Capacity", "Colonization", original_mle["Carrying_Capacity"], original_mle["Colonization"], log_x = TRUE)
  k_ana_plot <- create_scatter_plot(parameter_df, "Carrying_Capacity", "Anagenesis", original_mle["Carrying_Capacity"], original_mle["Anagenesis"], log_x = TRUE)
  col_ana_plot <- create_scatter_plot(parameter_df, "Colonization", "Anagenesis", original_mle["Colonization"], original_mle["Anagenesis"])
  
  # Apply label removal to all plots based on matrix position
  
  # Row 1
  clado_plot2 <- remove_axis_labels(clado_plot, remove_x = TRUE, remove_y = FALSE)
  
  # Row 2
  clado_ext_plot <- remove_axis_labels(clado_ext_plot, remove_x = TRUE, remove_y = FALSE)
  ext_plot2       <- remove_axis_labels(ext_plot, remove_x = TRUE, remove_y = TRUE)
  
  # Row 3
  clado_k_plot <- remove_axis_labels(clado_k_plot, remove_x = TRUE, remove_y = FALSE)
  ext_k_plot   <- remove_axis_labels(ext_k_plot, remove_x = TRUE, remove_y = TRUE)
  k_plot2       <- remove_axis_labels(k_plot, remove_x = TRUE, remove_y = TRUE)
  
  # Row 4
  clado_col_plot <- remove_axis_labels(clado_col_plot, remove_x = TRUE, remove_y = FALSE)
  ext_col_plot   <- remove_axis_labels(ext_col_plot, remove_x = TRUE, remove_y = TRUE)
  k_col_plot     <- remove_axis_labels(k_col_plot, remove_x = TRUE, remove_y = TRUE)
  col_plot2      <- remove_axis_labels(col_plot, remove_x = TRUE, remove_y = TRUE)
  
  # Row 5 (bottom row — keep x-axis)
  clado_ana_plot <- remove_axis_labels(clado_ana_plot, remove_x = FALSE, remove_y = FALSE)
  ext_ana_plot   <- remove_axis_labels(ext_ana_plot, remove_x = FALSE, remove_y = TRUE)
  k_ana_plot     <- remove_axis_labels(k_ana_plot, remove_x = FALSE, remove_y = TRUE)
  col_ana_plot   <- remove_axis_labels(col_ana_plot, remove_x = FALSE, remove_y = TRUE)
  ana_plot2       <- remove_axis_labels(ana_plot, remove_x = FALSE, remove_y = TRUE)
  
  blank_plot <- ggplot() + theme_void()
  
  parameter_matrix <- (
    clado_plot2 + blank_plot + blank_plot + blank_plot + blank_plot +
      clado_ext_plot + ext_plot2 + blank_plot + blank_plot + blank_plot +
      clado_k_plot + ext_k_plot + k_plot2 + blank_plot + blank_plot +
      clado_col_plot + ext_col_plot + k_col_plot + col_plot2 + blank_plot +
      clado_ana_plot + ext_ana_plot + k_ana_plot + col_ana_plot + ana_plot2
  ) + plot_layout(nrow = 5, ncol = 5)
  
  ggsave(paste0("../outputs/DAISIE5_parametric_bootstrap_results_", output_label, "_histo_scatter_raw.pdf"), 
         parameter_matrix, width = 11, height = 8, units = "in", device = pdf_device)
  
  # Confidence intervals
  ci_stats <- bootstrap_results %>%
    summarize(
      Cladogenesis_CI_lower = quantile(Cladogenesis, 0.025, na.rm = TRUE),
      Cladogenesis_CI_upper = quantile(Cladogenesis, 0.975, na.rm = TRUE),
      Extinction_CI_lower = quantile(Extinction, 0.025, na.rm = TRUE),
      Extinction_CI_upper = quantile(Extinction, 0.975, na.rm = TRUE),
      K_CI_lower = quantile(Carrying_Capacity, 0.025, na.rm = TRUE),
      K_CI_upper = quantile(Carrying_Capacity, 0.975, na.rm = TRUE),
      Colonization_CI_lower = quantile(Colonization, 0.025, na.rm = TRUE),
      Colonization_CI_upper = quantile(Colonization, 0.975, na.rm = TRUE),
      Anagenesis_CI_lower = quantile(Anagenesis, 0.025, na.rm = TRUE),
      Anagenesis_CI_upper = quantile(Anagenesis, 0.975, na.rm = TRUE)
    )
  
  summary_table <- data.frame(
    Parameter = c("Cladogenesis", "Extinction", "Carrying Capacity", "Colonization", "Anagenesis"),
    Original_MLE = original_mle,
    Bootstrap_Median = c(
      parameter_stats$Cladogenesis_median,
      parameter_stats$Extinction_median,
      parameter_stats$K_median,
      parameter_stats$Colonization_median,
      parameter_stats$Anagenesis_median
    ),
    CI_Lower = c(
      ci_stats$Cladogenesis_CI_lower,
      ci_stats$Extinction_CI_lower,
      ci_stats$K_CI_lower,
      ci_stats$Colonization_CI_lower,
      ci_stats$Anagenesis_CI_lower
    ),
    CI_Upper = c(
      ci_stats$Cladogenesis_CI_upper,
      ci_stats$Extinction_CI_upper,
      ci_stats$K_CI_upper,
      ci_stats$Colonization_CI_upper,
      ci_stats$Anagenesis_CI_upper
    )
  )
  
  write_xlsx(summary_table, paste0("../outputs/DAISIE5_parametric_bootstrap_results_", output_label, "_summary_table.xlsx"))
}

# Run for all four datasets (no replicate filter needed)

run_daisie_bootstrap_analysis(
  "../data/DAISIE5_sims/DAISIE5_parametric_bootstrap_results_maxage_5.14myr.txt", # input 1000 filtered bootstrap results
  "maxage_5.14myr"
)

run_daisie_bootstrap_analysis(
  "../data/DAISIE5_sims/DAISIE5_parametric_bootstrap_results_maxage_6.15myr.txt", # input 1000 filtered bootstrap results
  "maxage_6.15myr"
)

run_daisie_bootstrap_analysis(
  "../data/DAISIE5_sims/DAISIE5_parametric_bootstrap_results_precise_5.14myr.txt", # input 1000 filtered bootstrap results
  "precise_5.14myr"
)

run_daisie_bootstrap_analysis(
  "../data/DAISIE5_sims/DAISIE5_parametric_bootstrap_results_precise_6.15myr.txt", # input 1000 filtered bootstrap results
  "precise_6.15myr"
)
