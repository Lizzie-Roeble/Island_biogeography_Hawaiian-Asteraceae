#!/usr/bin/env Rscript

# Adjusted Robinson Foulds Distance

# Load necessary libraries
library(phangorn)
library(ape)
library(phytools)
library(ggplot2)
library(vegan)
library(reshape2)
library(viridis)
library(ggrepel)

# Read the trees
# Be sure path is correct
tr1 <- read.tree("../Hybpiper/Species_tree/Species_tree_hybpiper_wastral.tre")
tr2 <- read.tree("../Paragone/Species_tree/PF/Species_tree_paragonePF_wastral.tre")
tr3 <- read.tree("../Paragone/Species_tree/MO/Species_tree_paragoneMO_wastral.tre")
tr4 <- read.tree("../Paralogs/Species_tree/Species_tree_paralogs_astralpro.tre")

# 1. Preparation 

## Root the trees
outgroup <- c("Athrixia_phylicoides_R14373", "Phagnalon_sordidum_GC2743", "Leysera_gnaphalodes_R14525", "Oedera_pungens_R14489")

tr1 <- root(tr1, outgroup = outgroup, resolve.root = TRUE)
tr2 <- root(tr2, outgroup = outgroup, resolve.root = TRUE)
tr3 <- root(tr3, outgroup = outgroup, resolve.root = TRUE)
tr4 <- root(tr4, outgroup = outgroup, resolve.root = TRUE)

# Check if the trees are rooted
cat("Is tr1 rooted?", is.rooted(tr1), "\n")
cat("Is tr2 rooted?", is.rooted(tr2), "\n")
cat("Is tr3 rooted?", is.rooted(tr3), "\n")
cat("Is tr4 rooted?", is.rooted(tr4), "\n")

# 2. Make a table of RFadj values

# Create an empty matrix with 4 rows and 4 columns
matrixinp = matrix(data=NA, nrow=4, ncol=4) 

# Display the empty matrix 
print(matrixinp)

# Make list of trees
tree_list <- list(tr1, tr2, tr3, tr4)

for(j in 1:length(tree_list)){ 
  for (i in 1:length(tree_list)) {
    tree1 <- tree_list[[i]]  
    tree2 <- tree_list[[j]]
    matrixinp[i,j] = RF.dist(tree1, tree2, check.labels = TRUE, normalize = TRUE, rooted = TRUE)
  }
} 

# Give the matrix column and row names
colnames(matrixinp) <- c("HybPiper", "ParaGone-PF", "ParaGone-MO", "Paralogs")
rownames(matrixinp) <- c("HybPiper", "ParaGone-PF", "ParaGone-MO", "Paralogs")

# View the matrix
matrixinp

# export the matrix as a csv
write.csv(matrixinp, "RFDistance_Gnaphalieae_coalescent.csv")

# Plot the heatmap

# Melt the matrix to long format
matrix_long <- melt(matrixinp)

# Filter out the lower triangle and the diagonal by setting lower triangle values to NA
matrix_long <- matrix_long[as.numeric(matrix_long$Var1) <= as.numeric(matrix_long$Var2), ]

# Plot the heatmap with viridis color scale
heatmap_plot <- ggplot(data = matrix_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_viridis(name = "RF Distance", option = "D", direction = 1) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1),
        axis.text.y = element_text(size = 10),
        title = element_text(size = 10)) +
  labs(x = 'Tree Comparison', y = 'Tree Comparison', title = 'Heatmap of Adjusted Robinson-Foulds Distances')

# Print the heatmap
print(heatmap_plot)

ggsave(heatmap_plot,
       filename = "heatmap_RFDistance_Gnaphalieae_coalescent.pdf", 
       device = "pdf", width = 5, height = 4, units = "in")


# 3. Visualize Robinson-Foulds values as a PCA
# Classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis:

pca <- pcoa(matrixinp, correction = "none", rn = NULL)

# Convert PCoA results to a data frame for plotting
pca_df <- as.data.frame(pca$vectors)
pca_df$Sample <- rownames(pca_df)  # Add sample names for labeling

# Plot PCA
pca_plot <- ggplot(pca_df, aes(x = Axis.1, y = Axis.2)) +
  geom_point(size = 3, aes(color = Sample)) +
  geom_text_repel(aes(label = Sample), size = 4, box.padding = 0.5, point.padding = 0.5, max.overlaps = Inf) +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = length(unique(pca_df$Sample)), name = "Set1")) +  # Use a default ggplot2 color palette
  theme_bw() +  # Use theme_bw for a clean background
  theme(plot.margin = margin(10, 10, 10, 10),  # Ensure there's space around the plot
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none") +  # Remove legend
  labs(x = paste("PCoA Axis 1 (", round(pca$values$Relative_eigenvalues[1] * 100, 1), "%)", sep = ""),
       y = paste("PCoA Axis 2 (", round(pca$values$Relative_eigenvalues[2] * 100, 1), "%)", sep = ""),
       title = "PCA of Robinson-Foulds Distances")
pca_plot

# Save the PCA plot
ggsave(pca_plot, filename = "pca_RFDistance_Gnaphalieae_coalescent.pdf", device = "pdf", width = 6, height = 5, units = "in")
