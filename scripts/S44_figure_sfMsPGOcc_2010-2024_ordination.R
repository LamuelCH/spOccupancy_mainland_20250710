library(spOccupancy)
library(ggplot2)
library(reshape2)
library(corrplot)
library(vegan)
library(ggrepel)
library(igraph)
library(gridExtra)
library(grid)  # Added for textGrob
library(viridis)
library(dplyr)

# ============================================================================
# SPECIES RESIDUAL CORRELATION ANALYSIS AND VISUALIZATION
# ============================================================================

cat("============================================================\n")
cat("SPECIES RESIDUAL CORRELATION ANALYSIS\n")
cat("Script started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================\n\n")

# Create output directories
output_dir <- "output/sfMsPGOcc/residuals"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(output_dir, "/data"), recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# MAIN ANALYSIS FUNCTION
# ============================================================================

analyze_species_residual_correlations <- function(model_obj, output_path = output_dir) {
  
  cat("EXTRACTING MODEL INFORMATION:\n")
  cat("-----------------------------\n")
  
  # Extract lambda samples (factor loadings)
  lambda_samples <- model_obj$lambda.samples
  
  # Get dimensions from the column names
  colnames_lambda <- colnames(lambda_samples)
  
  # Extract species names and dimensions
  # Assuming column names are in format "species-factor"
  species_factor_splits <- strsplit(colnames_lambda, "-")
  species_names <- unique(sapply(species_factor_splits, function(x) paste(x[1:(length(x)-1)], collapse="-")))
  n_species <- length(species_names)
  n_factors <- model_obj$q
  n_samples <- dim(lambda_samples)[1]
  
  cat("Model parameters:\n")
  cat("- Number of species:", n_species, "\n")
  cat("- Species names:", paste(species_names, collapse=", "), "\n")
  cat("- Number of factors:", n_factors, "\n")
  cat("- Number of MCMC samples:", n_samples, "\n")
  
  # ============================================================================
  # CALCULATE CORRELATION MATRICES
  # ============================================================================
  
  cat("\nCALCULATING CORRELATION MATRICES:\n")
  cat("---------------------------------\n")
  
  # Initialize array to store correlation matrices for each MCMC sample
  cor_samples <- array(NA, dim = c(n_samples, n_species, n_species))
  
  # Progress bar for correlation calculation
  pb <- txtProgressBar(min = 0, max = n_samples, style = 3)
  
  for (i in 1:n_samples) {
    setTxtProgressBar(pb, i)
    
    # Extract lambda values for this MCMC sample
    lambda_i_vec <- lambda_samples[i,]
    
    # Reshape into matrix [n_species, n_factors]
    lambda_i <- matrix(NA, nrow = n_species, ncol = n_factors)
    
    # Fill the lambda matrix
    for (sp in 1:n_species) {
      for (fac in 1:n_factors) {
        col_name <- paste0(species_names[sp], "-", fac)
        col_idx <- which(colnames_lambda == col_name)
        if (length(col_idx) > 0) {
          lambda_i[sp, fac] <- lambda_i_vec[col_idx]
        }
      }
    }
    
    # Calculate covariance matrix: Lambda * Lambda'
    cov_i <- lambda_i %*% t(lambda_i)
    
    # Convert to correlation matrix
    sd_i <- sqrt(diag(cov_i))
    cor_i <- cov_i / (sd_i %o% sd_i)
    
    # Ensure diagonal is exactly 1
    diag(cor_i) <- 1
    
    # Store the correlation matrix
    cor_samples[i,,] <- cor_i
  }
  
  close(pb)
  
  # ============================================================================
  # CALCULATE SUMMARY STATISTICS
  # ============================================================================
  
  cat("\n\nCALCULATING SUMMARY STATISTICS:\n")
  cat("-------------------------------\n")
  
  # Calculate mean correlation
  mean_cor <- apply(cor_samples, c(2,3), mean)
  rownames(mean_cor) <- species_names
  colnames(mean_cor) <- species_names
  
  # Calculate standard deviation
  sd_cor <- apply(cor_samples, c(2,3), sd)
  
  # Calculate credible intervals
  lower_95 <- apply(cor_samples, c(2,3), function(x) quantile(x, 0.025))
  upper_95 <- apply(cor_samples, c(2,3), function(x) quantile(x, 0.975))
  lower_90 <- apply(cor_samples, c(2,3), function(x) quantile(x, 0.05))
  upper_90 <- apply(cor_samples, c(2,3), function(x) quantile(x, 0.95))
  
  # Identify significant associations
  sig_associations_95 <- (lower_95 > 0) | (upper_95 < 0)
  sig_associations_90 <- (lower_90 > 0) | (upper_90 < 0)
  
  # Create significant correlation matrix
  sig_cor <- mean_cor * sig_associations_95
  
  # Count significant correlations
  n_sig_pos <- sum(sig_associations_95 & mean_cor > 0 & upper.tri(sig_associations_95))
  n_sig_neg <- sum(sig_associations_95 & mean_cor < 0 & upper.tri(sig_associations_95))
  
  cat("Correlation summary:\n")
  cat("- Mean absolute correlation:", round(mean(abs(mean_cor[upper.tri(mean_cor)])), 3), "\n")
  cat("- Max positive correlation:", round(max(mean_cor[upper.tri(mean_cor)]), 3), "\n")
  cat("- Max negative correlation:", round(min(mean_cor[upper.tri(mean_cor)]), 3), "\n")
  cat("- Significant positive correlations (95% CI):", n_sig_pos, "\n")
  cat("- Significant negative correlations (95% CI):", n_sig_neg, "\n")
  
  # ============================================================================
  # SAVE CORRELATION DATA
  # ============================================================================
  
  cat("\nSAVING CORRELATION DATA:\n")
  cat("------------------------\n")
  
  # Save correlation matrices
  save(mean_cor, sd_cor, lower_95, upper_95, sig_associations_95, sig_cor,
       file = paste0(output_path, "/data/correlation_matrices.RData"))
  
  # Export as CSV for reference
  write.csv(mean_cor, file = paste0(output_path, "/data/mean_correlations.csv"))
  write.csv(sig_cor, file = paste0(output_path, "/data/significant_correlations.csv"))
  
  cat("- Saved correlation matrices to:", paste0(output_path, "/data/\n"))
  
  # ============================================================================
  # CREATE VISUALIZATIONS
  # ============================================================================
  
  cat("\nCREATING VISUALIZATIONS:\n")
  cat("------------------------\n")
  
  # 1. CORRELATION HEATMAP
  cat("1. Creating correlation heatmap...\n")
  
  # Prepare data for heatmap
  cor_melt <- melt(mean_cor)
  sig_melt <- melt(sig_associations_95)
  cor_melt$significant <- sig_melt$value
  cor_melt$label <- ifelse(cor_melt$significant & cor_melt$Var1 != cor_melt$Var2, 
                           sprintf("%.2f", cor_melt$value), "")
  
  p_heatmap <- ggplot(cor_melt, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    geom_text(aes(label = label), size = 3, color = "black") +
    scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", 
                         midpoint = 0, limits = c(-1, 1),
                         name = "Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(hjust = 1),
          plot.title = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 12),
          panel.grid = element_blank()) +
    labs(title = "Species Residual Correlations",
         subtitle = "Values shown for significant correlations (95% CI excludes 0)",
         x = "", y = "") +
    coord_fixed()
  
  ggsave(paste0(output_path, "/correlation_heatmap.png"), 
         p_heatmap, width = 10, height = 10, dpi = 300)
  ggsave(paste0(output_path, "/correlation_heatmap.pdf"), 
         p_heatmap, width = 10, height = 10)
  
  # 2. NMDS ORDINATION
  cat("2. Creating NMDS ordination...\n")
  
  # Convert correlation to distance
  dist_matrix <- as.dist(1 - abs(sig_cor))
  
  # Perform NMDS with error handling
  nmds_result <- tryCatch({
    metaMDS(dist_matrix, k = 2, trymax = 100, trace = 0)
  }, error = function(e) {
    cat("   NMDS failed, using PCA instead\n")
    NULL
  })
  
  if (!is.null(nmds_result)) {
    # NMDS succeeded
    species_coords <- as.data.frame(scores(nmds_result))
    colnames(species_coords) <- c("Dim1", "Dim2")
    species_coords$Species <- rownames(species_coords)
    method_name <- "NMDS"
    axis_labels <- list(x = "NMDS1", y = "NMDS2")
    stress_text <- paste0("Stress = ", round(nmds_result$stress, 3))
  } else {
    # Use PCA as fallback
    pca_result <- prcomp(sig_cor, scale = TRUE)
    species_coords <- as.data.frame(pca_result$x[,1:2])
    colnames(species_coords) <- c("Dim1", "Dim2")
    species_coords$Species <- rownames(species_coords)
    method_name <- "PCA"
    
    var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
    axis_labels <- list(
      x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
      y = paste0("PC2 (", round(var_explained[2], 1), "%)")
    )
    stress_text <- ""
  }
  
  # Prepare edge data for significant correlations
  edge_data <- melt(sig_cor)
  names(edge_data) <- c("Species1", "Species2", "Correlation")
  edge_data <- edge_data[edge_data$Species1 != edge_data$Species2 & 
                           abs(edge_data$Correlation) > 0.001, ]
  
  # Add coordinates for edges
  edge_data <- edge_data %>%
    left_join(species_coords, by = c("Species1" = "Species")) %>%
    rename(x = Dim1, y = Dim2) %>%
    left_join(species_coords, by = c("Species2" = "Species")) %>%
    rename(xend = Dim1, yend = Dim2)
  
  # Create ordination plot
  p_ordination <- ggplot(species_coords, aes(x = Dim1, y = Dim2)) +
    # Add edges for correlations
    geom_segment(data = edge_data,
                 aes(x = x, y = y, xend = xend, yend = yend,
                     size = abs(Correlation),
                     color = Correlation),
                 alpha = 0.5) +
    # Add species points
    geom_point(size = 4, color = "black", fill = "white", shape = 21) +
    # Add species labels
    geom_text_repel(aes(label = Species), size = 4, max.overlaps = 20) +
    # Formatting
    scale_size_continuous(range = c(0.5, 3), 
                          name = "Correlation\nStrength",
                          breaks = c(0.3, 0.5, 0.7),
                          labels = c("0.3", "0.5", "0.7")) +
    scale_color_gradient2(low = "#B2182B", mid = "gray80", high = "#2166AC",
                          midpoint = 0, 
                          name = "Correlation\nDirection") +
    theme_minimal() +
    theme(plot.title = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 12),
          legend.position = "right") +
    labs(title = paste0(method_name, " Ordination of Species Residual Correlations"),
         subtitle = paste("Lines show significant correlations (95% CI)", stress_text),
         x = axis_labels$x,
         y = axis_labels$y)
  
  ggsave(paste0(output_path, "/ordination_plot.png"), 
         p_ordination, width = 12, height = 10, dpi = 300)
  ggsave(paste0(output_path, "/ordination_plot.pdf"), 
         p_ordination, width = 12, height = 10)
  
  # 3. NETWORK VISUALIZATION
  cat("3. Creating network visualization...\n")
  
  # Filter edge data for significant correlations above threshold
  edge_data_filtered <- edge_data[abs(edge_data$Correlation) > 0.3, ]
  
  if (nrow(edge_data_filtered) > 0) {
    # Create edge list with just the necessary columns for igraph
    edge_list_for_graph <- edge_data_filtered[, c("Species1", "Species2", "Correlation")]
    
    # Remove any duplicate edges (keeping only unique pairs)
    edge_list_for_graph <- edge_list_for_graph[!duplicated(t(apply(edge_list_for_graph[,1:2], 1, sort))), ]
    
    # Create graph from edge list only (no vertices parameter)
    g <- graph_from_data_frame(edge_list_for_graph, directed = FALSE)
    
    # Get all unique species from the graph
    all_species_in_graph <- V(g)$name
    
    cat("   Network includes", length(all_species_in_graph), "species with", 
        ecount(g), "edges\n")
    
    # Set edge attributes
    E(g)$width <- abs(E(g)$Correlation) * 5
    E(g)$color <- ifelse(E(g)$Correlation > 0, "#2166AC", "#B2182B")
    
    # Calculate network metrics
    degree_centrality <- degree(g)
    betweenness_centrality <- betweenness(g)
    
    # Create layout
    set.seed(123) # For reproducibility
    layout <- layout_with_fr(g)
    layout_df <- as.data.frame(layout)
    names(layout_df) <- c("x", "y")
    layout_df$Species <- V(g)$name
    layout_df$Degree <- degree_centrality[V(g)$name]
    layout_df$Betweenness <- betweenness_centrality[V(g)$name]
    
    # Prepare edge data for plotting using the layout coordinates
    edge_plot_data <- edge_data_filtered
    edge_plot_data$x <- layout_df$x[match(edge_plot_data$Species1, layout_df$Species)]
    edge_plot_data$y <- layout_df$y[match(edge_plot_data$Species1, layout_df$Species)]
    edge_plot_data$xend <- layout_df$x[match(edge_plot_data$Species2, layout_df$Species)]
    edge_plot_data$yend <- layout_df$y[match(edge_plot_data$Species2, layout_df$Species)]
    
    # Create network plot
    p_network <- ggplot() +
      # Add edges
      geom_segment(data = edge_plot_data,
                   aes(x = x, y = y, xend = xend, yend = yend,
                       size = abs(Correlation),
                       color = ifelse(Correlation > 0, "Positive", "Negative")),
                   alpha = 0.7) +
      # Add nodes
      geom_point(data = layout_df, 
                 aes(x = x, y = y, size = Degree),
                 color = "black", fill = "white", shape = 21) +
      # Add labels
      geom_text_repel(data = layout_df,
                      aes(x = x, y = y, label = Species),
                      size = 4, max.overlaps = 20) +
      # Formatting
      scale_size_continuous(name = "Node Degree", range = c(3, 10)) +
      scale_color_manual(values = c("Positive" = "#2166AC", "Negative" = "#B2182B"),
                         name = "Correlation") +
      theme_void() +
      theme(plot.title = element_text(size = 16, face = "bold"),
            plot.subtitle = element_text(size = 12),
            legend.position = "right") +
      labs(title = "Network of Species Residual Correlations",
           subtitle = paste("Showing correlations with |r| > 0.3 (", ecount(g), " edges)"))
    
    ggsave(paste0(output_path, "/network_plot.png"), 
           p_network, width = 12, height = 10, dpi = 300)
    ggsave(paste0(output_path, "/network_plot.pdf"), 
           p_network, width = 12, height = 10)
    
    # Save network metrics
    network_metrics <- data.frame(
      Species = layout_df$Species,
      Degree = layout_df$Degree,
      Betweenness = layout_df$Betweenness
    )
    write.csv(network_metrics, 
              file = paste0(output_path, "/data/network_metrics.csv"),
              row.names = FALSE)
  } else {
    cat("   No correlations above threshold (|r| > 0.3) for network visualization\n")
    p_network <- NULL
  }
  
  # 4. CORRELATION DISTRIBUTION PLOT
  cat("4. Creating correlation distribution plot...\n")
  
  # Extract upper triangle correlations
  upper_tri_cor <- mean_cor[upper.tri(mean_cor)]
  upper_tri_sig <- sig_associations_95[upper.tri(sig_associations_95)]
  
  cor_dist_data <- data.frame(
    Correlation = upper_tri_cor,
    Significant = upper_tri_sig
  )
  
  p_distribution <- ggplot(cor_dist_data, aes(x = Correlation, fill = Significant)) +
    geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "#2166AC"),
                      labels = c("Non-significant", "Significant")) +
    theme_minimal() +
    theme(plot.title = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 12)) +
    labs(title = "Distribution of Species Residual Correlations",
         subtitle = "Shaded by significance (95% CI)",
         x = "Correlation Coefficient",
         y = "Count",
         fill = "")
  
  ggsave(paste0(output_path, "/correlation_distribution.png"), 
         p_distribution, width = 10, height = 6, dpi = 300)
  ggsave(paste0(output_path, "/correlation_distribution.pdf"), 
         p_distribution, width = 10, height = 6)
  
  # 5. COMBINED FIGURE FOR PUBLICATION
  cat("5. Creating combined publication figure...\n")
  
  # Create a combined figure with key plots
  combined_plot <- gridExtra::grid.arrange(
    p_heatmap, p_ordination,
    ncol = 2,
    top = grid::textGrob("Species Residual Correlation Analysis", 
                         gp = grid::gpar(fontsize = 18, fontface = "bold"))
  )
  
  ggsave(paste0(output_path, "/fig4_combined_residual_analysis.png"), 
         combined_plot, width = 20, height = 10, dpi = 300)
  ggsave(paste0(output_path, "/fig4_combined_residual_analysis.pdf"), 
         combined_plot, width = 20, height = 10)
  
  # ============================================================================
  # CREATE SUMMARY REPORT
  # ============================================================================
  
  cat("\nCREATING SUMMARY REPORT:\n")
  cat("------------------------\n")
  
  # Find strongest correlations
  cor_df <- melt(mean_cor)
  cor_df <- cor_df[cor_df$Var1 != cor_df$Var2, ]
  cor_df <- cor_df[order(abs(cor_df$value), decreasing = TRUE), ]
  
  # Get significant ones
  sig_df <- melt(sig_cor)
  sig_df <- sig_df[sig_df$Var1 != sig_df$Var2 & sig_df$value != 0, ]
  sig_df <- sig_df[order(abs(sig_df$value), decreasing = TRUE), ]
  
  # Write summary report
  report_file <- paste0(output_path, "/correlation_analysis_summary.txt")
  sink(report_file)
  
  cat("SPECIES RESIDUAL CORRELATION ANALYSIS SUMMARY\n")
  cat("=============================================\n")
  cat("Analysis date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  cat("Model Information:\n")
  cat("- Number of species:", n_species, "\n")
  cat("- Number of factors:", n_factors, "\n")
  cat("- Number of MCMC samples:", n_samples, "\n\n")
  
  cat("Correlation Summary:\n")
  cat("- Total pairwise correlations:", length(upper_tri_cor), "\n")
  cat("- Significant correlations (95% CI):", sum(upper_tri_sig), "\n")
  cat("- Proportion significant:", round(sum(upper_tri_sig)/length(upper_tri_sig)*100, 1), "%\n")
  cat("- Mean absolute correlation:", round(mean(abs(upper_tri_cor)), 3), "\n")
  cat("- Correlation range: [", round(min(upper_tri_cor), 3), ", ", 
      round(max(upper_tri_cor), 3), "]\n\n")
  
  cat("Top 10 Strongest Correlations:\n")
  cat("-----------------------------\n")
  for (i in 1:min(10, nrow(sig_df))) {
    cat(sprintf("%2d. %s - %s: %.3f\n", 
                i, sig_df$Var1[i], sig_df$Var2[i], sig_df$value[i]))
  }
  
  cat("\n\nOutput Files:\n")
  cat("-------------\n")
  cat("- Correlation heatmap: correlation_heatmap.png/pdf\n")
  cat("- Ordination plot: ordination_plot.png/pdf\n")
  cat("- Network plot: network_plot.png/pdf\n")
  cat("- Distribution plot: correlation_distribution.png/pdf\n")
  cat("- Combined figure: fig4_combined_residual_analysis.png/pdf\n")
  cat("- Data files: data/correlation_matrices.RData\n")
  cat("- CSV exports: data/mean_correlations.csv, data/significant_correlations.csv\n")
  
  sink()
  
  cat("Summary report saved to:", report_file, "\n")
  
  # ============================================================================
  # RETURN RESULTS
  # ============================================================================
  
  results <- list(
    correlation_matrix = mean_cor,
    significant_correlations = sig_cor,
    credible_intervals = list(lower_95 = lower_95, upper_95 = upper_95),
    species_coordinates = species_coords,
    plots = list(
      heatmap = p_heatmap,
      ordination = p_ordination,
      network = if(exists("p_network")) p_network else NULL,
      distribution = p_distribution
    ),
    summary_stats = list(
      n_species = n_species,
      n_factors = n_factors,
      n_significant = sum(upper_tri_sig),
      prop_significant = sum(upper_tri_sig)/length(upper_tri_sig)
    )
  )
  
  cat("\nAnalysis completed successfully!\n")
  cat("Results saved to:", output_path, "\n")
  
  return(results)
}

# ============================================================================
# EXAMPLE USAGE
# ============================================================================

# Example 1: If you have the model object already loaded
# results <- analyze_species_residual_correlations(out.sfMsPGOcc)

# Example 2: If you need to load the model from file
load("models/sfMsPGOcc/model_20250718_sfMsPGOcc_2010-2024_nthin500_nburn5e+05_nchain1_nsample6000_nfactors3_neighbors10.RData")
results <- analyze_species_residual_correlations(out.sfMsPGOcc)

# Example 3: Custom output directory
# results <- analyze_species_residual_correlations(out.sfMsPGOcc, 
#                                                  output_path = "output/custom_path")

cat("\n============================================================\n")
cat("Script ready. Use analyze_species_residual_correlations() function\n")
cat("with your fitted model object to generate all analyses.\n")
cat("============================================================\n")