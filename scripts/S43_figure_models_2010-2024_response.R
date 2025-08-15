library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
library(cowplot) # Added for get_legend()

# =============================================================================
# CONFIGURATION SECTION
# =============================================================================

# Define a color-blind-safe palette for better accessibility
okabe_ito_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")

# Define models to compare using the new palette
models_config <- list(
  list(
    name = "PGOcc",
    file = "output/PGOcc/predictions/marginal/PGOcc_2010-2024_response_curves_results.RData",
    type = "single",
    species = NULL,
    color = okabe_ito_palette[1]
  ),
  list(
    name = "spPGOcc",
    file = "output/spPGOcc/predictions/marginal/spPGOcc_2010-2024_response_curves_results.RData",
    type = "single",
    species = NULL,
    color = okabe_ito_palette[2]
  ),
  list(
    name = "sfMsPGOcc",
    file = "output/sfMsPGOcc/predictions/marginal/sfMsPGOcc_2010-2024_response_curves_results.RData",
    type = "multi",
    species = "Dasyurus",
    color = okabe_ito_palette[3]
  ),
  list(
    name = "msPGOcc",
    file = "output/msPGOcc/predictions/marginal/msPGOcc_2010-2024_response_curves_results.RData",
    type = "multi",
    species = "Dasyurus",
    color = okabe_ito_palette[4]
  )
)

# Define covariate labels for better plot annotation
covariate_labels <- c(
  "bio5" = "Maximum Temperature (°C)",
  "bio6" = "Minimum Temperature (°C)",
  "bio15" = "Precipitation Seasonality (%)",
  "dem" = "Elevation (m)",
  "tpi" = "Topographic Position Index",
  "fpc" = "Foliage Projective Cover (%)"
)

# Create a dedicated output directory for publication plots
output_dir <- "output/plots/figures"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# DATA LOADING AND PROCESSING (No changes made to this section)
# =============================================================================
cat("=== LOADING MODEL RESULTS ===\n")
extract_plot_data <- function(model_config) {
  cat(sprintf("Loading %s from %s\n", model_config$name, basename(model_config$file)))
  if (!file.exists(model_config$file)) {
    warning(sprintf("File not found: %s", model_config$file))
    return(NULL)
  }
  temp_env <- new.env()
  load(model_config$file, envir = temp_env)
  results <- if (exists("results_list", envir = temp_env)) {
    temp_env$results_list
  } else if (exists("results", envir = temp_env)) {
    temp_env$results
  } else {
    objs <- ls(envir = temp_env)
    res <- NULL
    for (obj in objs) {
      if (is.list(temp_env[[obj]])) {
        res <- temp_env[[obj]]
        break
      }
    }
    res
  }
  if (is.null(results)) {
    warning(sprintf("Could not find results in %s", model_config$file))
    return(NULL)
  }
  all_data <- list()
  for (cov_name in names(results)) {
    if (!results[[cov_name]]$success) {
      cat(sprintf("  Skipping %s (failed)\n", cov_name))
      next
    }
    if (!is.null(results[[cov_name]]$plot_data)) {
      plot_data <- results[[cov_name]]$plot_data
    } else {
      pred_obj <- results[[cov_name]]$prediction
      if (is.null(pred_obj)) next
      cov_values <- pred_obj$covariate_values
      if (model_config$type == "single") {
        psi_samples <- pred_obj$psi.0.samples
        occ_mean <- apply(psi_samples, 2, mean)
        occ_lower <- apply(psi_samples, 2, function(x) quantile(x, 0.025))
        occ_upper <- apply(psi_samples, 2, function(x) quantile(x, 0.975))
        plot_data <- data.frame(
          covariate = cov_values, occ_prob = occ_mean,
          lower = occ_lower, upper = occ_upper, sp = "Target_Species"
        )
      } else {
        psi_samples <- pred_obj$psi.0.samples
        if (is.character(model_config$species)) {
          sp_names <- pred_obj$sp.names
          if(is.null(sp_names) && !is.null(results[[cov_name]]$summary_stats)) {
            sp_names <- names(results[[cov_name]]$summary_stats)
          }
          sp_idx <- which(sp_names == model_config$species)
          if (length(sp_idx) == 0) {
            warning(sprintf("Species %s not found in %s", model_config$species, model_config$name))
            next
          }
        } else {
          sp_idx <- model_config$species
        }
        if (length(dim(psi_samples)) == 3) {
          sp_psi <- psi_samples[, sp_idx, ]
        } else {
          warning(sprintf("Unexpected dimension for multi-species model %s", model_config$name))
          next
        }
        occ_mean <- apply(sp_psi, 2, mean)
        occ_lower <- apply(sp_psi, 2, function(x) quantile(x, 0.025))
        occ_upper <- apply(sp_psi, 2, function(x) quantile(x, 0.975))
        plot_data <- data.frame(
          covariate = cov_values, occ_prob = occ_mean,
          lower = occ_lower, upper = occ_upper, sp = model_config$species
        )
      }
    }
    plot_data$model <- model_config$name
    all_data[[cov_name]] <- plot_data
  }
  return(all_data)
}
all_model_data <- list()
for (i in seq_along(models_config)) {
  model_data <- extract_plot_data(models_config[[i]])
  if (!is.null(model_data)) {
    all_model_data[[models_config[[i]]$name]] <- model_data
  }
}
cat(sprintf("\nSuccessfully loaded data from %d models\n", length(all_model_data)))
all_covariates <- Reduce(intersect, lapply(all_model_data, names))
cat(sprintf("Common covariates across all models: %s\n", paste(all_covariates, collapse = ", ")))

# =============================================================================
# REFINED PLOTTING FUNCTION
# =============================================================================
create_comparison_plot <- function(covariate_name, model_data_list, model_configs) {
  combined_data <- bind_rows(lapply(names(model_data_list), function(m) model_data_list[[m]][[covariate_name]]))
  if (is.null(combined_data) || nrow(combined_data) == 0) {
    warning(sprintf("No data for covariate %s", covariate_name))
    return(NULL)
  }
  model_colors <- setNames(sapply(model_configs, `[[`, "color"), sapply(model_configs, `[[`, "name"))
  x_label <- covariate_labels[covariate_name]
  if (is.na(x_label)) x_label <- covariate_name
  
  # Create comparison plot with publication-ready theme
  p <- ggplot(combined_data, aes(x = covariate, y = occ_prob, color = model, fill = model)) +
    geom_line(linewidth = 1.0) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = 0) +
    # --- FIX: Use guides() to properly merge legends ---
    scale_color_manual(values = model_colors) +
    scale_fill_manual(values = model_colors) +
    guides(color = guide_legend("Model"), fill = guide_legend("Model")) +
    scale_y_continuous(labels = percent_format(), limits = c(0, NA)) +
    labs(
      x = x_label, 
      y = "Occurrence Probability"
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 11, color = "black"),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    )
  return(p)
}

# =============================================================================
# CREATE COMPARISON PLOTS
# =============================================================================
cat("\n=== CREATING PUBLICATION-QUALITY PLOTS ===\n")
comparison_plots <- list()
for (cov in all_covariates) {
  cat(sprintf("Creating plot for %s\n", cov))
  p <- create_comparison_plot(cov, all_model_data, models_config)
  if (!is.null(p)) {
    comparison_plots[[cov]] <- p
    # Save individual plots in multiple high-quality formats
    ggsave(
      file.path(output_dir, paste0("comparison_", cov, ".pdf")), p,
      width = 7, height = 5, device = cairo_pdf
    )
    ggsave(
      file.path(output_dir, paste0("comparison_", cov, ".tiff")), p,
      width = 7, height = 5, dpi = 600
    )
  }
}

# =============================================================================
# CREATE COMBINED PLOTS - FIXED VERSION
# =============================================================================
cat("\n=== CREATING COMBINED PLOTS ===\n")
create_plot_grid <- function(plot_list, ncol = 3, title = NULL, subtitle = NULL, legend_position = "bottom") {
  # Remove legends from all plots
  plots_no_legend <- lapply(plot_list, function(p) p + theme(legend.position = "none"))
  
  # Extract legend from the first plot using a more robust approach
  # Create a simplified version of the first plot just for legend extraction
  first_plot_for_legend <- plot_list[[1]] + 
    theme(
      legend.position = legend_position,
      legend.box.margin = margin(t = 10, b = 5),
      legend.box = "horizontal",  # Ensure horizontal layout
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      legend.key.size = unit(1.2, "cm"),  # Make legend keys larger
      legend.spacing.x = unit(0.5, "cm")   # Add spacing between legend items
    )
  
  # Extract legend using cowplot's approach but with error handling
  legend <- tryCatch({
    get_legend(first_plot_for_legend)
  }, error = function(e) {
    # Fallback: create a simpler legend extraction
    warning("Standard legend extraction failed, using alternative method")
    # Alternative approach using ggplotGrob
    g <- ggplotGrob(first_plot_for_legend)
    legend_idx <- which(g$layout$name == "guide-box")
    if(length(legend_idx) > 0) {
      legend_grob <- g$grobs[[legend_idx[1]]]
      return(legend_grob)
    } else {
      warning("Could not extract legend")
      return(NULL)
    }
  })
  
  # Create the plot grid
  plot_grid <- wrap_plots(plots_no_legend, ncol = ncol)
  
  # Combine with legend (only if legend extraction was successful)
  if (!is.null(legend)) {
    if (legend_position == "bottom") {
      final_plot <- plot_grid / legend + plot_layout(heights = c(1, 0.15))  # Increased legend space
    } else if (legend_position == "right") {
      final_plot <- plot_grid | legend + plot_layout(widths = c(1, 0.2))
    } else {
      final_plot <- plot_grid / legend + plot_layout(heights = c(1, 0.15))
    }
  } else {
    # If legend extraction failed, create a manual legend
    warning("Creating manual legend as fallback")
    final_plot <- create_manual_legend_plot(plot_grid, plot_list)
  }
  
  # Add title and subtitle if provided
  if (!is.null(title) || !is.null(subtitle)) {
    final_plot <- final_plot +
      plot_annotation(
        title = title,
        subtitle = subtitle,
        theme = theme(
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 14, hjust = 0.5)
        )
      )
  }
  return(final_plot)
}

# Helper function to create manual legend if extraction fails
create_manual_legend_plot <- function(plot_grid, plot_list) {
  # Extract model names and colors from the first plot
  first_plot_data <- layer_data(plot_list[[1]])
  if ("colour" %in% names(first_plot_data)) {
    # Create a simple legend plot
    legend_data <- data.frame(
      model = unique(first_plot_data$colour),
      x = 1:length(unique(first_plot_data$colour)),
      y = 1
    )
    
    legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = model)) +
      geom_point(size = 4) +
      scale_color_identity() +
      theme_void() +
      theme(
        legend.position = "bottom",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12)
      ) +
      labs(color = "Model")
    
    return(plot_grid / get_legend(legend_plot) + plot_layout(heights = c(1, 0.15)))
  }
  
  # If all else fails, return just the plot grid
  return(plot_grid)
}

if (length(comparison_plots) > 0) {
  # Add panel labels (a, b, c...) to each plot for the grid
  plots_for_grid <- mapply(function(p, label) {
    p + labs(tag = label) + theme(plot.tag = element_text(face = "bold", size=16))
  }, comparison_plots, paste0("(", letters[1:length(comparison_plots)], ")"), SIMPLIFY = FALSE)
  
  full_grid <- create_plot_grid(
    plots_for_grid,
    ncol = 3,
    title = "Environmental Response Curves: Multi-Model Comparison",
    subtitle = "Spotted-tailed Quoll (Dasyurus maculatus)"
  )
  n_rows <- ceiling(length(comparison_plots) / 3)
  plot_height <- 4 * n_rows + 2
  ggsave(
    file.path(output_dir, "all_response_curves_comparison.pdf"),
    full_grid, width = 15, height = plot_height, device = cairo_pdf
  )
  ggsave(
    file.path(output_dir, "all_response_curves_comparison.tiff"),
    full_grid, width = 15, height = plot_height, dpi = 600
  )
}

# =============================================================================
# MODEL COMPARISON SUMMARY (No changes made to this section)
# =============================================================================
cat("\n=== CREATING MODEL COMPARISON SUMMARY ===\n")
calculate_model_divergence <- function(covariate_name, model_data_list) {
  model_names <- names(model_data_list)
  n_models <- length(model_names)
  if (n_models < 2) return(NULL)
  predictions <- list()
  for (model in model_names) {
    if (covariate_name %in% names(model_data_list[[model]])) {
      predictions[[model]] <- model_data_list[[model]][[covariate_name]]$occ_prob
    }
  }
  if (length(predictions) < 2) return(NULL)
  divergence_matrix <- matrix(NA, n_models, n_models, dimnames = list(model_names, model_names))
  for (i in 1:(n_models - 1)) {
    for (j in (i + 1):n_models) {
      rmsd <- sqrt(mean((predictions[[i]] - predictions[[j]])^2))
      divergence_matrix[i, j] <- rmsd
      divergence_matrix[j, i] <- rmsd
    }
  }
  diag(divergence_matrix) <- 0
  return(divergence_matrix)
}
divergence_summary <- list()
for (cov in all_covariates) {
  div_matrix <- calculate_model_divergence(cov, all_model_data)
  if (!is.null(div_matrix)) {
    divergence_summary[[cov]] <- div_matrix
  }
}
save(divergence_summary, file = file.path(output_dir, "model_divergence_summary.RData"))
cat("\nModel Divergence Summary (RMSD):\n")
for (cov in names(divergence_summary)) {
  cat(sprintf("\n%s:\n", covariate_labels[cov]))
  print(round(divergence_summary[[cov]], 3))
}

cat(sprintf("\n=== COMPLETED ===\n"))
cat(sprintf("Publication-quality plots saved to: %s\n", output_dir))