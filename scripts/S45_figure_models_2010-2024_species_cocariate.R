# Load necessary libraries
library(dplyr)
library(ggplot2)
library(spOccupancy)

# Load all four models
# Adjust file paths as needed
load("models/PGOcc/model_20250710_PGOcc_2010-2024_5r_nthin100_nburn1e+05_nchain4_nsample250000.RData")          # Single species quoll
load("models/msPGOcc/model_20250710_msPGOcc_2010-2024_nthin100_nburn1e+05_nchain4_nsample250000.RData")      # Multi-species
load("models/spPGOcc/model_20250710_spPGOcc_2010-2024_nthin1000_nburn1e+06_nchain1_nsample6000_neighbors10.RData")      # Single species spatial quoll  
load("models/sfMsPGOcc/model_20250718_sfMsPGOcc_2010-2024_nthin500_nburn5e+05_nchain1_nsample6000_nfactors3_neighbors10.RData")  # Multi-species spatial

# Function to extract quoll-specific effects from different model types
extract_quoll_effects <- function(model_obj, model_type, target_species = "Dasyurus") {
  
  process_single_species_samples <- function(samples, param_type) {
    if (is.null(samples) || ncol(samples) == 0) {
      message(paste("No samples provided for param_type:", param_type))
      return(NULL)
    }
    
    all_param_names <- colnames(samples)
    results_list <- list()
    
    message(paste("\n--- Processing single species samples for type:", param_type, "---"))
    message("Parameter names found:")
    print(all_param_names)
    
    # Process each parameter
    for (param_name in all_param_names) {
      message(paste("  Processing parameter:", param_name))
      param_samples <- samples[, param_name]
      
      results_list[[length(results_list) + 1]] <- data.frame(
        Parameter = param_name,
        Mean = mean(param_samples),
        Lower = quantile(param_samples, 0.025, type = 8),
        Upper = quantile(param_samples, 0.975, type = 8),
        Type = param_type,
        Model = model_type
      )
    }
    
    if (length(results_list) == 0) return(NULL)
    return(do.call(rbind, results_list))
  }
  
  process_multi_species_samples <- function(samples, param_type, species_name) {
    if (is.null(samples) || ncol(samples) == 0) {
      message(paste("No samples provided for param_type:", param_type))
      return(NULL)
    }
    
    all_param_names <- colnames(samples)
    species_params <- all_param_names[grepl(paste0("-", species_name, "$"), all_param_names)]
    
    message(paste("\n--- Processing", species_name, "samples for type:", param_type, "---"))
    message(paste("Species-specific parameters found for", species_name, ":"))
    print(species_params)
    
    if (length(species_params) == 0) {
      message(paste("No parameters found for species:", species_name))
      return(NULL)
    }
    
    results_list <- list()
    
    # Process each species-specific parameter
    for (param_name in species_params) {
      # Remove species suffix to get base parameter name
      base_param <- gsub(paste0("-", species_name, "$"), "", param_name)
      message(paste("  Processing species parameter:", param_name, "->", base_param))
      param_samples <- samples[, param_name]
      
      results_list[[length(results_list) + 1]] <- data.frame(
        Parameter = base_param,
        Mean = mean(param_samples),
        Lower = quantile(param_samples, 0.025, type = 8),
        Upper = quantile(param_samples, 0.975, type = 8),
        Type = param_type,
        Model = model_type
      )
    }
    
    if (length(results_list) == 0) return(NULL)
    return(do.call(rbind, results_list))
  }
  
  # Extract effects based on model type
  if (model_type %in% c("PGOcc", "spPGOcc")) {
    # Single species models - extract directly from beta and alpha samples
    beta_df <- process_single_species_samples(model_obj$beta.samples, "Occurrence")
    alpha_df <- process_single_species_samples(model_obj$alpha.samples, "Detection")
  } else {
    # Multi-species models - extract species-specific effects
    beta_df <- process_multi_species_samples(model_obj$beta.samples, "Occurrence", target_species)
    alpha_df <- process_multi_species_samples(model_obj$alpha.samples, "Detection", target_species)
  }
  
  # Combine results
  result <- rbind(beta_df, alpha_df)
  
  if (is.null(result) || nrow(result) == 0) {
    warning(paste("No parameters were processed for model type:", model_type))
    return(data.frame())
  }
  
  return(result)
}

# Extract effects from all four models
message("Extracting effects from PGOcc model...")
df_PGOcc <- extract_quoll_effects(out.PGOcc, "PGOcc")

message("Extracting effects from msPGOcc model...")
df_msPGOcc <- extract_quoll_effects(out.msPGOcc, "msPGOcc")

message("Extracting effects from spPGOcc model...")
df_spPGOcc <- extract_quoll_effects(out.spPGOcc, "spPGOcc")

message("Extracting effects from sfMsPGOcc model...")
df_sfMsPGOcc <- extract_quoll_effects(out.sfMsPGOcc, "sfMsPGOcc")

# Combine all results
df_all_models <- rbind(df_PGOcc, df_msPGOcc, df_spPGOcc, df_sfMsPGOcc)

# Check if data was generated
if (nrow(df_all_models) > 0) {
  
  # Create custom y-axis labels
  df_all_models$ParameterDisplay <- sapply(df_all_models$Parameter, function(param_name) {
    # Custom label mapping
    custom_labels <- c(
      "(Intercept)" = "Baseline",
      "scale(effort)" = "Survey Effort",
      "scale(bio5)" = "Max Temperature",
      "scale(bio6)" = "Min Temperature", 
      "scale(bio15)" = "Precipitation Seasonality",
      "scale(dem)" = "Elevation",
      "scale(tpi)" = "Topographic Position Index",
      "scale(fpc)" = "Forest Cover"
    )
    
    # Check if we have a custom label
    if (param_name %in% names(custom_labels)) {
      return(custom_labels[[param_name]])
    }
    
    # For parameters not in custom list, apply the original logic
    if (param_name == "(Intercept)") {
      return("Intercept")
    }
    
    # Check if it's an explicit quadratic term I(...^2)
    is_explicit_quadratic <- grepl("I\\(.*\\^2\\)", param_name)
    
    # Create a base name: remove scale(), I(), parentheses, and ^2 for cleaner display
    base_name <- param_name
    base_name <- gsub("scale\\(", "", base_name) # Remove "scale("
    base_name <- gsub("I\\(", "", base_name)     # Remove "I("
    base_name <- gsub("\\^2", "", base_name)     # Remove "^2"
    base_name <- gsub("\\)", "", base_name)     # Remove all remaining parentheses
    base_name <- trimws(base_name)              # Clean up any leading/trailing spaces
    
    if (is_explicit_quadratic) {
      return(paste0(base_name, "²")) # Append a superscript 2
    } else {
      return(base_name)
    }
  })
  
  # Define the IsQuadratic column for the plot's shape aesthetic
  df_all_models$IsQuadratic <- grepl("I\\(.*\\^2\\)", df_all_models$Parameter)
  
  # Add formatting columns
  df_all_models$MeanLabel <- format(round(df_all_models$Mean, 2), nsmall = 2)
  df_all_models$CrossesZero <- (df_all_models$Lower < 0 & df_all_models$Upper > 0)
  
  # Create parameter groupings
  df_all_models <- df_all_models %>%
    mutate(Group = case_when(
      Parameter == "(Intercept)" ~ "Baseline",
      grepl("effort", Parameter) ~ "Effort",
      grepl("bio5|bio6|bio15", Parameter) ~ "Climate",
      grepl("dem|tpi", Parameter) ~ "Terrain", 
      grepl("fpc", Parameter) ~ "Forest",
      TRUE ~ "Other"
    ))
  
  # Set factor levels to control order (fixed the typo "Basline" -> "Baseline")
  df_all_models$Group <- factor(df_all_models$Group, 
                                levels = c("Baseline", "Effort", "Climate", "Terrain", "Forest", "Other"))
  
  # Set model factor levels for consistent ordering
  df_all_models$Model <- factor(df_all_models$Model, 
                                levels = c("PGOcc", "msPGOcc", "spPGOcc", "sfMsPGOcc"))
  
  # OPTION: Set to TRUE to exclude quadratic terms from the plot
  exclude_quadratic <- TRUE  # Change to FALSE to include quadratic terms
  
  # Filter data based on quadratic option
  plot_data <- if (exclude_quadratic) {
    df_all_models[!df_all_models$IsQuadratic, ]
  } else {
    df_all_models
  }
  
  # Define colors for different models
  model_colors <- c("PGOcc" = "#E69F00",        # Blue
                    "msPGOcc" = "#F0E442",      # Orange  
                    "spPGOcc" = "#56B4E9",      # Green
                    "sfMsPGOcc" = "#009E73")    # Red
  
  # Create the forest plot comparing all four models
  p_model_comparison <- ggplot(plot_data, 
                               aes(x = Mean, y = ParameterDisplay, color = Model, shape = IsQuadratic)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(aes(xmin = Lower, xmax = Upper, alpha = CrossesZero),
                    position = position_dodge(width = 0.8), size = 0.7) +
    scale_alpha_manual(values = c(`TRUE` = 0.5, `FALSE` = 1), guide = "none") +
    facet_grid(Group ~ Type, scales = "free_y", space = "free_y") +
    scale_color_manual(values = model_colors) +
    scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16)) +
    labs(
      title = "Spotted-tailed Quoll Covariate Effects: Model Comparison",
      subtitle = if (exclude_quadratic) "Occupancy Model Results (Linear Terms Only)" else "Occupancy Model Results",
      x = "Effect Size (logit scale)",
      y = "",
      color = "Model Type",
      shape = if (exclude_quadratic) NULL else "Quadratic Term"
    ) +
    theme_bw() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(hjust = 0),
      # Transparency settings
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.box.background = element_rect(fill = "transparent", color = NA)
    )
  
  # Calculate dynamic dimensions
  n_parameters <- length(unique(plot_data$ParameterDisplay))
  n_types <- length(unique(plot_data$Type))
  plot_height <- max(8, n_parameters * 0.4 + 3)
  plot_width <- max(10, n_types * 5)
  
  # Define filename
  filename <- paste0("output/plots/quoll_model_comparison_", 
                     if (exclude_quadratic) "linear_only" else "all_terms", 
                     ".png")
  
  # Save with dynamic sizing and transparent background
  ggsave(
    filename = filename,
    plot = p_model_comparison,
    width = plot_width,
    height = plot_height,
    dpi = 300,
    bg = "transparent",
    limitsize = FALSE
  )
  
  # Print summary statistics
  cat("\n=== MODEL COMPARISON SUMMARY ===\n")
  cat("Number of parameters per model:\n")
  print(table(plot_data$Model, plot_data$Type))
  
  cat("\nParameters that cross zero by model:\n")
  crossing_zero <- plot_data %>%
    group_by(Model, Type) %>%
    summarise(
      total_params = n(),
      crossing_zero = sum(CrossesZero),
      prop_crossing_zero = round(crossing_zero / total_params, 2),
      .groups = 'drop'
    )
  print(crossing_zero)
  
  # Display the plot
  print(p_model_comparison)
  
} else {
  warning("No data was extracted from any models. Check model objects and species name.")
}
