# =============================================================================
# Posterior Predictive Checks for spOccupancy Models
# Enhanced for four model types: sfMsPGOcc, spPGOcc, msPGOcc, PGOcc
# =============================================================================

# Set up logging
log_file <- file.path("logs", paste0("ppc_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)

# Function to log messages both to console and log file
log_message <- function(msg) {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  full_msg <- paste(timestamp, msg)
  cat(full_msg, "\n")
  cat(full_msg, "\n", file = log_file, append = TRUE)
}

# Function to handle errors
handle_error <- function(expr, step_name) {
  tryCatch(
    expr,
    error = function(e) {
      log_message(paste("ERROR in", step_name, ":", e$message))
      return(NULL)
    },
    warning = function(w) {
      log_message(paste("WARNING in", step_name, ":", w$message))
      # Continue execution despite warning
      invokeRestart("muffleWarning")
    }
  )
}

# Timer function
start_timer <- function() {
  return(Sys.time())
}

end_timer <- function(start_time, step_name) {
  duration <- Sys.time() - start_time
  log_message(paste(step_name, "completed in", round(as.numeric(duration), 2), attr(duration, "units")))
}

# Main execution
log_message("=== Starting Posterior Predictive Checks for Four spOccupancy Models ===")

# Load required packages
start_time <- start_timer()
if (!require(spOccupancy)) {
  log_message("ERROR: spOccupancy package not installed. Please install it first.")
  quit(status = 1)
}
log_message("â spOccupancy package loaded successfully")
end_timer(start_time, "Loading packages")

# Create output directory
if (!dir.exists("output")) {
  dir.create("output", recursive = TRUE)
  log_message("Created output directory")
}

# Define model files to process
# Update these paths according to your actual file locations
model_files <- list(
  sfMsPGOcc = "models/sfMsPGOcc/model_20250718_sfMsPGOcc_2010-2024_nthin500_nburn5e+05_nchain1_nsample6000_nfactors3_neighbors10.RData",
  spPGOcc = "models/spPGOcc/model_20250710_spPGOcc_2010-2024_nthin1000_nburn1e+06_nchain1_nsample6000_neighbors10.RData",
  msPGOcc = "models/msPGOcc/model_20250710_msPGOcc_2010-2024_nthin100_nburn1e+05_nchain4_nsample250000.RData",
  PGOcc = "models/PGOcc/model_20250710_PGOcc_2010-2024_5r_nthin100_nburn1e+05_nchain4_nsample250000.RData"
)

# List to store all results
all_results <- list()

# Function to run all PPCs for a given model
run_ppc_for_model <- function(model_obj, model_name, model_type) {
  log_message(paste("\n--- Running PPCs for", model_name, "---"))
  results <- list()
  
  # Freeman-Tukey PPC grouped by sites
  log_message("Running Freeman-Tukey PPC (group by sites)...")
  start_ppc_time <- start_timer()
  ft1_name <- paste0("ppc.", model_type, ".ft1")
  results[[ft1_name]] <- handle_error(
    ppcOcc(model_obj, fit.stat = 'freeman-tukey', group = 1),
    paste("Freeman-Tukey PPC by sites for", model_name)
  )
  end_timer(start_ppc_time, "Freeman-Tukey PPC by sites")
  
  if (!is.null(results[[ft1_name]])) {
    log_message(paste("Summary of Freeman-Tukey PPC by sites for", model_name, ":"))
    print(summary(results[[ft1_name]]))
  }
  
  # Freeman-Tukey PPC grouped by replicates
  log_message("Running Freeman-Tukey PPC (group by replicates)...")
  start_ppc_time <- start_timer()
  ft2_name <- paste0("ppc.", model_type, ".ft2")
  results[[ft2_name]] <- handle_error(
    ppcOcc(model_obj, fit.stat = 'freeman-tukey', group = 2),
    paste("Freeman-Tukey PPC by replicates for", model_name)
  )
  end_timer(start_ppc_time, "Freeman-Tukey PPC by replicates")
  
  if (!is.null(results[[ft2_name]])) {
    log_message(paste("Summary of Freeman-Tukey PPC by replicates for", model_name, ":"))
    print(summary(results[[ft2_name]]))
  }
  
  # Chi-Squared PPC grouped by sites
  log_message("Running Chi-Squared PPC (group by sites)...")
  start_ppc_time <- start_timer()
  cq1_name <- paste0("ppc.", model_type, ".cq1")
  results[[cq1_name]] <- handle_error(
    ppcOcc(model_obj, fit.stat = 'chi-squared', group = 1),
    paste("Chi-Squared PPC by sites for", model_name)
  )
  end_timer(start_ppc_time, "Chi-Squared PPC by sites")
  
  if (!is.null(results[[cq1_name]])) {
    log_message(paste("Summary of Chi-Squared PPC by sites for", model_name, ":"))
    print(summary(results[[cq1_name]]))
  }
  
  # Chi-Squared PPC grouped by replicates
  log_message("Running Chi-Squared PPC (group by replicates)...")
  start_ppc_time <- start_timer()
  cq2_name <- paste0("ppc.", model_type, ".cq2")
  results[[cq2_name]] <- handle_error(
    ppcOcc(model_obj, fit.stat = 'chi-squared', group = 2),
    paste("Chi-Squared PPC by replicates for", model_name)
  )
  end_timer(start_ppc_time, "Chi-Squared PPC by replicates")
  
  if (!is.null(results[[cq2_name]])) {
    log_message(paste("Summary of Chi-Squared PPC by replicates for", model_name, ":"))
    print(summary(results[[cq2_name]]))
  }
  
  # Calculate WAIC
  log_message(paste("Calculating WAIC for", model_name, "..."))
  start_waic_time <- start_timer()
  waic_name <- paste0("waic.", model_type)
  
  # Determine if model is multi-species or single-species
  if (model_type %in% c("sfMsPGOcc", "msPGOcc")) {
    # Multi-species models
    results[[waic_name]] <- handle_error(
      waicOcc(model_obj, by.sp = TRUE),
      paste("WAIC calculation for", model_name)
    )
  } else {
    # Single-species models
    results[[waic_name]] <- handle_error(
      waicOcc(model_obj),
      paste("WAIC calculation for", model_name)
    )
  }
  
  end_timer(start_waic_time, paste("WAIC calculation for", model_name))
  
  if (!is.null(results[[waic_name]])) {
    log_message(paste("WAIC for", model_name, "calculated successfully"))
    print(results[[waic_name]])
  }
  
  return(results)
}

# Process each model
for (model_type in names(model_files)) {
  model_file <- model_files[[model_type]]
  
  if (!is.na(model_file) && file.exists(model_file)) {
    start_time <- start_timer()
    log_message(paste("\n=== Processing", model_type, "model ==="))
    log_message(paste("Loading model from:", model_file))
    
    # Load model into a new environment
    model_env <- new.env()
    load_result <- handle_error(
      load(model_file, envir = model_env),
      paste("loading", model_type, "model")
    )
    
    if (!is.null(load_result)) {
      # Find the model object in the environment
      obj_name <- paste0("out.", model_type)
      
      if (exists(obj_name, envir = model_env)) {
        model_obj <- get(obj_name, envir = model_env)
        log_message(paste(model_type, "model loaded successfully"))
        
        # Run all PPCs for this model
        model_results <- run_ppc_for_model(model_obj, model_type, model_type)
        
        # Add results to the main list
        all_results <- c(all_results, model_results)
        
      } else {
        log_message(paste("ERROR: Could not find", obj_name, "in the loaded model file"))
        log_message(paste("Available objects:", paste(ls(envir = model_env), collapse = ", ")))
      }
    }
    
    end_timer(start_time, paste(model_type, "model processing"))
    
  } else {
    log_message(paste("WARNING:", model_type, "model file not found:", model_file))
    log_message(paste("Creating NULL placeholders for", model_type, "results"))
    
    # Create NULL placeholders for missing models
    for (suffix in c(".ft1", ".ft2", ".cq1", ".cq2")) {
      all_results[[paste0("ppc.", model_type, suffix)]] <- NULL
    }
    all_results[[paste0("waic.", model_type)]] <- NULL
  }
}

# Save all results
start_time <- start_timer()
output_file <- "output/ppc_results_four_models.RData"
log_message(paste("\n=== Saving results to", output_file, "==="))

save_result <- handle_error({
  # Extract all objects from the results list to the current environment
  list2env(all_results, environment())
  
  # Save all objects
  save(
    # Spatial Factor Multi-Species Model results
    ppc.sfMsPGOcc.ft1, ppc.sfMsPGOcc.ft2,
    ppc.sfMsPGOcc.cq1, ppc.sfMsPGOcc.cq2,
    waic.sfMsPGOcc,
    
    # Spatial Single-Species Model results
    ppc.spPGOcc.ft1, ppc.spPGOcc.ft2,
    ppc.spPGOcc.cq1, ppc.spPGOcc.cq2,
    waic.spPGOcc,
    
    # Non-spatial Multi-Species Model results
    ppc.msPGOcc.ft1, ppc.msPGOcc.ft2,
    ppc.msPGOcc.cq1, ppc.msPGOcc.cq2,
    waic.msPGOcc,
    
    # Non-spatial Single-Species Model results
    ppc.PGOcc.ft1, ppc.PGOcc.ft2,
    ppc.PGOcc.cq1, ppc.PGOcc.cq2,
    waic.PGOcc,
    
    file = output_file
  )
}, "saving results")

if (!is.null(save_result)) {
  log_message("â Results saved successfully")
} else {
  log_message("â Failed to save results")
}
end_timer(start_time, "Saving results")

# Final summary
log_message("\n=== FINAL SUMMARY ===")
log_message("Posterior predictive check script completed")
log_message(paste("See log file for details:", log_file))

# Print a summary of which checks were successful
successful_checks <- names(all_results)[!sapply(all_results, is.null)]
if (length(successful_checks) > 0) {
  log_message(paste("\nSuccessful checks:", length(successful_checks), "out of", length(all_results)))
  log_message("Completed tests:")
  for (check in successful_checks) {
    log_message(paste("  â", check))
  }
} else {
  log_message("No checks completed successfully")
}

# Print WAIC comparison if multiple models are available
waic_results <- all_results[grep("^waic\\.", names(all_results))]
waic_results <- waic_results[!sapply(waic_results, is.null)]

if (length(waic_results) > 1) {
  log_message("\n=== WAIC COMPARISON ===")
  for (model in names(waic_results)) {
    waic_val <- waic_results[[model]]
    if (is.list(waic_val) && "WAIC" %in% names(waic_val)) {
      log_message(paste(model, "- WAIC:", round(waic_val$WAIC, 2)))
    }
  }
}

log_message("\n=== SCRIPT EXECUTION COMPLETED ===")