# =============================================================================
# Multi-Species Occupancy Model (msPGOcc) Fitting Script
# Cleaned and Enhanced Version with Detailed Progress Reporting
# =============================================================================

# STEP 1: INITIAL SETUP AND ERROR HANDLING ===================================
tryCatch({
  
  cat("=== STEP 1: INITIAL SETUP ===\n")
  
  # Set up logging infrastructure
  log_file <- file.path("logs", paste0("job_", Sys.getenv("SLURM_ARRAY_TASK_ID"), "_", 
                                       format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
  dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
  cat("Log file created:", log_file, "\n")
  
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  if(length(args) == 0) {
    job_index <- 0  
    cat("No job index provided, using default value 0 for local testing\n")
  } else {
    job_index <- as.integer(args[1])
    if(is.na(job_index)) {
      stop("Invalid job index: must be an integer")
    }
    cat("Job index from SLURM:", job_index, "\n")
  }
  
  # Load required packages
  cat("Loading required packages...\n")
  if(!require(spOccupancy, quietly = TRUE)) {
    stop("Package 'spOccupancy' not available. Please install it first.")
  }
  cat("✓ spOccupancy package loaded successfully\n")
  
  # Set random seed for reproducibility
  set.seed(123)
  cat("✓ Random seed set to 123\n")
  
  # STEP 2: MCMC PARAMETER CONFIGURATION ====================================
  cat("\n=== STEP 2: MCMC PARAMETER CONFIGURATION ===\n")
  
  # Define parameter arrays
  thinning_factors <- c(1, 10, 25, 50, 100)
  cat("Available thinning factors:", paste(thinning_factors, collapse = ", "), "\n")
  
  # Validate job index
  if(job_index + 1 > length(thinning_factors)) {
    stop(sprintf("Job index %d out of range (max: %d)", job_index, length(thinning_factors) - 1))
  }
  
  # Set MCMC parameters
  n.thin <- thinning_factors[job_index + 1]
  n.total <- 6000
  n.chains <- 4
  n.burn <- 1000 * n.thin
  n.samples <- ((n.total * n.thin) / n.chains) + n.burn
  
  # Calculate effective samples
  effective_samples <- (n.samples - n.burn) / n.thin
  
  cat("MCMC Configuration:\n")
  cat("  - Thinning factor (n.thin):", n.thin, "\n")
  cat("  - Total iterations per chain (n.samples):", n.samples, "\n")
  cat("  - Burn-in iterations (n.burn):", n.burn, "\n")
  cat("  - Number of chains (n.chains):", n.chains, "\n")
  cat("  - Effective samples per chain:", effective_samples, "\n")
  cat("  - Total effective samples:", effective_samples * n.chains, "\n")
  
  # System resources
  n.omp.threads <- 8
  cat("  - OpenMP threads:", n.omp.threads, "\n")
  
  # STEP 3: DATA LOADING AND VALIDATION =====================================
  cat("\n=== STEP 3: DATA LOADING AND VALIDATION ===\n")
  
  # Find and set working directory
  data_file <- "input/data_top12_2010-2024_5r.RData"
  if(!file.exists(data_file)) {
    possible_dirs <- c(".", "..", "../..")
    found <- FALSE
    for(dir in possible_dirs) {
      if(file.exists(file.path(dir, data_file))) {
        setwd(dir)
        found <- TRUE
        cat("Working directory set to:", getwd(), "\n")
        break
      }
    }
    if(!found) {
      stop("Input data file not found. Check working directory.")
    }
  }
  
  # Load data with validation
  cat("Loading data file:", data_file, "\n")
  tryCatch({
    load(data_file)
    if(!exists("data")) {
      stop("Data object not found in loaded file")
    }
  }, error = function(e) {
    stop(paste("Error loading data file:", e$message))
  })
  cat("✓ Data loaded successfully\n")
  
  # Display data structure
  cat("\nData Structure Summary:\n")
  str(data)
  
  # Validate required components
  required_components <- c("y", "occ.covs", "det.covs")
  missing_components <- required_components[!required_components %in% names(data)]
  if(length(missing_components) > 0) {
    stop(paste("Missing required data components:", paste(missing_components, collapse=", ")))
  }
  cat("✓ All required data components present:", paste(required_components, collapse = ", "), "\n")
  
  # STEP 4: SPECIES ANALYSIS AND ORDERING ===================================
  cat("\n=== STEP 4: SPECIES ANALYSIS AND ORDERING ===\n")
  
  # Extract and analyze species
  sp.names <- rownames(data$y)
  n.species <- length(sp.names)
  n.sites <- ncol(data$y)
  n.visits <- dim(data$y)[3]
  
  cat("Species Information:\n")
  cat("  - Number of species:", n.species, "\n")
  cat("  - Number of sites:", n.sites, "\n")
  cat("  - Number of visits:", n.visits, "\n")
  cat("  - Original species order:", paste(sp.names, collapse = ", "), "\n")
  
  # Calculate and display occurrence probabilities
  occ_probs <- apply(data$y, 1, mean, na.rm = TRUE)
  cat("\nRaw occurrence probabilities:\n")
  for(i in 1:length(sp.names)) {
    cat(sprintf("  %s: %.3f\n", sp.names[i], occ_probs[i]))
  }
  
  # Define new species ordering (from most to least detectable/common)
  sp.ordered <- c(
    "Rattus",    # Factor 1
    "Dasyurus",    # Factor 2 
    "Vombatus",    # Factor 3
    "Tachyglossus",      
    "Canis",
    "Felis",
    "Vulpes",
    "Perameles",
    "Antechinus", 
    "Trichosurus",
    "Oryctolagus",
    "Wallabia"
  )
  
  
  cat("\nNew species ordering (ecological rationale):\n")
  for(i in 1:length(sp.ordered)) {
    if(sp.ordered[i] %in% sp.names) {
      cat(sprintf("  %d. %s (prob: %.3f)\n", i, sp.ordered[i], occ_probs[sp.ordered[i]]))
    }
  }
  
  # Validate species ordering
  missing_sp <- sp.ordered[!sp.ordered %in% sp.names]
  if(length(missing_sp) > 0) {
    stop(paste("Species in ordering not found in data:", paste(missing_sp, collapse=", ")))
  }
  
  # Reorder detection data
  y.new <- data$y[sp.ordered, ,]
  data.ordered <- data
  data.ordered$y <- y.new
  cat("✓ Species detection data reordered successfully\n")
  
  # STEP 5: MODEL FORMULATION ===============================================
  cat("\n=== STEP 5: MODEL FORMULATION ===\n")
  
  # Define model formulas
  occ.formula <- ~ scale(bio5) + I(scale(bio5)^2) + 
    scale(bio6) + I(scale(bio6)^2) + 
    scale(bio15) +
    scale(dem) + I(scale(dem)^2) +
    scale(tpi) + I(scale(tpi)^2) +
    scale(fpc) + I(scale(fpc)^2)
  
  det.formula <- ~ scale(effort) + (1|project)
  
  cat("Occupancy formula:\n")
  cat("  ", deparse(occ.formula), "\n")
  cat("Detection formula:\n")
  cat("  ", deparse(det.formula), "\n")
  
  # Count covariates
  occ_terms <- attr(terms(occ.formula), "term.labels")
  det_terms <- attr(terms(det.formula), "term.labels")
  cat("  - Occupancy covariates:", length(occ_terms), "\n")
  cat("  - Detection covariates:", length(det_terms), "\n")
  
  # STEP 6: INITIAL VALUES AND PRIORS ======================================
  cat("\n=== STEP 6: INITIAL VALUES AND PRIORS ===\n")
  
  # Set initial values
  z_inits <- apply(data.ordered$y, c(1, 2), max, na.rm = TRUE)
  
  inits <- list(
    alpha.comm = 0,      # Community-level detection intercept
    beta.comm = 0,       # Community-level occupancy intercept
    beta = 0,            # Species-specific occupancy effects
    alpha = 0,           # Species-specific detection effects
    tau.sq.beta = 1,     # Occupancy variance
    tau.sq.alpha = 1,    # Detection variance
    z = z_inits          # Occupancy states
  )
  
  cat("Initial values configured:\n")
  cat("  - Community intercepts set to 0\n")
  cat("  - Variance parameters set to 1\n")
  cat("  - Occupancy states initialized from detection data\n")
  cat("  - Initial z matrix dimensions:", dim(z_inits), "\n")
  
  # Set priors
  priors <- list(
    beta.comm.normal = list(mean = 0, var = 2.72),    # Weakly informative
    alpha.comm.normal = list(mean = 0, var = 2.72),   # Weakly informative
    tau.sq.beta.ig = list(a = 0.1, b = 0.1),         # Inverse gamma
    tau.sq.alpha.ig = list(a = 0.1, b = 0.1)         # Inverse gamma
  )
  
  cat("Prior distributions:\n")
  cat("  - Community effects: Normal(0, 2.72)\n")
  cat("  - Variance parameters: InverseGamma(0.1, 0.1)\n")
  
  # STEP 7: MODEL FITTING ===================================================
  cat("\n=== STEP 7: MODEL FITTING ===\n")
  
  # Prepare output directory
  output_dir <- "models/msPGOcc"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  cat("Output directory:", output_dir, "\n")
  
  # Configure reporting
  verbose <- TRUE
  n.report <- 250
  cat("Progress will be reported every", n.report, "iterations\n")
  
  # Record start time
  start_time <- Sys.time()
  cat("Model fitting started at:", format(start_time), "\n")
  cat("Estimated runtime: ~", round((n.samples * n.chains) / 400), "minutes\n")
  
  # Fit the model
  cat("\n--- MODEL FITTING IN PROGRESS ---\n")
  model_result <- tryCatch({
    msPGOcc(
      occ.formula = occ.formula,
      det.formula = det.formula,
      data = data.ordered,
      inits = inits,
      priors = priors,
      n.omp.threads = n.omp.threads,
      verbose = verbose,
      n.samples = n.samples,
      n.report = n.report,
      n.burn = n.burn,
      n.thin = n.thin,
      n.chains = n.chains
    )
  }, error = function(e) {
    cat("ERROR in model fitting:", e$message, "\n")
    return(NULL)
  }, warning = function(w) {
    cat("WARNING in model fitting:", w$message, "\n")
    invokeRestart("muffleWarning")
  })
  
  # Calculate runtime
  end_time <- Sys.time()
  runtime <- end_time - start_time
  cat("\n--- MODEL FITTING COMPLETED ---\n")
  cat("Start time:", format(start_time), "\n")
  cat("End time:", format(end_time), "\n")
  cat("Total runtime:", format(runtime), "\n")
  
  # STEP 8: MODEL SAVING AND SUMMARY ========================================
  cat("\n=== STEP 8: MODEL SAVING AND SUMMARY ===\n")
  
  if(!is.null(model_result)) {
    # Assign result
    out.msPGOcc <- model_result
    
    # Generate output filename
    output_file <- file.path(output_dir, paste0(
      "model_", format(start_time, "%Y%m%d"), "_msPGOcc_2010-2024",
      "_nthin", n.thin,
      "_nburn", n.burn, 
      "_nchain", n.chains,
      "_nsample", n.samples,
      ".RData"
    ))
    
    # Save model
    tryCatch({
      save(out.msPGOcc, file = output_file)
      cat("✓ Model saved successfully to:", output_file, "\n")
      
      # Display file size
      file_size <- file.size(output_file)
      cat("  File size:", round(file_size / (1024^2), 2), "MB\n")
      
    }, error = function(e) {
      cat("ERROR saving model:", e$message, "\n")
      # Try backup location
      backup_file <- file.path(tempdir(), basename(output_file))
      cat("Attempting backup save to:", backup_file, "\n")
      save(out.msPGOcc, file = backup_file)
      cat("✓ Model saved to backup location\n")
    })
    
    # Display basic model summary
    cat("\nModel Summary:\n")
    cat("  - Species fitted:", nrow(data.ordered$y), "\n")
    cat("  - Sites analyzed:", ncol(data.ordered$y), "\n")
    cat("  - Effective samples per chain:", (n.samples - n.burn) / n.thin, "\n")
    cat("  - Total effective samples:", ((n.samples - n.burn) / n.thin) * n.chains, "\n")
    
    cat("\n✓ MODEL FITTING SUCCESSFUL!\n")
    
  } else {
    cat("✗ MODEL FITTING FAILED - Check error messages above\n")
    quit(status = 1)
  }
  
}, error = function(e) {
  # Global error handler
  cat("\n=== CRITICAL ERROR ===\n")
  cat("Error message:", e$message, "\n")
  cat("Call stack:\n")
  print(sys.calls())
  
  # Log error to file
  error_file <- file.path("logs", paste0("error_", Sys.getenv("SLURM_ARRAY_TASK_ID"), "_", 
                                         format(start_time, "%Y%m%d_%H%M%S"), ".txt"))
  dir.create(dirname(error_file), showWarnings = FALSE, recursive = TRUE)
  
  writeLines(c(
    paste("Timestamp:", Sys.time()),
    paste("Error:", e$message),
    "Call stack:",
    paste(capture.output(print(sys.calls())), collapse = "\n")
  ), error_file)
  
  cat("Error logged to:", error_file, "\n")
  quit(status = 1)
})

cat("\n=== SCRIPT EXECUTION COMPLETED ===\n")
cat("Session info:\n")
print(sessionInfo())