# =============================================================================
# Spatial Single-Species Occupancy Model (spPGOcc) Fitting Script
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
  
  # Parse command line arguments for SLURM array jobs
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
  
  # Define parameter arrays for batch processing
  thinning_factors <- c(250, 500, 1000)
  cat("Available thinning factors:", paste(thinning_factors, collapse = ", "), "\n")
  
  # Validate job index
  if(job_index + 1 > length(thinning_factors)) {
    stop(sprintf("Job index %d out of range (max: %d)", job_index, length(thinning_factors) - 1))
  }
  
  # Set MCMC parameters for spatial model
  n.thin <- thinning_factors[job_index + 1]
  n.sample <- 6000                        # Samples per chain after burn-in
  n.burn <- 1000 * n.thin                # Burn-in period
  batch.length <- 25                      # Batch length for adaptive sampling
  n.batch <- ceiling((n.burn/n.thin + n.sample) * n.thin / batch.length)
  n.chains <- 1                          # Multiple chains for convergence
  n.neighbors <- 10                      # NNGP neighbors for spatial approximation
  
  cat("MCMC Configuration for Spatial Model:\n")
  cat("  - Thinning factor (n.thin):", n.thin, "\n")
  cat("  - Samples per chain (n.sample):", n.sample, "\n")
  cat("  - Burn-in iterations (n.burn):", n.burn, "\n")
  cat("  - Batch length:", batch.length, "\n")
  cat("  - Number of batches (n.batch):", n.batch, "\n")
  cat("  - Number of chains (n.chains):", n.chains, "\n")
  cat("  - NNGP neighbors:", n.neighbors, "\n")
  cat("  - Total iterations per chain:", n.batch * batch.length, "\n")
  cat("  - Effective samples per chain:", n.sample, "\n")
  cat("  - Total effective samples:", n.sample * n.chains, "\n")
  
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
  
  # Scale coordinates to kilometers for spatial modeling
  cat("Scaling coordinates to kilometers for spatial analysis...\n")
  data$coords <- data$coords / 1000
  coord_range <- apply(data$coords, 2, range)
  cat("Coordinate ranges (km):\n")
  cat("  X: [", coord_range[1,1], ", ", coord_range[2,1], "]\n")
  cat("  Y: [", coord_range[1,2], ", ", coord_range[2,2], "]\n")
  
  # Display data structure
  cat("\nSpatial Data Structure Summary:\n")
  str(data)
  
  # Validate required spatial components
  required_components <- c("y", "occ.covs", "det.covs", "coords")
  missing_components <- required_components[!required_components %in% names(data)]
  if(length(missing_components) > 0) {
    stop(paste("Missing required data components:", paste(missing_components, collapse=", ")))
  }
  cat("✓ All required spatial data components present:", paste(required_components, collapse = ", "), "\n")
  
  # STEP 4: SPECIES ANALYSIS AND ORDERING ===================================
  cat("\n=== STEP 4: SPECIES ANALYSIS AND ORDERING ===\n")
  
  # Extract and analyze species information
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
  
  # Define new species ordering
  sp.ordered <- c("Dasyurus")   # For single-species model\n    
  
  
  cat("\nSpecies ordering:\n")
  for(i in 1:length(sp.ordered)) {
    if(sp.ordered[i] %in% sp.names) {
      cat(sprintf("  %d. %s (prob: %.3f)\n", i, sp.ordered[i], occ_probs[sp.ordered[i]]))
    }
  }
  
  # Validate and reorder species data
  missing_sp <- sp.ordered[!sp.ordered %in% sp.names]
  if(length(missing_sp) > 0) {
    stop(paste("Species in ordering not found in data:", paste(missing_sp, collapse=", ")))
  }
  
  y.new <- data$y[sp.ordered, ,]
  data.ordered <- data
  data.ordered$y <- y.new
  cat("✓ Species detection data reordered successfully\n")
  
  # STEP 5: SPATIAL CORRELATION SETUP ======================================
  cat("\n=== STEP 5: SPATIAL CORRELATION SETUP ===\n")
  
  # Calculate spatial distances for correlation modeling
  cat("Calculating pairwise distances between sites...\n")
  dist.data <- try(dist(data.ordered$coords))
  if(inherits(dist.data, "try-error")) {
    stop("Error calculating distances between sites. Check coords data.")
  }
  
  # Analyze spatial scale
  min.dist <- min(dist.data, na.rm = TRUE)
  max.dist <- max(dist.data, na.rm = TRUE)
  mean.dist <- mean(dist.data, na.rm = TRUE)
  
  cat("Spatial Distance Summary (km):\n")
  cat("  - Minimum distance:", round(min.dist, 2), "\n")
  cat("  - Maximum distance:", round(max.dist, 2), "\n")
  cat("  - Mean distance:", round(mean.dist, 2), "\n")
  
  # Set spatial correlation model
  cov.model <- "exponential"
  cat("  - Correlation model:", cov.model, "\n")
  
  # STEP 6: MODEL FORMULATION ===============================================
  cat("\n=== STEP 6: MODEL FORMULATION ===\n")
  
  # Define model formulas with comprehensive covariate set
  occ.formula <- ~ scale(bio5) + I(scale(bio5)^2) + 
    scale(bio6) + I(scale(bio6)^2) + 
    scale(bio15) +
    scale(dem) + I(scale(dem)^2) +
    scale(tpi) + I(scale(tpi)^2) +
    scale(fpc) + I(scale(fpc)^2)
  
  det.formula <- ~ scale(effort) + (1 | project)
  
  cat("Model Formulas:\n")
  cat("Occupancy formula:\n")
  cat("  ", deparse(occ.formula), "\n")
  cat("Detection formula:\n")
  cat("  ", deparse(det.formula), "\n")
  
  # Count model terms
  occ_terms <- attr(terms(occ.formula), "term.labels")
  det_terms <- attr(terms(det.formula), "term.labels")
  cat("  - Occupancy covariates:", length(occ_terms), "\n")
  cat("  - Detection covariates:", length(det_terms), "\n")
  
  # STEP 7: INITIAL VALUES AND PRIORS ======================================
  cat("\n=== STEP 7: INITIAL VALUES AND PRIORS ===\n")
  
  # Set initial values for spatial model
  cat("Setting initial values...\n")
  z_inits <- apply(data.ordered$y, 1, max, na.rm = TRUE)
  
  # Set reasonable spatial decay parameter initial value
  # Based on effective range of ~30km for mammal communities
  phi_init <- round(3 / mean.dist, 3) 
  
  # Model design matrices
  X <- model.matrix(occ.formula, data.ordered$occ.covs)
  
  inits <- list(
    beta = 0,                 # Occupancy effects
    alpha = 0,                # Detection effects
    sigma.sq = 2,             # Spatial variance
    w = rep(0, nrow(X)),      # Spatial random effects
    phi = phi_init,           # Spatial decay parameter
    z = z_inits               # Occupancy states
  )
  
  cat("Initial values configured:\n")
  cat("  - Occupancy coefficients: 0\n")
  cat("  - Detection coefficients: 0\n")
  cat("  - Spatial variance (sigma.sq): 2\n")
  cat("  - Spatial decay (phi):", phi_init, "(effective range)\n")
  cat("  - Spatial random effects: 0\n")
  cat("  - Occupancy states: From detection data\n")
  
  # Set prior distributions
  cat("\nSetting prior distributions...\n")
  
  # Validate distance bounds for spatial priors
  if(min.dist <= 0 || is.infinite(min.dist) || is.na(min.dist)) {
    warning("Invalid minimum distance detected. Using default value of 50km.")
    min.dist <- 50
  }
  
  # min.dist = 3 # 3 km
  # max.dist.prior = 200 # 40, 10, 200 km respectively
  min.dist = min(dist.data, na.rm = TRUE)
  max.dist.prior = max(dist.data, na.rm = TRUE)
  
  priors <- list(
    beta.normal = list(mean = 0, var = 2.72),           # Weakly informative
    alpha.normal = list(mean = 0, var = 2.72),          # Weakly informative
    sigma.sq.ig = c(2, 2),                              # Spatial variance
    phi.unif = c(3 / max.dist.prior, 3 / min.dist)  # Spatial decay bounds
  )
  
  cat("Prior distributions:\n")
  cat("  - Occupancy effects: Normal(0, 2.72)\n")
  cat("  - Detection effects: Normal(0, 2.72)\n")
  cat("  - Spatial variance: InverseGamma(2, 2)\n")
  cat("  - Spatial decay (phi): Uniform(", round(3/max.dist.prior, 4), ", ", 
      round(3/min.dist, 4), ")\n")
  cat("  - Effective range bounds: [", round(3/(3/min.dist), 1), ", ", 
      round(3/(3/max.dist.prior), 1), "] km\n")
  
  # Set tuning parameters for adaptive MCMC
  tuning <- list(phi = 1)
  cat("  - Tuning parameter (phi):", tuning$phi, "\n")
  
  # STEP 8: MODEL FITTING ===================================================
  cat("\n=== STEP 8: SPATIAL MODEL FITTING ===\n")
  
  # Prepare output directory
  output_dir <- "models/spPGOcc"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  cat("Output directory:", output_dir, "\n")
  
  # Configure progress reporting
  verbose <- TRUE
  n.report <- 250
  accept.rate <- 0.43
  cat("Model fitting configuration:\n")
  cat("  - Progress reporting every:", n.report, "iterations\n")
  cat("  - Target acceptance rate:", accept.rate, "\n")
  cat("  - NNGP approximation: TRUE\n")
  cat("  - Spatial correlation model:", cov.model, "\n")
  
  # Record start time and estimate duration
  start_time <- Sys.time()
  cat("Model fitting started at:", format(start_time), "\n")
  total_iterations <- n.batch * batch.length * n.chains
  estimated_minutes <- total_iterations / 250  # Rough estimate
  cat("Estimated runtime: ~", round(estimated_minutes), "minutes\n")
  cat("Total iterations across all chains:", total_iterations, "\n")
  
  # Fit the spatial model
  cat("\n--- SPATIAL SINGLE-SPECIES MODEL FITTING IN PROGRESS ---\n")
  model_result <- tryCatch({
    spPGOcc(
      occ.formula = occ.formula,
      det.formula = det.formula,
      data = data.ordered,
      inits = inits,
      n.batch = n.batch,
      batch.length = batch.length,
      accept.rate = accept.rate,
      priors = priors,
      cov.model = cov.model,
      tuning = tuning,
      n.omp.threads = n.omp.threads,
      verbose = verbose,
      NNGP = TRUE,
      n.neighbors = n.neighbors,
      n.report = n.report,
      n.burn = n.burn,
      n.thin = n.thin,
      n.chains = n.chains
    )
  }, error = function(e) {
    cat("ERROR in spatial model fitting:", e$message, "\n")
    return(NULL)
  }, warning = function(w) {
    cat("WARNING in spatial model fitting:", w$message, "\n")
    invokeRestart("muffleWarning")
  })
  
  # Calculate runtime
  end_time <- Sys.time()
  runtime <- end_time - start_time
  cat("\n--- SPATIAL SINGLE-SPECIES MODEL FITTING COMPLETED ---\n")
  cat("Start time:", format(start_time), "\n")
  cat("End time:", format(end_time), "\n")
  cat("Total runtime:", format(runtime), "\n")
  cat("Runtime per 1000 iterations:", format(runtime / (total_iterations/1000)), "\n")
  
  # STEP 9: MODEL SAVING AND SUMMARY =======================================
  cat("\n=== STEP 9: MODEL SAVING AND SUMMARY ===\n")
  
  if(!is.null(model_result)) {
    # Assign result
    out.spPGOcc <- model_result
    
    # Generate comprehensive output filename
    output_file <- file.path(output_dir, paste0(
      "model_", format(start_time, "%Y%m%d"), "_spPGOcc_2010-2024",
      "_nthin", n.thin,
      "_nburn", n.burn,
      "_nchain", n.chains,
      "_nsample", n.sample,
      "_neighbors", n.neighbors,
      ".RData"
    ))
    
    # Save model with error handling
    tryCatch({
      save(out.spPGOcc, file = output_file)
      cat("✓ Spatial model saved successfully to:", output_file, "\n")
      
      # Display file information
      file_size <- file.size(output_file)
      cat("  File size:", round(file_size / (1024^2), 2), "MB\n")
      
    }, error = function(e) {
      cat("ERROR saving model:", e$message, "\n")
      # Try backup location
      backup_file <- file.path(tempdir(), basename(output_file))
      cat("Attempting backup save to:", backup_file, "\n")
      save(out.spPGOcc, file = backup_file)
      cat("✓ Model saved to backup location\n")
    })
    
    # Display comprehensive model summary
    cat("\nSpatial Model Summary:\n")
    cat("  - Model type: spPGOcc (Spatial Single-Species Occupancy)\n")
    cat("  - Species fitted:", n.species, "\n")
    cat("  - Sites analyzed:", nrow(data.ordered$y), "\n")
    cat("  - Spatial correlation model:", cov.model, "\n")
    cat("  - NNGP neighbors:", n.neighbors, "\n")
    cat("  - Effective samples per chain:", (n.batch * batch.length - n.burn) / n.thin, "\n")
    cat("  - Total effective samples:", ((n.batch * batch.length - n.burn) / n.thin) * n.chains, "\n")
    cat("  - Spatial extent (km):", round(max.dist - min.dist, 1), "\n")
    cat("  - Computational time:", format(runtime), "\n")
    
    # Basic convergence information
    if("rhat" %in% names(out.spPGOcc)) {
      cat("\nConvergence Information Available:\n")
      cat("  - Rhat values computed: ✓\n")
      cat("  - ESS values computed: ✓\n")
      cat("  - Use summary(out.spPGOcc) for detailed diagnostics\n")
    }
    
    cat("\n✓ SPATIAL SINGLE-SPECIES MODEL FITTING SUCCESSFUL!\n")
    cat("Next steps:\n")
    cat("  1. Check convergence: summary(out.spPGOcc)\n")
    cat("  2. Plot diagnostics: plot(out.spPGOcc, 'beta', density=FALSE)\n")
    cat("  3. Examine spatial parameters: plot(out.spPGOcc, 'theta', density=FALSE)\n")
    
  } else {
    cat("✗ SPATIAL SINGLE-SPECIES MODEL FITTING FAILED - Check error messages above\n")
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
    paste("Model Type: spPGOcc (Spatial Single-Species)"),
    paste("Error:", e$message),
    "Call stack:",
    paste(capture.output(print(sys.calls())), collapse = "\n")
  ), error_file)
  
  cat("Error logged to:", error_file, "\n")
  quit(status = 1)
})

cat("\n=== SPATIAL SINGLE-SPECIES MODEL SCRIPT EXECUTION COMPLETED ===\n")
cat("Session info:\n")
print(sessionInfo())