# =============================================================================
# Spatial Factor Multi-Species Occupancy Model (sfMsPGOcc) Fitting Script
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
  thinning_factors <- c(1, 10, 25, 50, 100)
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
  
  # Define new species ordering for factor model identifiability
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
  
  
  cat("\nNew species ordering (factor model anchoring):\n")
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
  
  # STEP 6: FACTOR MODEL SETUP ==============================================
  cat("\n=== STEP 6: FACTOR MODEL SETUP ===\n")
  
  # Set number of latent factors
  n.factors <- 3
  N <- nrow(data.ordered$y)
  
  cat("Factor Model Configuration:\n")
  cat("  - Number of species (N):", N, "\n")
  cat("  - Number of latent factors:", n.factors, "\n")
  cat("  - Factor model type: Spatial with NNGP approximation\n")
  
  # Initialize factor loadings matrix with careful constraints
  cat("\nInitializing factor loadings matrix...\n")
  lambda.inits <- matrix(0, nrow = N, ncol = n.factors)
  
  # Set diagonal elements to 1
  diag(lambda.inits) <- 1
  # Set lower triangular elements to random values from a standard normal dist
  lambda.inits[lower.tri(lambda.inits)] <- rnorm(sum(lower.tri(lambda.inits)))
  
  # # Set constraint values for model identifiability
  # # Based on Carvalho et al. (2008) approach - reordering and transforming loadings
  # # NEW ORDER: Vombatus (1st), Macropus (2nd), Dasyurus (3rd), then others
  # 
  # # FACTOR 1 LOADINGS (Vombatus anchor - transformed from converged results)
  # # Vombatus had highest magnitude (2.2393), now scaled to 1.0 (scale factor: 0.446568)
  # factor1_values <- c(
  #   1.0000,    # Vombatus - anchor species (fixed at 1.0) [was 2.2393]
  #   0.0979,    # Macropus - large herbivore (slightly similar) [was 0.2193, scaled]
  #   -0.6013,   # Dasyurus - native predator (different from Vombatus) [was -1.3465, scaled]
  #   0.4466,    # Rattus - commensal rodent (similar to Vombatus) [was 1.0000, scaled]
  #   0.2476,    # Felis - introduced predator (somewhat similar) [was 0.5545, scaled]
  #   0.2052,    # Vulpes - introduced predator (somewhat similar) [was 0.4594, scaled]
  #   0.0484,    # Perameles - native small mammal (slightly similar) [was 0.1083, scaled]
  #   0.9134,    # Antechinus - native small mammal (very similar pattern) [was 2.0454, scaled]
  #   0.1926,    # Trichosurus - arboreal marsupial (somewhat similar) [was 0.4313, scaled]
  #   0.3589,    # Oryctolagus - introduced herbivore (similar) [was 0.8036, scaled]
  #   0.1307,    # Wallabia - medium herbivore (slightly similar) [was 0.2927, scaled]
  #   -0.0748    # Canis - large competitor (slightly different) [was -0.1676, scaled]
  # )
  # 
  # # FACTOR 2 LOADINGS (Macropus anchor - transformed from converged results)
  # # Macropus had high magnitude (-2.6642), scaled to 1.0 with sign flip (scale factor: 0.375347)
  # factor2_values <- c(
  #   0.0000,    # Vombatus - fixed constraint (upper triangular)
  #   1.0000,    # Macropus - anchor species (fixed at 1.0) [was -2.6642, now positive]
  #   0.0000,    # Dasyurus - fixed constraint (for Factor 3)
  #   0.0000,    # Rattus - fixed constraint (was original anchor)
  #   0.3511,    # Felis - different from Macropus pattern [was -0.9354, transformed]
  #   0.9476,    # Vulpes - very different from Macropus [was -2.5247, transformed]
  #   -0.0793,   # Perameles - slightly similar to Macropus [was 0.2112, transformed]
  #   0.0517,    # Antechinus - slightly different from Macropus [was -0.1378, transformed]
  #   0.4970,    # Trichosurus - quite different from Macropus [was -1.3240, transformed]
  #   0.9044,    # Oryctolagus - very different from Macropus [was -2.4096, transformed]
  #   0.4080,    # Wallabia - fellow macropod but different size [was -1.0869, transformed]
  #   0.1223     # Canis - somewhat different from Macropus [was -0.3259, transformed]
  # )
  # 
  # # FACTOR 3 LOADINGS (Dasyurus anchor - from original Factor 2 relationships)
  # # Uses original Factor 2 loadings when Dasyurus was anchor (1.0) - NO scaling applied
  # # Reordered to match new species order with proper constraints
  # factor3_values <- c(
  #   0.0000,    # Vombatus - fixed constraint (upper triangular)
  #   0.0000,    # Macropus - fixed constraint (upper triangular)
  #   1.0000,    # Dasyurus - anchor species (fixed at 1.0)
  #   0.0000,    # Rattus - constraint (was 0.0000 in original Factor 2)
  #   -0.9354,   # Felis - original relationship to Dasyurus from Factor 2
  #   -2.5247,   # Vulpes - original relationship to Dasyurus from Factor 2
  #   0.2112,    # Perameles - original relationship to Dasyurus from Factor 2
  #   -0.1378,   # Antechinus - original relationship to Dasyurus from Factor 2
  #   -1.3240,   # Trichosurus - original relationship to Dasyurus from Factor 2
  #   -2.4096,   # Oryctolagus - original relationship to Dasyurus from Factor 2
  #   -1.0869,   # Wallabia - original relationship to Dasyurus from Factor 2
  #   -0.3259    # Canis - original relationship to Dasyurus from Factor 2
  # )
  
  # if(length(factor1_values) != N) {
  #   stop(sprintf("Loading values length (%d) doesn't match species count (%d)", 
  #                length(factor1_values), N))
  # }
  # 
  # if(length(factor2_values) != N) {
  #   stop(sprintf("Loading values length (%d) doesn't match species count (%d)", 
  #                length(factor2_values), N))
  # }
  # 
  # if(length(factor3_values) != N) {
  #   stop(sprintf("Loading values length (%d) doesn't match species count (%d)", 
  #                length(factor3_values), N))
  # }
  
  # # Assign to matrix
  # lambda.inits[,1] <- factor1_values
  # lambda.inits[,2] <- factor2_values  
  # lambda.inits[,3] <- factor3_values
  
  rownames(lambda.inits) <- sp.ordered
  colnames(lambda.inits) <- paste0("Factor", 1:n.factors)
  
  cat("Factor loadings initialization matrix:\n")
  print(round(lambda.inits, 3))
  
  # Verify constraints for identifiability
  cat("\nConstraint verification for 3-factor model:\n")
  cat("  - Factor 1 anchor (Vombatus, pos 1):", lambda.inits[1,1], "\n")
  cat("  - Factor 2 anchor (Macropus, pos 2):", lambda.inits[2,2], "\n")
  cat("  - Factor 3 anchor (Dasyurus, pos 3):", lambda.inits[3,3], "\n")
  cat("  - Upper triangular constraints:\n")
  cat("    Rattus Factor 2:", lambda.inits[1,2], "\n")
  cat("    Rattus Factor 3:", lambda.inits[1,3], "\n") 
  cat("    Dasyurus Factor 3:", lambda.inits[2,3], "\n")
  cat("  - All upper triangle elements sum:", sum(lambda.inits[upper.tri(lambda.inits)]), "(should be 0)\n")
  
  # STEP 7: MODEL FORMULATION ===============================================
  cat("\n=== STEP 7: MODEL FORMULATION ===\n")
  
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
  
  # STEP 8: INITIAL VALUES AND PRIORS ======================================
  cat("\n=== STEP 8: INITIAL VALUES AND PRIORS ===\n")
  
  # Set initial values for spatial factor model
  cat("Setting initial values...\n")
  z_inits <- apply(data.ordered$y, c(1, 2), max, na.rm = TRUE)
  
  # Set reasonable spatial decay parameter initial value
  # Based on effective range of ~30km for mammal communities
  # phi_init <- c(3/12, 3/5, 3/4)  # median distance
  phi_init = round(3/mean.dist, 3)
    
  inits <- list(
    alpha.comm = 0,           # Community-level detection intercept
    beta.comm = 0,            # Community-level occupancy intercept
    beta = 0,                 # Species-specific occupancy effects
    alpha = 0,                # Species-specific detection effects
    tau.sq.beta = 1,          # Occupancy variance
    tau.sq.alpha = 1,         # Detection variance
    phi = phi_init,           # Spatial decay parameter
    lambda = lambda.inits,    # Factor loadings
    z = z_inits               # Occupancy states
  )
  
  cat("Initial values configured:\n")
  cat("  - Community intercepts: 0\n")
  cat("  - Variance parameters: 1\n")
  cat("  - Spatial decay (phi):", phi_init, "(effective range)\n")
  cat("  - Factor loadings: From previous model estimates\n")
  cat("  - Occupancy states: From detection data\n")
  
  # Set prior distributions
  cat("\nSetting prior distributions...\n")
  
  # Validate distance bounds for spatial priors
  if(min.dist <= 0 || is.infinite(min.dist) || is.na(min.dist)) {
    warning("Invalid minimum distance detected. Using default value of 50km.")
    min.dist <- 100
  }
  
  # min.dist = c(3, 3, 3) # 3 km
  # max.dist.prior = c(40, 10, 200) # 40, 10, 200 km respectively
  min.dist = min(dist.data, na.rm = TRUE)
  max.dist.prior = max(dist.data, na.rm = TRUE)
  
  priors <- list(
    beta.comm.normal = list(mean = 0, var = 2.72),           # Weakly informative
    alpha.comm.normal = list(mean = 0, var = 2.72),          # Weakly informative
    tau.sq.beta.ig = list(a = 0.1, b = 0.1),                # Inverse gamma
    tau.sq.alpha.ig = list(a = 0.1, b = 0.1),               # Inverse gamma
    rho.unif = list(a = -1, b = 1),                         # Correlation bounds
    sigma.sq.t.ig = list(a = 0.1, b = 0.1),                 # Temporal variance
    phi.unif = list(3 / max.dist.prior, 3 / min.dist)       # Spatial decay bounds
  )
  
  cat("Prior distributions:\n")
  cat("  - Community effects: Normal(0, 2.72)\n")
  cat("  - Variance parameters: InverseGamma(0.1, 0.1)\n")
  cat("  - Spatial correlation: Uniform(-1, 1)\n")
  cat("  - Spatial decay (phi): Uniform(", round(3/max.dist.prior, 4), ", ", 
      round(3/min.dist, 4), ")\n")
  cat("  - Effective range bounds: [", round(3/(3/min.dist), 1), ", ", 
      round(3/(3/max.dist.prior), 1), "] km\n")
  
  # Set tuning parameters for adaptive MCMC
  tuning <- list(phi = 1)
  cat("  - Tuning parameter (phi):", tuning$phi, "\n")
  
  # STEP 9: MODEL FITTING ===================================================
  cat("\n=== STEP 9: SPATIAL FACTOR MODEL FITTING ===\n")
  
  # Prepare output directory
  output_dir <- "models/sfMsPGOcc"
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
  estimated_minutes <- total_iterations / 150  # Rough estimate
  cat("Estimated runtime: ~", round(estimated_minutes), "minutes\n")
  cat("Total iterations across all chains:", total_iterations, "\n")
  
  # Fit the spatial factor model
  cat("\n--- SPATIAL FACTOR MODEL FITTING IN PROGRESS ---\n")
  model_result <- tryCatch({
    sfMsPGOcc(
      occ.formula = occ.formula,
      det.formula = det.formula,
      data = data.ordered,
      inits = inits,
      n.batch = n.batch,
      batch.length = batch.length,
      accept.rate = accept.rate,
      priors = priors,
      n.factors = n.factors,
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
  cat("\n--- SPATIAL FACTOR MODEL FITTING COMPLETED ---\n")
  cat("Start time:", format(start_time), "\n")
  cat("End time:", format(end_time), "\n")
  cat("Total runtime:", format(runtime), "\n")
  cat("Runtime per 1000 iterations:", format(runtime / (total_iterations/1000)), "\n")
  
  # STEP 10: MODEL SAVING AND SUMMARY =======================================
  cat("\n=== STEP 10: MODEL SAVING AND SUMMARY ===\n")
  
  if(!is.null(model_result)) {
    # Assign result
    out.sfMsPGOcc <- model_result
    
    # Generate comprehensive output filename
    output_file <- file.path(output_dir, paste0(
      "model_", format(start_time, "%Y%m%d"), "_sfMsPGOcc_2010-2024",
      "_nthin", n.thin,
      "_nburn", n.burn,
      "_nchain", n.chains,
      "_nsample", n.sample,
      "_nfactors", n.factors,
      "_neighbors", n.neighbors,
      ".RData"
    ))
    
    # Save model with error handling
    tryCatch({
      save(out.sfMsPGOcc, file = output_file)
      cat("✓ Spatial factor model saved successfully to:", output_file, "\n")
      
      # Display file information
      file_size <- file.size(output_file)
      cat("  File size:", round(file_size / (1024^2), 2), "MB\n")
      
    }, error = function(e) {
      cat("ERROR saving model:", e$message, "\n")
      # Try backup location
      backup_file <- file.path(tempdir(), basename(output_file))
      cat("Attempting backup save to:", backup_file, "\n")
      save(out.sfMsPGOcc, file = backup_file)
      cat("✓ Model saved to backup location\n")
    })
    
    # Display comprehensive model summary
    cat("\nSpatial Factor Model Summary:\n")
    cat("  - Model type: sfMsPGOcc (Spatial Factor Multi-Species)\n")
    cat("  - Species fitted:", nrow(data.ordered$y), "\n")
    cat("  - Sites analyzed:", ncol(data.ordered$y), "\n")
    cat("  - Latent factors:", n.factors, "\n")
    cat("  - Spatial correlation model:", cov.model, "\n")
    cat("  - NNGP neighbors:", n.neighbors, "\n")
    cat("  - Effective samples per chain:", (n.batch * batch.length - n.burn) / n.thin, "\n")
    cat("  - Total effective samples:", ((n.batch * batch.length - n.burn) / n.thin) * n.chains, "\n")
    cat("  - Spatial extent (km):", round(max.dist - min.dist, 1), "\n")
    cat("  - Computational time:", format(runtime), "\n")
    
    # Basic convergence information
    if("rhat" %in% names(out.sfMsPGOcc)) {
      cat("\nConvergence Information Available:\n")
      cat("  - Rhat values computed: ✓\n")
      cat("  - ESS values computed: ✓\n")
      cat("  - Use summary(out.sfMsPGOcc) for detailed diagnostics\n")
    }
    
    cat("\n✓ SPATIAL FACTOR MODEL FITTING SUCCESSFUL!\n")
    cat("Next steps:\n")
    cat("  1. Check convergence: summary(out.sfMsPGOcc)\n")
    cat("  2. Plot diagnostics: plot(out.sfMsPGOcc, 'lambda', density=FALSE)\n")
    cat("  3. Examine spatial parameters: plot(out.sfMsPGOcc, 'theta', density=FALSE)\n")
    
  } else {
    cat("✗ SPATIAL FACTOR MODEL FITTING FAILED - Check error messages above\n")
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
    paste("Model Type: sfMsPGOcc (Spatial Factor)"),
    paste("Error:", e$message),
    "Call stack:",
    paste(capture.output(print(sys.calls())), collapse = "\n")
  ), error_file)
  
  cat("Error logged to:", error_file, "\n")
  quit(status = 1)
})

cat("\n=== SPATIAL FACTOR MODEL SCRIPT EXECUTION COMPLETED ===\n")
cat("Session info:\n")
print(sessionInfo())