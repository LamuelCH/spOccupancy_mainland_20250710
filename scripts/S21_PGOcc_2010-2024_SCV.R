# =============================================================================
# Single-Species Occupancy Model (PGOcc) Fitting Script
# Enhanced Version with Configurable K-Fold spatial blocking Cross-Validation
# =============================================================================

# ===========================================================================
# CONFIGURATION SECTION - MODIFY THESE PARAMETERS AS NEEDED
# ===========================================================================
K_FOLDS <- 2  # <-- CHANGE THIS VALUE TO SET NUMBER OF FOLDS (2, 3, 5, 10, etc.)
CV_TYPE <- "SCV"  # spatial blocking Cross-Validation
# ===========================================================================

# STEP 1: INITIAL SETUP AND ERROR HANDLING ===================================
tryCatch({
  
  cat("=== STEP 1: INITIAL SETUP ===\n")
  cat("*** ", K_FOLDS, "-FOLD spatial blocking CROSS-VALIDATION CONFIGURED ***\n\n", sep = "")
  
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
  cat("  spOccupancy package loaded successfully\n")
  
  if(!require(blockCV, quietly = TRUE)) {
    stop("Package 'blockCV' not available. Please install it first.")
  }
  cat("  blockCV package loaded successfully\n")
  
  if(!require(sf, quietly = TRUE)) {
    stop("Package 'sf' not available. Please install it first.")
  }
  cat("  sf package loaded successfully\n")
  
  if(!require(terra, quietly = TRUE)) {
    stop("Package 'terra' not available. Please install it first.")
  }
  cat("  terra package loaded successfully\n")
  
  # Set random seed for reproducibility
  set.seed(123)
  cat("  Random seed set to 123\n")
  
  # STEP 2: MCMC PARAMETER CONFIGURATION ====================================
  cat("\n=== STEP 2: MCMC PARAMETER CONFIGURATION ===\n")
  
  # Define parameter arrays for batch processing
  thinning_factors <- c(100)
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
  n.omp.threads <- 1
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
  cat("  Data loaded successfully\n")
  
  # Display data structure
  cat("\nData Structure Summary:\n")
  str(data)
  
  # Validate required components
  required_components <- c("y", "occ.covs", "det.covs", "coords")
  missing_components <- required_components[!required_components %in% names(data)]
  if(length(missing_components) > 0) {
    stop(paste("Missing required data components:", paste(missing_components, collapse=", ")))
  }
  cat("  All required data components present:", paste(required_components, collapse = ", "), "\n")
  
  # STEP 4: SPECIES SELECTION AND ANALYSIS ==================================
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
  sp.ordered <- c("Dasyurus")
  
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
  cat("â Species detection data reordered successfully\n")
  
  # STEP 5: spatial blocking CROSS-VALIDATION SETUP ==================
  cat("\n=== STEP 5: spatial blocking CROSS-VALIDATION SETUP ===\n")
  
  # Convert coordinates to sf object
  cat("Converting coordinates to sf object for spatial analysis...\n")
  coords_sf <- try(st_as_sf(as.data.frame(data.ordered$coords), 
                            coords = c("X", "Y"), 
                            crs = "EPSG:3577"))
  
  if(inherits(coords_sf, "try-error")) {
    stop("Error converting coordinates to sf object. Check coords data.")
  }
  cat("  Coordinates converted to sf object successfully\n")
  
  # Load environmental data for clustering
  cat("Loading environmental data for clustering...\n")
  bio <- rast("input/CHELSA_bio_1981-2010_V.2.1_EPSG3577.tif")
  bio <- bio[[c("bio5","bio6", "bio15")]]
  terrain <- rast("input/unmasked_env_terrain_EPSG3577.tif")
  terrain <- terrain[[c("dem", "tpi")]]
  fpc <- rast("input/env_foliage_EPSG3577.tif")
  env_stack <- c(bio, terrain, fpc)
  cat("  Environmental stack created with", nlyr(env_stack), "layers\n")
  
  # Create spatial blocks using spatial blocking
  cat("Creating spatial blocks using spatial blocking...\n")
  set.seed(123)
  sb <- cv_cluster(x = coords_sf,
                   # r = env_stack,
                   k = K_FOLDS, 
                   scale = TRUE)
  
  cat("  Spatial blocks created successfully\n")
  
  # Analyze fold assignments
  fold_ids <- sb$folds_ids
  fold_table <- table(fold_ids)
  cat("  Fold size distribution:\n")
  print(fold_table)
  
  # STEP 6: MODEL FORMULATION ===============================================
  cat("\n=== STEP 6: MODEL FORMULATION ===\n")
  
  # Define model formulas
  occ.formula <- ~ scale(bio5) + I(scale(bio5)^2) + 
    scale(bio6) + I(scale(bio6)^2) + 
    scale(bio15) +
    scale(dem) + I(scale(dem)^2) +
    scale(tpi) + I(scale(tpi)^2) +
    scale(fpc) + I(scale(fpc)^2)
  
  det.formula <- ~ scale(effort) + (1 | project)
  
  cat("Model Formulas:\n")
  cat("Occupancy formula: ", deparse(occ.formula), "\n")
  cat("Detection formula: ", deparse(det.formula), "\n")
  
  # STEP 7: INITIAL VALUES AND PRIORS ======================================
  cat("\n=== STEP 7: INITIAL VALUES AND PRIORS ===\n")
  
  # Set initial values
  z_inits <- apply(data.ordered$y, 1, max, na.rm = TRUE)
  inits <- list(beta = 0, alpha = 0, z = z_inits)
  cat("Initial values set.\n")
  
  # Set priors
  priors <- list(beta.normal = list(mean = 0, var = 2.72),
                 alpha.normal = list(mean = 0, var = 2.72))
  cat("Priors set.\n")
  
  # STEP 8: CROSS-VALIDATION MODEL FITTING =================================
  cat("\n=== STEP 8: spatial blocking CROSS-VALIDATION MODEL FITTING ===\n")
  
  output_dir <- "models/PGOcc"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  cat("Output directory:", output_dir, "\n")
  
  verbose <- TRUE
  n.report <- 250
  cv.threads <- min(K_FOLDS, 8)
  
  start_time <- Sys.time()
  cat("Model fitting started at:", format(start_time), "\n")
  
  model_result <- tryCatch({
    PGOcc(
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
      n.chains = n.chains,
      k.fold = K_FOLDS,
      k.fold.threads = cv.threads,
      k.fold.seed = 123,
      k.fold.only = TRUE,
      custom.cluster = sb$folds_ids
    )
  }, error = function(e) {
    cat("\n  ERROR in model fitting:", e$message, "\n")
    return(NULL)
  })
  
  end_time <- Sys.time()
  runtime <- end_time - start_time
  cat("\n--- MODEL FITTING COMPLETED ---\n")
  cat("Total runtime:", format(runtime), "\n")
  
  # STEP 9: MODEL SAVING AND SUMMARY ========================================
  cat("\n=== STEP 9: MODEL SAVING AND SUMMARY ===\n")
  
  if(!is.null(model_result)) {
    out.PGOcc <- model_result
    
    output_file <- file.path(output_dir, paste0(
      format(start_time, "%Y%m%d"), "_", CV_TYPE, "_", K_FOLDS, "fold_model_PGOcc_2010-2024",
      "_nthin", n.thin,
      "_nburn", n.burn, 
      "_nchain", n.chains,
      "_nsample", n.samples,
      ".RData"
    ))
    
    save(out.PGOcc, file = output_file)
    cat("  Model saved successfully to:", output_file, "\n")
    
    cat("\nCross-Validation Performance:\n")
    cat("  - Mean deviance:", round(mean(out.PGOcc$k.fold.deviance), 2), "\n")
    
  } else {
    cat("  MODEL FITTING FAILED\n")
  }
  
}, error = function(e) {
  cat("\n=== CRITICAL ERROR ===\n")
  cat("Error message:", e$message, "\n")
  quit(status = 1)
})

cat("\n=== SCRIPT EXECUTION COMPLETED ===\n")