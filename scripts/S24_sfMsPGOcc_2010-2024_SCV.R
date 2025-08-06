# =============================================================================
# Spatial Factor Multi-Species Occupancy Model (sfMsPGOcc) Fitting Script
# Enhanced Version with Configurable Spatial Block Cross-Validation
# =============================================================================

# ===========================================================================
# CONFIGURATION SECTION - MODIFY THESE PARAMETERS AS NEEDED
# ===========================================================================
K_FOLDS <- 10  # <-- CHANGE THIS VALUE TO SET NUMBER OF FOLDS (2, 3, 5, 10, etc.)
CV_TYPE <- "SCV"  # Spatial Clustering Cross-Validation
# ===========================================================================

# STEP 1: INITIAL SETUP AND ERROR HANDLING ===================================
tryCatch({
  
  cat("=== STEP 1: INITIAL SETUP ===\n")
  cat("*** ", K_FOLDS, "-FOLD Spatial Clustering CROSS-VALIDATION CONFIGURED ***\n\n", sep = "")
  
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
  cat(" spOccupancy package loaded successfully\n")
  
  if(!require(blockCV, quietly = TRUE)) {
    stop("Package 'blockCV' not available. Please install it first.")
  }
  cat(" blockCV package loaded successfully\n")
  
  if(!require(sf, quietly = TRUE)) {
    stop("Package 'sf' not available. Please install it first.")
  }
  cat(" sf package loaded successfully\n")
  
  if(!require(terra, quietly = TRUE)) {
    stop("Package 'terra' not available. Please install it first.")
  }
  cat(" terra package loaded successfully\n")
  
  # Set random seed for reproducibility
  set.seed(123)
  cat(" Random seed set to 123\n")
  
  # STEP 2: MCMC PARAMETER CONFIGURATION ====================================
  cat("\n=== STEP 2: MCMC PARAMETER CONFIGURATION ===\n")
  
  # Define parameter arrays for batch processing
  thinning_factors <- c(500)
  cat("Available thinning factors:", paste(thinning_factors, collapse = ", "), "\n")
  
  # Validate job index
  if(job_index + 1 > length(thinning_factors)) {
    stop(sprintf("Job index %d out of range (max: %d)", job_index, length(thinning_factors) - 1))
  }
  
  # Set MCMC parameters for spatial factor model
  n.thin <- thinning_factors[job_index + 1]
  n.sample <- 6000                        # Samples per chain after burn-in
  n.burn <- 1000 * n.thin                # Burn-in period
  batch.length <- 25                      # Batch length for adaptive sampling
  n.batch <- ceiling((n.burn/n.thin + n.sample) * n.thin / batch.length)
  n.chains <- 1                          # Multiple chains for convergence
  n.neighbors <- 10                      # NNGP neighbors for spatial approximation
  n.factors <- 3                         # Number of latent factors
  
  cat("MCMC Configuration for Spatial Factor Model:\n")
  cat("  - Thinning factor (n.thin):", n.thin, "\n")
  cat("  - Samples per chain (n.sample):", n.sample, "\n")
  cat("  - Burn-in iterations (n.burn):", n.burn, "\n")
  cat("  - Batch length:", batch.length, "\n")
  cat("  - Number of batches (n.batch):", n.batch, "\n")
  cat("  - Number of chains (n.chains):", n.chains, "\n")
  cat("  - NNGP neighbors:", n.neighbors, "\n")
  cat("  - Latent factors:", n.factors, "\n")
  cat("  - Total iterations per chain:", n.batch * batch.length, "\n")
  cat("  - Effective samples per chain:", n.sample, "\n")
  cat("  - Total effective samples:", n.sample * n.chains, "\n")
  
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
  cat(" Data loaded successfully\n")
  
  # Scale coordinates to kilometers for spatial modeling
  cat("Scaling coordinates to kilometers for spatial analysis...\n")
  data$coords <- data$coords / 1000
  coord_range <- apply(data$coords, 2, range)
  cat("Coordinate ranges (km):\n")
  cat("  X: [", coord_range[1,1], ", ", coord_range[2,1], "]\n")
  cat("  Y: [", coord_range[1,2], ", ", coord_range[2,2], "]\n")
  
  # Display data structure
  cat("\nMulti-Species Spatial Data Structure Summary:\n")
  str(data)
  
  # Validate required spatial components
  required_components <- c("y", "occ.covs", "det.covs", "coords")
  missing_components <- required_components[!required_components %in% names(data)]
  if(length(missing_components) > 0) {
    stop(paste("Missing required data components:", paste(missing_components, collapse=", ")))
  }
  cat(" All required spatial data components present:", paste(required_components, collapse = ", "), "\n")
  
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
  
  # Define ecological ordering for factor model identifiability
  sp.ordered <- c(
    "Rattus",        # Factor 1 anchor
    "Dasyurus",      # Factor 2 anchor
    "Vombatus",      # Factor 3 anchor
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
  
  cat("\nEcological species ordering for factor model:\n")
  for(i in 1:length(sp.ordered)) {
    if(sp.ordered[i] %in% sp.names) {
      cat(sprintf("  %2d. %s (prob: %.3f)", i, sp.ordered[i], occ_probs[sp.ordered[i]]))
      if(i <= n.factors) cat(" [Factor ", i, " anchor]")
      cat("\n")
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
  cat(" Species detection data reordered successfully\n")
  
  # STEP 5: SPATIAL BLOCK CROSS-VALIDATION SETUP ============================
  cat("\n=== STEP 5: SPATIAL BLOCK CROSS-VALIDATION SETUP ===\n")
  
  # Convert coordinates to sf object
  cat("\nConverting coordinates to sf object for spatial analysis...\n")
  coords_sf <- try(st_as_sf(as.data.frame(data.ordered$coords), 
                            coords = c("X", "Y"), 
                            crs = "EPSG:3577"))
  
  if(inherits(coords_sf, "try-error")) {
    stop("Error converting coordinates to sf object. Check coords data.")
  }
  cat(" Coordinates converted to sf object successfully\n")
  cat("  - CRS: EPSG:3577 (GDA94 / Australian Albers)\n")
  cat("  - Number of points:", nrow(coords_sf), "\n")
  
  # Load environmental data for clustering
  cat("\nLoading environmental data for spatial clustering...\n")
  
  # Load bioclimatic variables
  cat("  Loading CHELSA bioclimatic variables...\n")
  bio <- rast("input/CHELSA_bio_1981-2010_V.2.1_EPSG3577.tif")
  bio <- bio[[c("bio5","bio6", "bio15")]]
  cat("     Selected: bio5 (Max temp warmest month), bio6 (Min temp coldest month), bio15 (Precipitation seasonality)\n")
  
  # Load terrain variables
  cat("  Loading terrain variables...\n")
  terrain <- rast("input/unmasked_env_terrain_EPSG3577.tif")
  terrain <- terrain[[c("dem", "tpi")]]
  cat("     Selected: dem (Digital Elevation Model), tpi (Topographic Position Index)\n")
  
  # Load vegetation variable
  cat("  Loading vegetation variable...\n")
  fpc <- rast("input/env_foliage_EPSG3577.tif")
  cat("     Selected: fpc (Foliage Projective Cover)\n")
  
  # Stack environmental layers
  env_stack <- c(bio, terrain, fpc)
  cat("\n Environmental stack created with", nlyr(env_stack), "layers\n")
  
  # Convert environmental stack coordinates to kilometers
  cat("\nConverting environmental raster coordinates to kilometers...\n")
  
  # Get current extent and resolution
  ext_m <- ext(env_stack)
  res_m <- res(env_stack)
  
  cat("  Original extent (m):\n")
  cat("    X: [", ext_m[1], ", ", ext_m[2], "]\n")
  cat("    Y: [", ext_m[3], ", ", ext_m[4], "]\n")
  cat("  Original resolution (m):", res_m[1], "x", res_m[2], "\n")
  
  # Create new extent in kilometers
  ext_km <- ext(ext_m[1]/1000, ext_m[2]/1000, ext_m[3]/1000, ext_m[4]/1000)
  
  # Create new raster with km coordinates
  env_stack_km <- rast(nrows = nrow(env_stack), 
                       ncols = ncol(env_stack),
                       nlyrs = nlyr(env_stack),
                       extent = ext_km,
                       crs = crs(env_stack))
  
  # Copy values and layer names
  values(env_stack_km) <- values(env_stack)
  names(env_stack_km) <- names(env_stack)
  
  cat("\n  Converted extent (km):\n")
  cat("    X: [", ext_km[1], ", ", ext_km[2], "]\n")
  cat("    Y: [", ext_km[3], ", ", ext_km[4], "]\n")
  cat("  Converted resolution (km):", res(env_stack_km)[1], "x", res(env_stack_km)[2], "\n")
  cat(" Environmental raster coordinates converted to kilometers\n")
  
  # Create spatial blocks using Spatial Clustering
  cat("\nCreating spatial blocks using Spatial Clustering...\n")
  cat("  Method: cv_cluster (Spatial Clustering)\n")
  cat("  Number of folds:", K_FOLDS, "\n")
  cat("  Scaling: TRUE (standardize environmental variables)\n")
  
  # Perform clustering
  set.seed(123)  # For reproducibility
  sb <- cv_cluster(x = coords_sf,
                   k = K_FOLDS, 
                   scale = TRUE)
  
  cat("\n Spatial blocks created successfully\n")
  
  # Analyze spatial block results
  cat("\nSpatial Block Analysis:\n")
  
  # Get fold assignments
  fold_ids <- sb$folds_ids
  fold_table <- table(fold_ids)
  
  cat("  Fold size distribution:\n")
  for(i in 1:K_FOLDS) {
    if(i %in% names(fold_table)) {
      cat(sprintf("    Fold %2d: %4d sites (%.1f%%)\n", 
                  i, fold_table[as.character(i)], fold_table[as.character(i)]/sum(fold_table)*100))
    }
  }
  
  cat("\n  Fold statistics:\n")
  cat("    - Min fold size:", min(fold_table), "\n")
  cat("    - Max fold size:", max(fold_table), "\n")
  cat("    - Mean fold size:", round(mean(fold_table), 1), "\n")
  cat("    - SD fold size:", round(sd(fold_table), 1), "\n")
  
  # Calculate spatial autocorrelation of folds
  cat("\n  Spatial characteristics:\n")
  cat("    - Clustering method: K-means on environmental space\n")
  cat("    - Environmental variables used:", paste(names(env_stack_km), collapse = ", "), "\n")
  
  # Save CV plot
  cat("\nGenerating spatial block visualization...\n")
  diagnostics_dir <- paste0("output/sfMsPGOcc/diagnostics/", CV_TYPE, "_", K_FOLDS, "fold")
  dir.create(diagnostics_dir, showWarnings = FALSE, recursive = TRUE)
  plot_file <- file.path(diagnostics_dir, paste0("spatial_blocks_", K_FOLDS, "fold_", format(Sys.Date(), "%Y%m%d"), ".png"))
  
  png(plot_file, width = 1200, height = 800, res = 150)
  cv_plot(cv = sb, 
          x = coords_sf,
          r = env_stack_km)
  dev.off()
  
  cat(" Spatial block plot saved to:", plot_file, "\n")
  
  # STEP 6: SPATIAL CORRELATION ANALYSIS ====================================
  cat("\n=== STEP 6: SPATIAL CORRELATION ANALYSIS ===\n")
  
  # Calculate spatial distances
  cat("Calculating pairwise distances between sites...\n")
  dist.data <- try(dist(data.ordered$coords))
  if(inherits(dist.data, "try-error")) {
    stop("Error calculating distances between sites. Check coords data.")
  }
  
  # Analyze spatial scale
  min.dist <- min(dist.data, na.rm = TRUE)
  max.dist <- max(dist.data, na.rm = TRUE)
  mean.dist <- mean(dist.data, na.rm = TRUE)
  median.dist <- median(as.matrix(dist.data)[upper.tri(as.matrix(dist.data))], na.rm = TRUE)
  
  cat("Spatial Distance Summary (km):\n")
  cat("  - Minimum distance:", round(min.dist, 2), "\n")
  cat("  - Maximum distance:", round(max.dist, 2), "\n")
  cat("  - Mean distance:", round(mean.dist, 2), "\n")
  cat("  - Median distance:", round(median.dist, 2), "\n")
  
  # Analyze within-fold distances
  cat("\nWithin-fold distance analysis:\n")
  for(fold in 1:K_FOLDS) {
    fold_sites <- which(fold_ids == fold)
    if(length(fold_sites) > 1) {
      fold_coords <- data.ordered$coords[fold_sites, ]
      fold_dist <- dist(fold_coords)
      cat(sprintf("  Fold %2d - Mean distance: %.1f km, Max: %.1f km\n", 
                  fold, mean(fold_dist), max(fold_dist)))
    }
  }
  
  # Set spatial correlation model
  cov.model <- "exponential"
  cat("\n  - Spatial correlation model:", cov.model, "\n")
  
  # STEP 7: FACTOR MODEL SETUP ==============================================
  cat("\n=== STEP 7: FACTOR MODEL SETUP ===\n")
  
  # Set number of species
  N <- nrow(data.ordered$y)
  
  cat("Factor Model Configuration:\n")
  cat("  - Number of species (N):", N, "\n")
  cat("  - Number of latent factors:", n.factors, "\n")
  cat("  - Factor model type: Spatial with NNGP approximation\n")
  cat("  - Identifiability constraints: Upper triangular zero constraints\n")
  
  # Initialize factor loadings matrix (matching RCV script)
  cat("\nInitializing factor loadings matrix...\n")
  lambda.inits <- matrix(0, nrow = N, ncol = n.factors)
  
  # Set diagonal elements to 1
  diag(lambda.inits) <- 1
  # Set lower triangular elements to random values from a standard normal dist
  lambda.inits[lower.tri(lambda.inits)] <- rnorm(sum(lower.tri(lambda.inits)))
  
  rownames(lambda.inits) <- sp.ordered
  colnames(lambda.inits) <- paste0("Factor", 1:n.factors)
  
  cat("Factor loadings initialization matrix:\n")
  print(round(lambda.inits, 3))
  
  # Verify constraints
  cat("\nConstraint verification for factor model:\n")
  cat("  - Factor 1 anchor (", sp.ordered[1], "):", lambda.inits[1,1], "\n")
  cat("  - Factor 2 anchor (", sp.ordered[2], "):", lambda.inits[2,2], "\n")
  cat("  - Factor 3 anchor (", sp.ordered[3], "):", lambda.inits[3,3], "\n")
  cat("  - Upper triangle sum:", sum(lambda.inits[upper.tri(lambda.inits)]), "(should be 0)\n")
  
  # STEP 8: MODEL FORMULATION ===============================================
  cat("\n=== STEP 8: MODEL FORMULATION ===\n")
  
  # Define model formulas
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
  cat("\nModel complexity:\n")
  cat("  - Occupancy covariates:", length(occ_terms), "\n")
  cat("  - Detection covariates:", length(det_terms), "\n")
  cat("  - Latent factors:", n.factors, "\n")
  cat("  - Total species modeled:", N, "\n")
  
  # STEP 9: INITIAL VALUES AND PRIORS ======================================
  cat("\n=== STEP 9: INITIAL VALUES AND PRIORS ===\n")
  
  # Set initial values (matching RCV script)
  cat("Setting initial values...\n")
  z_inits <- apply(data.ordered$y, c(1, 2), max, na.rm = TRUE)
  
  # Set spatial decay parameter based on effective range
  phi_init <- round(3 / mean.dist, 3)
  
  inits <- list(
    alpha.comm = 0,
    beta.comm = 0,
    beta = 0,
    alpha = 0,
    tau.sq.beta = 1,
    tau.sq.alpha = 1,
    phi = phi_init,
    lambda = lambda.inits,
    z = z_inits
  )
  
  cat("Initial values configured:\n")
  cat("  - Community-level coefficients: 0\n")
  cat("  - Species-specific coefficients: 0\n")
  cat("  - Variance parameters: 1\n")
  cat("  - Spatial decay (phi):", phi_init, "\n")
  cat("  - Effective range: ~", round(3/phi_init, 1), "km\n")
  cat("  - Factor loadings: Random with constraints\n")
  cat("  - Initial occupancy states: From detection data\n")
  
  # Set prior distributions (matching RCV script)
  cat("\nSetting prior distributions...\n")
  
  # Validate distance bounds
  if(min.dist <= 0 || is.infinite(min.dist) || is.na(min.dist)) {
    warning("Invalid minimum distance detected. Using default value of 3km.")
    min.dist <- 3
  }
  
  # Set prior bounds based on data
  max.dist.prior <- max(dist.data, na.rm = TRUE)
  
  priors <- list(
    beta.comm.normal = list(mean = 0, var = 2.72),
    alpha.comm.normal = list(mean = 0, var = 2.72),
    tau.sq.beta.ig = list(a = 0.1, b = 0.1),
    tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
    phi.unif = list(3 / max.dist.prior, 3 / min.dist)
  )
  
  cat("Prior distributions:\n")
  cat("  - Community effects (beta, alpha): Normal(0, 2.72)\n")
  cat("  - Variance parameters (tau.sq): InverseGamma(0.1, 0.1)\n")
  cat("  - Spatial decay (phi): Uniform(", round(3/max.dist.prior, 4), ", ", 
      round(3/min.dist, 4), ")\n")
  cat("  - Effective range bounds: [", round(3/(3/min.dist), 1), ", ", 
      round(3/(3/max.dist.prior), 1), "] km\n")
  
  # Set tuning parameters
  tuning <- list(phi = 1)
  cat("\nMCMC tuning:\n")
  cat("  - Phi tuning parameter:", tuning$phi, "\n")
  cat("  - Target acceptance rate: 0.43\n")
  
  # STEP 10: CROSS-VALIDATION MODEL FITTING =================================
  cat("\n=== STEP 10: SPATIAL CROSS-VALIDATION MODEL FITTING ===\n")
  
  # Prepare output directory
  output_dir <- "models/sfMsPGOcc"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  cat("Output directory:", output_dir, "\n")
  
  # Configure model fitting
  verbose <- TRUE
  n.report <- 250
  accept.rate <- 0.43
  
  # Set CV threads based on number of folds
  cv.threads <- min(K_FOLDS, 10)  # Use up to 8 threads or number of folds, whichever is smaller
  
  cat("\nCross-validation configuration:\n")
  cat("  - K-fold cross-validation:", K_FOLDS, "folds\n")
  cat("  - Fold assignment: Spatial Clustering (spatial blocks)\n")
  cat("  - CV threads:", cv.threads, "\n")
  cat("  - CV seed: 123\n")
  cat("  - Mode: k.fold.only = TRUE (CV predictions only)\n")
  
  cat("\nModel fitting configuration:\n")
  cat("  - Progress reporting every:", n.report, "iterations\n")
  cat("  - Target acceptance rate:", accept.rate, "\n")
  cat("  - NNGP approximation: TRUE\n")
  cat("  - NNGP neighbors:", n.neighbors, "\n")
  cat("  - Spatial correlation model:", cov.model, "\n")
  
  # Record start time
  start_time <- Sys.time()
  cat("\nModel fitting started at:", format(start_time), "\n")
  total_iterations <- n.batch * batch.length * n.chains
  estimated_minutes <- total_iterations / 150
  cat("Estimated runtime: ~", round(estimated_minutes), "minutes\n")
  cat("Total iterations across all chains:", format(total_iterations, big.mark = ","), "\n")
  
  # Fit the spatial factor model with cross-validation
  cat("\n--- SPATIAL FACTOR", K_FOLDS, "FOLD Spatial Clustering CV MODEL FITTING IN PROGRESS ---\n")
  cat("Each of the", K_FOLDS, "folds will be held out and predicted using the remaining data...\n\n")
  
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
      n.chains = n.chains,
      k.fold = K_FOLDS,              # Using the configurable parameter
      k.fold.threads = cv.threads,   # Adjusted based on fold number
      k.fold.seed = 123,
      k.fold.only = TRUE,
      custom.cluster = sb$folds_ids  # Spatial Clustering fold assignments
    )
  }, error = function(e) {
    cat("\nÂ  ERROR in spatial factor model fitting:", e$message, "\n")
    return(NULL)
  }, warning = function(w) {
    cat("\nÂ   WARNING in spatial factor model fitting:", w$message, "\n")
    invokeRestart("muffleWarning")
  })
  
  # Calculate runtime
  end_time <- Sys.time()
  runtime <- end_time - start_time
  cat("\n--- SPATIAL FACTOR", K_FOLDS, "FOLD Spatial Clustering CV MODEL FITTING COMPLETED ---\n")
  cat("Start time:", format(start_time), "\n")
  cat("End time:", format(end_time), "\n")
  cat("Total runtime:", format(runtime), "\n")
  cat("Runtime per 1000 iterations:", format(runtime / (total_iterations/1000)), "\n")
  
  # STEP 11: MODEL SAVING AND RESULTS SUMMARY ==============================
  cat("\n=== STEP 11: MODEL SAVING AND RESULTS SUMMARY ===\n")
  
  if(!is.null(model_result)) {
    # Assign result
    out.sfMsPGOcc <- model_result
    
    # Generate output filename with execution date at the beginning
    output_file <- file.path(output_dir, paste0(
      format(start_time, "%Y%m%d"), "_", CV_TYPE, "_", K_FOLDS, "fold_model_sfMsPGOcc_2010-2024",
      "_nthin", n.thin,
      "_nburn", n.burn,
      "_nchain", n.chains,
      "_nsample", n.sample,
      "_nfactors", n.factors,
      "_neighbors", n.neighbors,
      ".RData"
    ))
    
    # Save model
    tryCatch({
      save(out.sfMsPGOcc, file = output_file)
      cat(" Cross-validation model saved successfully to:\n")
      cat("  ", output_file, "\n")
      
      # Display file information
      file_size <- file.size(output_file)
      cat("  File size:", round(file_size / (1024^2), 2), "MB\n")
      
    }, error = function(e) {
      cat("Â  ERROR saving model:", e$message, "\n")
      backup_file <- file.path(tempdir(), basename(output_file))
      cat("Attempting backup save to:", backup_file, "\n")
      save(out.sfMsPGOcc, file = backup_file)
      cat(" Model saved to backup location\n")
    })
    
    # Display model summary
    cat("\nSpatial Factor Cross-Validation Model Summary:\n")
    cat("  - Model type: sfMsPGOcc (Spatial Factor Multi-Species Occupancy)\n")
    cat("  - Cross-validation:", K_FOLDS, "fold Spatial Clustering\n")
    cat("  - Species fitted:", N, "species\n")
    cat("  - Sites analyzed:", n.sites, "\n")
    cat("  - Latent factors:", n.factors, "\n")
    cat("  - Spatial correlation model:", cov.model, "\n")
    cat("  - NNGP neighbors:", n.neighbors, "\n")
    cat("  - Effective samples per chain:", n.sample, "\n")
    cat("  - Total effective samples:", n.sample * n.chains, "\n")
    cat("  - Spatial extent (km):", round(max.dist, 1), "\n")
    cat("  - Computational time:", format(runtime), "\n")
    
    # Cross-validation results
    if("k.fold.deviance" %in% names(out.sfMsPGOcc)) {
      cat("\n", K_FOLDS, "Fold Cross-validation Results:\n", sep = "")
      cat("  - Mean deviance:", round(mean(out.sfMsPGOcc$k.fold.deviance), 2), "\n")
      cat("  - SD deviance:", round(sd(out.sfMsPGOcc$k.fold.deviance), 2), "\n")
      cat("  - Fold deviances:\n")
      for(i in 1:length(out.sfMsPGOcc$k.fold.deviance)) {
        cat(sprintf("    Fold %2d: %.2f\n", i, out.sfMsPGOcc$k.fold.deviance[i]))
      }
    }
    
    cat("\n SPATIAL FACTOR", K_FOLDS, "FOLD Spatial Clustering CV MODEL FITTING SUCCESSFUL!\n")
    
  } else {
    cat("Â  SPATIAL FACTOR", K_FOLDS, "FOLD Spatial Clustering CV MODEL FITTING FAILED\n")
    cat("Check error messages above for details\n")
    quit(status = 1)
  }
  
}, error = function(e) {
  # Global error handler
  cat("\n=== CRITICAL ERROR ===\n")
  cat("Error message:", e$message, "\n")
  cat("Call stack:\n")
  print(sys.calls())
  
  # Log error
  error_file <- file.path("logs", paste0("error_", Sys.getenv("SLURM_ARRAY_TASK_ID"), "_", 
                                         format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
  dir.create(dirname(error_file), showWarnings = FALSE, recursive = TRUE)
  
  writeLines(c(
    paste("Timestamp:", Sys.time()),
    paste("Model Type: sfMsPGOcc (Spatial Factor Multi-Species with", K_FOLDS, "fold", CV_TYPE, ")"),
    paste("Error:", e$message),
    "Call stack:",
    paste(capture.output(print(sys.calls())), collapse = "\n")
  ), error_file)
  
  cat("Error logged to:", error_file, "\n")
  quit(status = 1)
})

# STEP 12: PERFORMANCE METRICS AND DIAGNOSTICS ================================
cat("\n=== STEP 12: PERFORMANCE METRICS AND DIAGNOSTICS ===\n")

if(exists("out.sfMsPGOcc") && !is.null(out.sfMsPGOcc)) {
  
  # Save diagnostic plots
  cat("Generating diagnostic plots...\n")
  diagnostics_dir <- paste0("output/sfMsPGOcc/diagnostics/", CV_TYPE, "_", K_FOLDS, "fold")
  dir.create(diagnostics_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Save fold assignment summary
  fold_summary <- data.frame(
    fold = 1:K_FOLDS,
    n_sites = as.vector(table(sb$folds_ids)),
    prop_sites = as.vector(table(sb$folds_ids)) / length(sb$folds_ids)
  )
  
  write.csv(fold_summary, 
            file.path(diagnostics_dir, paste0(K_FOLDS, "fold_summary_", format(Sys.Date(), "%Y%m%d"), ".csv")),
            row.names = FALSE)
  
  cat(" Fold summary saved\n")
  
  # Save species ordering and factor model structure
  species_summary <- data.frame(
    species = sp.ordered,
    order = 1:length(sp.ordered),
    raw_occurrence = occ_probs[sp.ordered],
    factor_anchor = c(rep("yes", n.factors), rep("no", length(sp.ordered) - n.factors))
  )
  
  write.csv(species_summary,
            file.path(diagnostics_dir, paste0("species_summary_", format(Sys.Date(), "%Y%m%d"), ".csv")),
            row.names = FALSE)
  
  cat(" Species summary saved\n")
  
  # Create diagnostic summary file
  diag_file <- file.path(diagnostics_dir, paste0(
    K_FOLDS, "fold_", CV_TYPE, "_diagnostics_", 
    format(start_time, "%Y%m%d_%H%M%S"), ".txt"
  ))
  
  sink(diag_file)
  cat("=== ", K_FOLDS, " FOLD Spatial Clustering CV DIAGNOSTIC SUMMARY ===\n", sep = "")
  cat("Generated:", format(Sys.time()), "\n\n")
  
  cat("Model Configuration:\n")
  cat("  - K-fold CV:", K_FOLDS, "folds\n")
  cat("  - CV Type: Spatial Clustering\n")
  cat("  - Total sites:", n.sites, "\n")
  cat("  - Sites per fold: ~", round(n.sites/K_FOLDS), "\n")
  cat("  - Species:", N, "\n")
  cat("  - Factors:", n.factors, "\n")
  cat("  - MCMC samples:", n.sample, "\n")
  cat("  - Runtime:", format(runtime), "\n\n")
  
  if("k.fold.deviance" %in% names(out.sfMsPGOcc)) {
    cat("Cross-Validation Performance:\n")
    cat("  - Mean deviance:", round(mean(out.sfMsPGOcc$k.fold.deviance), 2), "\n")
    cat("  - SD deviance:", round(sd(out.sfMsPGOcc$k.fold.deviance), 2), "\n")
    cat("  - Min deviance:", round(min(out.sfMsPGOcc$k.fold.deviance), 2), "\n")
    cat("  - Max deviance:", round(max(out.sfMsPGOcc$k.fold.deviance), 2), "\n\n")
    
    cat("Fold-specific deviances:\n")
    for(i in 1:length(out.sfMsPGOcc$k.fold.deviance)) {
      cat(sprintf("  Fold %2d: %.2f\n", i, out.sfMsPGOcc$k.fold.deviance[i]))
    }
  }
  
  sink()
  cat(" Diagnostics saved to:", diag_file, "\n")
  
  # System information
  cat("\nSystem Information:\n")
  cat("  - R version:", R.version$version.string, "\n")
  cat("  - Platform:", R.version$platform, "\n")
  
  # Memory usage
  gc_info <- gc()
  cat("\nMemory Usage:\n")
  cat("  - Used (Mb):", round(sum(gc_info[, 2]), 1), "\n")
  cat("  - Max used (Mb):", round(sum(gc_info[, 6]), 1), "\n")
}

cat("\n=== ", K_FOLDS, " FOLD Spatial Clustering CV SCRIPT EXECUTION COMPLETED ===\n", sep = "")
cat("Session info:\n")
print(sessionInfo())