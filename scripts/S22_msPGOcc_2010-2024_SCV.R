# =============================================================================
# Non-Spatial Multi-Species Occupancy Model (msPGOcc) Fitting Script
# Enhanced Version with Configurable Spatial Blocking K-Fold Cross-Validation
# =============================================================================

# ===========================================================================
# CONFIGURATION SECTION - MODIFY THESE PARAMETERS AS NEEDED
# ===========================================================================
K_FOLDS <- 2  # <-- CHANGE THIS VALUE TO SET NUMBER OF FOLDS (2, 3, 5, 10, etc.)
CV_TYPE <- "SCV"  # Spatial Blocking Cross-Validation
# ===========================================================================

# STEP 1: INITIAL SETUP AND ERROR HANDLING ===================================
tryCatch({
  
  cat("=== STEP 1: INITIAL SETUP ===\n")
  cat("*** ", K_FOLDS, "-FOLD Spatial Blocking CROSS-VALIDATION CONFIGURED ***\n\n", sep = "")
  
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
  cat(" Data loaded successfully\n")
  
  # Scale coordinates to kilometers for spatial clustering
  cat("Scaling coordinates to kilometers for Spatial Blocking...\n")
  data$coords <- data$coords / 1000
  coord_range <- apply(data$coords, 2, range)
  cat("Coordinate ranges (km):\n")
  cat("  X: [", coord_range[1,1], ", ", coord_range[2,1], "]\n")
  cat("  Y: [", coord_range[1,2], ", ", coord_range[2,2], "]\n")
  
  # Display data structure
  cat("\nData Structure Summary:\n")
  str(data)
  
  # Validate required components
  required_components <- c("y", "occ.covs", "det.covs", "coords")
  missing_components <- required_components[!required_components %in% names(data)]
  if(length(missing_components) > 0) {
    stop(paste("Missing required data components:", paste(missing_components, collapse=", ")))
  }
  cat(" All required data components present:", paste(required_components, collapse = ", "), "\n")
  
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
  
  # Define new species ordering (ecological rationale)
  sp.ordered <- c(
    "Rattus",        # Most common
    "Dasyurus",      # Carnivore
    "Vombatus",      # Large herbivore
    "Tachyglossus",  # Monotreme     
    "Canis",         # Canid
    "Felis",         # Feline
    "Vulpes",        # Fox
    "Perameles",     # Bandicoot
    "Antechinus",    # Small marsupial
    "Trichosurus",   # Possum
    "Oryctolagus",   # Rabbit
    "Wallabia"       # Wallaby
  )
  
  cat("\nNew species ordering (ecological rationale):\n")
  for(i in 1:length(sp.ordered)) {
    if(sp.ordered[i] %in% sp.names) {
      cat(sprintf("  %2d. %s (prob: %.3f)\n", i, sp.ordered[i], occ_probs[sp.ordered[i]]))
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
  cat(" Species detection data reordered successfully\n")
  
  # STEP 5: Spatial Blocking CROSS-VALIDATION SETUP ==================
  cat("\n=== STEP 5: Spatial Blocking CROSS-VALIDATION SETUP ===\n")
  
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
  
  # Create spatial blocks using Spatial Blocking
  cat("\nCreating spatial blocks using Spatial Blocking...\n")
  cat("  Method: cv_cluster (Spatial Blocking)\n")
  cat("  Number of folds:", K_FOLDS, "\n")
  cat("  Scaling: TRUE (standardize environmental variables)\n")
  
  # Perform clustering
  set.seed(123)  # For reproducibility
  sb <- cv_cluster(x = coords_sf,
                   # r = env_stack_km,
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
  
  # Calculate spatial characteristics
  cat("\n  Spatial characteristics:\n")
  cat("    - Clustering method: K-means on environmental space\n")
  cat("    - Environmental variables used:", paste(names(env_stack_km), collapse = ", "), "\n")
  
  # Save CV plot
  cat("\nGenerating spatial block visualization...\n")
  diagnostics_dir <- paste0("output/msPGOcc/diagnostics/", CV_TYPE, "_", K_FOLDS, "fold")
  dir.create(diagnostics_dir, showWarnings = FALSE, recursive = TRUE)
  plot_file <- file.path(diagnostics_dir, paste0("spatial_blocks_", K_FOLDS, "fold_", format(Sys.Date(), "%Y%m%d"), ".png"))
  
  png(plot_file, width = 1200, height = 800, res = 150)
  cv_plot(cv = sb, 
          x = coords_sf,
          r = env_stack_km)
  dev.off()
  
  cat(" Spatial block plot saved to:", plot_file, "\n")
  
  # STEP 6: SPATIAL ANALYSIS (for context) ==================================
  cat("\n=== STEP 6: SPATIAL ANALYSIS (for context) ===\n")
  
  # Calculate spatial distances
  cat("Calculating pairwise distances between sites...\n")
  dist.data <- try(dist(data.ordered$coords))
  if(!inherits(dist.data, "try-error")) {
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
    cat("  Note: This is a non-spatial model; spatial structure used for CV only\n")
  }
  
  # STEP 7: MODEL FORMULATION ===============================================
  cat("\n=== STEP 7: MODEL FORMULATION ===\n")
  
  # Define model formulas
  occ.formula <- ~ scale(bio5) + I(scale(bio5)^2) + 
    scale(bio6) + I(scale(bio6)^2) + 
    scale(bio15) +
    scale(dem) + I(scale(dem)^2) +
    scale(tpi) + I(scale(tpi)^2) +
    scale(fpc) + I(scale(fpc)^2)
  
  det.formula <- ~ scale(effort) + (1|project)
  
  cat("Model Formulas:\n")
  cat("Occupancy formula:\n")
  cat("  ", deparse(occ.formula), "\n")
  cat("Detection formula:\n")
  cat("  ", deparse(det.formula), "\n")
  
  # Count covariates
  occ_terms <- attr(terms(occ.formula), "term.labels")
  det_terms <- attr(terms(det.formula), "term.labels")
  cat("\nModel complexity:\n")
  cat("  - Occupancy covariates:", length(occ_terms), "\n")
  cat("  - Detection covariates:", length(det_terms), "\n")
  cat("  - Total species modeled:", n.species, "\n")
  
  # STEP 8: INITIAL VALUES AND PRIORS ======================================
  cat("\n=== STEP 8: INITIAL VALUES AND PRIORS ===\n")
  
  # Set initial values
  cat("Setting initial values...\n")
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
  cat("  - Community intercepts: 0\n")
  cat("  - Species-specific effects: 0\n")
  cat("  - Variance parameters: 1\n")
  cat("  - Occupancy states: Initialized from detection data\n")
  cat("  - Initial z matrix dimensions:", dim(z_inits)[1], "species x", dim(z_inits)[2], "sites\n")
  
  # Set priors (matching spatial model scripts)
  cat("\nSetting prior distributions...\n")
  priors <- list(
    beta.comm.normal = list(mean = 0, var = 2.72),    # Weakly informative
    alpha.comm.normal = list(mean = 0, var = 2.72),   # Weakly informative
    tau.sq.beta.ig = list(a = 0.1, b = 0.1),         # Inverse gamma
    tau.sq.alpha.ig = list(a = 0.1, b = 0.1)         # Inverse gamma
  )
  
  cat("Prior distributions:\n")
  cat("  - Community effects (beta, alpha): Normal(0, 2.72)\n")
  cat("  - Variance parameters (tau.sq): InverseGamma(0.1, 0.1)\n")
  
  # STEP 9: CROSS-VALIDATION MODEL FITTING =================================
  cat("\n=== STEP 9: Spatial Blocking CROSS-VALIDATION MODEL FITTING ===\n")
  
  # Prepare output directory
  output_dir <- "models/msPGOcc"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  cat("Output directory:", output_dir, "\n")
  
  # Configure reporting
  verbose <- TRUE
  n.report <- 250
  
  # Set CV threads based on number of folds
  cv.threads <- min(K_FOLDS, 8)  # Use up to 8 threads or number of folds, whichever is smaller
  
  cat("\nCross-validation configuration:\n")
  cat("  - K-fold cross-validation:", K_FOLDS, "folds\n")
  cat("  - Fold assignment: Spatial Blocking (spatial blocks)\n")
  cat("  - CV threads:", cv.threads, "\n")
  cat("  - CV seed: 123\n")
  cat("  - Mode: k.fold.only = TRUE (CV predictions only)\n")
  
  cat("\nModel fitting configuration:\n")
  cat("  - Progress reporting every:", n.report, "iterations\n")
  cat("  - Verbose output:", verbose, "\n")
  
  # Record start time
  start_time <- Sys.time()
  cat("\nModel fitting started at:", format(start_time), "\n")
  total_iterations <- n.samples * n.chains
  estimated_minutes <- round(total_iterations / 400)
  cat("Estimated runtime: ~", estimated_minutes, "minutes\n")
  cat("Total iterations across all chains:", format(total_iterations, big.mark = ","), "\n")
  
  # Fit the non-spatial multi-species model with Spatial Blocking CV
  cat("\n--- NON-SPATIAL", K_FOLDS, "FOLD Spatial Blocking CV MODEL FITTING IN PROGRESS ---\n")
  cat("Each of the", K_FOLDS, "folds will be held out and predicted using the remaining data...\n\n")
  
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
      n.chains = n.chains,
      k.fold = K_FOLDS,              # Using the configurable parameter
      k.fold.threads = cv.threads,   # Adjusted based on fold number
      k.fold.seed = 123,
      k.fold.only = TRUE,
      custom.cluster = sb$folds_ids  # Spatial Blocking fold assignments
    )
  }, error = function(e) {
    cat("\nÂ  ERROR in non-spatial multi-species model fitting:", e$message, "\n")
    return(NULL)
  }, warning = function(w) {
    cat("\nÂ  WARNING in non-spatial multi-species model fitting:", w$message, "\n")
    invokeRestart("muffleWarning")
  })
  
  # Calculate runtime
  end_time <- Sys.time()
  runtime <- end_time - start_time
  cat("\n--- NON-SPATIAL", K_FOLDS, "FOLD Spatial Blocking CV MODEL FITTING COMPLETED ---\n")
  cat("Start time:", format(start_time), "\n")
  cat("End time:", format(end_time), "\n")
  cat("Total runtime:", format(runtime), "\n")
  cat("Runtime per 1000 iterations:", format(runtime / (total_iterations/1000)), "\n")
  
  # STEP 10: MODEL SAVING AND RESULTS SUMMARY ==============================
  cat("\n=== STEP 10: MODEL SAVING AND RESULTS SUMMARY ===\n")
  
  if(!is.null(model_result)) {
    # Assign result
    out.msPGOcc <- model_result
    
    # Generate output filename with execution date at the beginning
    output_file <- file.path(output_dir, paste0(
      format(start_time, "%Y%m%d"), "_", CV_TYPE, "_", K_FOLDS, "fold_model_msPGOcc_2010-2024",
      "_nthin", n.thin,
      "_nburn", n.burn,
      "_nchain", n.chains,
      "_nsample", n.samples,
      ".RData"
    ))
    
    # Save model
    tryCatch({
      save(out.msPGOcc, file = output_file)
      cat(" Cross-validation model saved successfully to:\n")
      cat("  ", output_file, "\n")
      
      # Display file size
      file_size <- file.size(output_file)
      cat("  File size:", round(file_size / (1024^2), 2), "MB\n")
      
    }, error = function(e) {
      cat("Â  ERROR saving model:", e$message, "\n")
      # Try backup location
      backup_file <- file.path(tempdir(), basename(output_file))
      cat("Attempting backup save to:", backup_file, "\n")
      save(out.msPGOcc, file = backup_file)
      cat(" Model saved to backup location\n")
    })
    
    # Display model summary
    cat("\nNon-Spatial Multi-Species Model Summary:\n")
    cat("  - Model type: msPGOcc (Non-Spatial Multi-Species Occupancy)\n")
    cat("  - Cross-validation:", K_FOLDS, "fold Spatial Blocking\n")
    cat("  - Species fitted:", n.species, "\n")
    cat("  - Sites analyzed:", n.sites, "\n")
    cat("  - Effective samples per chain:", effective_samples, "\n")
    cat("  - Total effective samples:", effective_samples * n.chains, "\n")
    cat("  - Computational time:", format(runtime), "\n")
    
    # Cross-validation results
    if("k.fold.deviance" %in% names(out.msPGOcc)) {
      cat("\n", K_FOLDS, "Fold Cross-validation Results:\n", sep = "")
      cat("  - Mean deviance:", round(mean(out.msPGOcc$k.fold.deviance), 2), "\n")
      cat("  - SD deviance:", round(sd(out.msPGOcc$k.fold.deviance), 2), "\n")
      cat("  - Fold deviances:\n")
      for(i in 1:length(out.msPGOcc$k.fold.deviance)) {
        cat(sprintf("    Fold %2d: %.2f\n", i, out.msPGOcc$k.fold.deviance[i]))
      }
    }
    
    cat("\n NON-SPATIAL", K_FOLDS, "FOLD Spatial Blocking CV MODEL FITTING SUCCESSFUL!\n")
   
    
  } else {
    cat("Â  NON-SPATIAL", K_FOLDS, "FOLD Spatial Blocking CV MODEL FITTING FAILED\n")
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
    paste("Model Type: msPGOcc (Non-Spatial Multi-Species with", K_FOLDS, "fold", CV_TYPE, ")"),
    paste("Error:", e$message),
    "Call stack:",
    paste(capture.output(print(sys.calls())), collapse = "\n")
  ), error_file)
  
  cat("Error logged to:", error_file, "\n")
  quit(status = 1)
})

# STEP 11: PERFORMANCE METRICS AND DIAGNOSTICS ================================
cat("\n=== STEP 11: PERFORMANCE METRICS AND DIAGNOSTICS ===\n")

if(exists("out.msPGOcc") && !is.null(out.msPGOcc)) {
  
  # Create diagnostic directory
  diagnostics_dir <- paste0("output/msPGOcc/diagnostics/", CV_TYPE, "_", K_FOLDS, "fold")
  dir.create(diagnostics_dir, showWarnings = FALSE, recursive = TRUE)
  cat("Generating diagnostic summaries...\n")
  
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
  
  # Save species ordering summary
  species_summary <- data.frame(
    species = sp.ordered,
    order = 1:length(sp.ordered),
    raw_occurrence = occ_probs[sp.ordered]
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
  cat("=== ", K_FOLDS, " FOLD Spatial Blocking CV DIAGNOSTIC SUMMARY ===\n", sep = "")
  cat("Generated:", format(Sys.time()), "\n\n")
  
  cat("Model Configuration:\n")
  cat("  - Model type: msPGOcc (Non-Spatial)\n")
  cat("  - K-fold CV:", K_FOLDS, "folds\n")
  cat("  - CV Type: Spatial Blocking\n")
  cat("  - Total sites:", n.sites, "\n")
  cat("  - Sites per fold: ~", round(n.sites/K_FOLDS), "\n")
  cat("  - Species:", n.species, "\n")
  cat("  - MCMC samples:", n.samples, "\n")
  cat("  - Chains:", n.chains, "\n")
  cat("  - Thinning:", n.thin, "\n")
  cat("  - Runtime:", format(runtime), "\n\n")
  
  cat("Spatial Blocking Details:\n")
  cat("  - Method: K-means clustering in environmental space\n")
  cat("  - Environmental variables used:\n")
  for(var in names(env_stack_km)) {
    cat("    -", var, "\n")
  }
  cat("\n")
  
  if("k.fold.deviance" %in% names(out.msPGOcc)) {
    cat("Cross-Validation Performance:\n")
    cat("  - Mean deviance:", round(mean(out.msPGOcc$k.fold.deviance), 2), "\n")
    cat("  - SD deviance:", round(sd(out.msPGOcc$k.fold.deviance), 2), "\n")
    cat("  - Min deviance:", round(min(out.msPGOcc$k.fold.deviance), 2), "\n")
    cat("  - Max deviance:", round(max(out.msPGOcc$k.fold.deviance), 2), "\n\n")
    
    cat("Fold-specific deviances:\n")
    for(i in 1:length(out.msPGOcc$k.fold.deviance)) {
      cat(sprintf("  Fold %2d: %.2f\n", i, out.msPGOcc$k.fold.deviance[i]))
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

cat("\n=== ", K_FOLDS, " FOLD Spatial Blocking CV SCRIPT EXECUTION COMPLETED ===\n", sep = "")
cat("Session info:\n")
print(sessionInfo())