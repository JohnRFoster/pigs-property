
prep_and_run_mcmc <- function(informed, post_path, df, monitors_add, custom_samplers){

  source("R/nimble_removal_model.R")
  source("R/functions_prep_nimble.R")
  source("R/functions_nimble.R")
  source("R/mcmc_parallel.R")

  # ===================================================
  # Prepare data for NIMBLE ----
  # ===================================================

  # mean litter size year from VerCauteren et al. 2019 pg 63
  data_litter_size <- round(
    c(
      5.6, 6.1, 5.6, 6.1, 4.2, 5.0, 5.0, 6.5, 5.5, 6.8,
      5.6, 5.9, 4.9, 5.1, 4.5, 4.7, 5.3, 5.7, 7.4, 8.4,
      4.7, 4.9, 3.0, 3.0, 4.8, 4.8, 4.2, 5.4, 4.7, 5.2, 5.4
    )
  )

  constants <- nimble_constants(
    df,
    data_litter_size,
    interval,
    data_repo,
    informed,
    post_path
  )

  data <- nimble_data(df, data_litter_size)

  n_chains <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
  if(is.na(n_chains)) n_chains <- 3

  message(n_chains, " chains")

  inits <- list()
  for(i in seq_len(n_chains)){
    set.seed(i)
    inits[[i]] <- nimble_inits(constants, data)
  }

  test_build(modelCode, constants, data, inits[[1]])

  # ===================================================
  # Run Nimble in parallel ----
  # ===================================================

  np_dir <- length(unique(df$agrp_prp_id))
  dest_mcmc <- file.path(out_dir, paste0(np_dir, "_mcmc"))
  dest_posterior <- file.path(out_dir, paste0(np_dir, "_posterior"))
  if(!dir.exists(dest_mcmc)) dir.create(dest_mcmc, recursive = TRUE, showWarnings = FALSE)
  if(!dir.exists(dest_posterior)) dir.create(dest_posterior, recursive = TRUE, showWarnings = FALSE)

  write_rds(df, file.path(dest_posterior, "modelData.rds"))

  message("Begin Parallel sampling and make cluster")
  cl <- makeCluster(n_chains)

  mcmc_parallel(
    cl = cl,
    model_code = modelCode,
    model_constants = constants,
    model_data = data,
    model_inits = inits,
    params_check = params_check,
    n_iters = config$n_iter,
    dest = dest,
    monitors_add = monitors_add,
    custom_samplers = NULL
  )

  stopCluster(cl)
}

