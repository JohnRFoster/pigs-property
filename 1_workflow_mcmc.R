#---------
#
# Workflow for fitting property-level Bayes model to MIS data
#
#---------

library(dplyr)
library(tidyr)
library(readr)
library(parallel)

config_name <- "hpc_dev"
config <- config::get(config = config_name)

source("R/functions_data.R")
source("R/functions_prep_nimble.R")
source("R/functions_nimble.R")
source("R/nimble_removal_model.R")
source("R/mcmc_parallel.R")

# ===================================================
# Data ingest ----
# ===================================================

data_repo <- config$data_repo

## MIS data ----
file <- file.path(data_repo, config$file_mis)
interval <- config$interval
dev <- config$dev
n <- config$np
data_mis <- get_data(file, interval, dev, n)

## observation covariates ----
file <- file.path(data_repo, config$file_land)
data_obs <- get_obs_covars(file)

## join MIS with observation covariates ----
data_join <- left_join(data_mis, data_obs,
                       by = join_by(county_code))

## filter missing states ----
data_final <- data_join |>
  filter(!st_name %in% c("CALIFORNIA", "ALABAMA", "ARIZONA", "ARKANSAS"))

targets::tar_assert_true(!any(is.na(data_final$c_road_den)))
targets::tar_assert_true(!any(is.na(data_final$c_rugged)))
targets::tar_assert_true(!any(is.na(data_final$c_canopy)))

# mean litter size year from VerCauteren et al. 2019 pg 63
data_litter_size <- round(
  c(
    5.6, 6.1, 5.6, 6.1, 4.2, 5.0, 5.0, 6.5, 5.5, 6.8,
    5.6, 5.9, 4.9, 5.1, 4.5, 4.7, 5.3, 5.7, 7.4, 8.4,
    4.7, 4.9, 3.0, 3.0, 4.8, 4.8, 4.2, 5.4, 4.7, 5.2, 5.4
  )
)



# ===================================================
# Prepare data for NIMBLE ----
# ===================================================

constants <- nimble_constants(data_final, data_litter_size, interval, data_repo)
data <- nimble_data(data_final, data_litter_size)

inits <- list()
n_chains <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
if(is.na(n_chains)) n_chains <- 3

message(n_chains, " chains")

for(i in seq_len(n_chains)){
  set.seed(i)
  inits[[i]] <- nimble_inits(constants, data)
}

params_check <- c(
  "beta_p",
  "beta1",
  "log_gamma",
  "log_rho",
  "phi_mu",
  "psi_phi",
  "log_nu",
  "p_mu"
)

monitors_add <- "N"

# ===================================================
# Run Nimble in parallel ----
# ===================================================

cl <- makeCluster(n_chains)

out_dir <- "/lustrefs/ceah/feral-swine/property-fits"
np_dir <- paste0("dev", config$np)
dest <- file.path(out_dir, np_dir)

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

