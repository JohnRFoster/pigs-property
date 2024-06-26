#---------
#
# Workflow for fitting property-level Bayes model to MIS data
#
#---------

library(dplyr)
library(tidyr)
library(readr)
library(compareMCMCs)
library(foreach)
library(parallel)
library(purrr)

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
# itest <- nimble_inits(constants, data)
itest <- nimble_inits_sample(config$file_init, constants, data, buffer = 600)

test_build(modelCode, constants, data, itest)

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

monitors_add <- NULL

source("R/mcmc_configs.R")
cl <- makeCluster(length(nimbleMCMCdefs))
doParallel::registerDoParallel(cl)


message("Benchmarking MCMCs...")
out <- foreach(
  i = 1:length(cl),
  .inorder = TRUE
) %dopar% {
  library(compareMCMCs)
  library(nimble)
  source("R/functions_nimble.R")

  mcmc_name <- names(nimbleMCMCdefs)[i]
  mcmc_engine <- list(nimbleMCMCdefs[[i]])
  names(mcmc_engine) <- mcmc_name

  mcmc_out <- compareMCMCs(
    modelInfo = list(code = modelCode,
                     data = data,
                     constants = constants,
                     inits = itest),
    MCMCs = mcmc_name,
    nimbleMCMCdefs = mcmc_engine,
    monitors = params_check,
    MCMCcontrol = list(niter = 10000, burnin = 2500)
  )

  mcmc_out
}

stopCluster(cl)
message("  done")

str(out)

# MCMC efficiency for a parameter is defined as the effective sample size divided by
# the computation time in seconds
# i.e. It is the number of effectively independent samples generated per second.

message("\n\n====== Efficiency by MCMC engine ======")
get_byMCMC <- function(x) combineMetrics(x)$byMCMC
map(out, get_byMCMC) |> list_rbind() |> arrange(desc(min_efficiency)) |> as.matrix()

message("\n\n====== Least efficient node from each MCMC engine ======")
get_byParameter <- function(x) combineMetrics(x)$byParameter
map(out, get_byParameter) |>
  list_rbind() |>
  as_tibble() |>
  group_by(MCMC) |>
  filter(efficiency == min(efficiency)) |> as.matrix()


