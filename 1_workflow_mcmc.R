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
source("R/prep_and_run_mcmc.R")
source("R/iterative_fitting.R")

# ===================================================
# Data ingest ----
# ===================================================

data_repo <- config$data_repo

## MIS data ----
file <- file.path(data_repo, config$file_mis)
interval <- config$interval
dev <- config$dev

# data_farm_bill <- read_csv(file.path(data_repo, "All_FB_Agreements.csv"))
# farm_bill_properties <- data_farm_bill |> pull(`AGR ID`)

data_mis <- get_data(file, interval)

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

print_info <- function(df){
  message("\n=======================================")
  fit_properties <- unique(df$agrp_prp_id)
  message("Total properties in data: ", length(fit_properties))
  message("Total counties in data: ", length(unique(df$county_code)))
  message("=======================================\n")
}

params_check <- config$params_check
out_dir <- config$out_dir
files_in_out_dir <- list.files(out_dir)
first_fit <- length(files_in_out_dir) == 0

monitors_add <- "N"
custom_samplers <- NULL

if(first_fit){ # run first fit

  informed <- FALSE

  data_for_nimble <- subset_data_for_development(
    df = data_final,
    max_length = 1000,          # maximum time series length (includes unsampled PPs)
    min_sampled_pp = 0.5,      # minimum number of sampled PPs in time series
    n_strata = 15,             # number of samples per strata (decile) of environmental covaraites
    properties_include = NULL # properties we want to make sure are in development data
  )

  print_info(data_for_nimble)
  prep_and_run_mcmc(informed, out_dir, data_for_nimble, monitors_add, custom_samplers)

} else { # run iterative fitting

  informed <- TRUE

  n_props_fit <- as.numeric(stringr::str_extract(files_in_out_dir, "\\d*(?=\\D)"))
  last_fit <- max(n_props_fit)
  n_total_properties <- length(unique(data_final$agrp_prp_id))

  n_properties_to_fit <- n_total_properties - last_fit

  for(i in seq_len(n_properties_to_fit)){

    n_props_fit <- as.numeric(stringr::str_extract(files_in_out_dir, "\\d*(?=\\D)"))
    last_fit <- max(n_props_fit)
    post_path <- file.path(out_dir, paste0(last_fit, "_posterior"))
    path_last_fit <- file.path(post_path, "modelData.rds")
    data_last_fit <- read_rds(path_last_fit)
    data_for_nimble <- get_next_property(data_final, data_last_fit)

    message("\n==========================================================")
    print_info(data_for_nimble)
    message("==========================================================\n")

    prep_and_run_mcmc(informed, post_path, data_for_nimble, monitors_add, custom_samplers)

  }

}










