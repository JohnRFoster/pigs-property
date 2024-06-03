#---------
#
# Workflow for fitting property-level Bayes model to MIS data
#
#---------

library(dplyr)
library(tidyr)
library(readr)
library(parallel)

config_name <- "default"
config <- config::get(config = config_name)

source("R/functions_data.R")
source("R/functions_prep_nimble.R")
source("R/functions_nimble.R")
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

if(dev){

} else {
  data_for_nimble <- data_final
}



params_check <- config$params_check
out_dir <- config$out_dir
np <- constants$n_property

files_in_out_dir <- list.files(out_dir)
first_fit <- length(files_in_out_dir) == 0

if(first_fit){
  informed <- FALSE

  data_for_nimble <- subset_data_for_development(
    df = data_final,
    max_length = 1000,          # maximum time series length (includes unsampled PPs)
    min_sampled_pp = 0.5,      # minimum number of sampled PPs in time series
    n_strata = 15,             # number of samples per strata (decile) of environmental covaraites
    properties_include = NULL # properties we want to make sure are in development data
  )

} else {
  informed <- TRUE
}

fit_properties <- unique(data_for_nimble$agrp_prp_id)
message("\nTotal properties in data: ", length(fit_properties))
message("\nTotal counties in data: ", length(unique(data_for_nimble$county_code)))
message("\nMethods in data:")
table(data_for_nimble$method)


monitors_add <- "N"
custom_samplers <- NULL
post_path <- out_dir

if(first_fit){ # run first fit
  prep_and_run_mcmc(informed, out_dir, data_for_nimble, monitors_add, custom_samplers)
} else { # run iterative fitting

  n_props_fit <- as.numeric(grep("\\d*", files_in_out_dir, value = TRUE))
  last_fit <- max(n_props_fit)
  n_total_properties <- length(unique(data_final$agrp_prp_id))

  n_properties_to_fit <- n_total_properties - last_fit

  for(i in seq_len(n_properties_to_fit)){

    n_props_fit <- as.numeric(grep("\\d*", files_in_out_dir, value = TRUE))
    last_fit <- max(n_props_fit)
    path_last_fit <- file.path(out_dir, last_fit, "modelData.rds")
    data_last_fit <- read_rds(path_last_fit)
    data_for_nimble <- get_next_property(data_final, data_last_fit)

    prep_and_run_mcmc(informed, out_dir, data_for_nimble, monitors_add, custom_samplers)

  }

}










