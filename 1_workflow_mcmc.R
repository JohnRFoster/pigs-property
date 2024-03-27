#---------
#
# Workflow for fitting property-level Bayes model to MIS data
#
#---------

library(dplyr)
library(tidyr)
library(readr)

config_name <- "default"
config <- config::get(config = config_name)
interval <- config$interval

source("R/functions_data.R")
source("R/functions_prep_nimble.R")

# ===================================================
# Data ingest ----
# ===================================================

data_repo <- config$data_repo

## MIS data ----
file <- file.path(data_repo, config$file_mis)
data_mis <- get_data(file, interval)

## observation covariates ----
file <- file.path(data_repo, config$file_land)
data_obs <- get_obs_covars(file)

## join MIS with observation covariates ----
data_join <- left_join(data_mis, data_obs)

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

constants <- nimble_constants(data_final, data_litter_size, interval)
data <- nimble_data(data_final, data_litter_size)



