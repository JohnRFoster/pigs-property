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
# config_name <- "default"
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

data_mis <- get_data(file, interval, data_repo)

## observation covariates ----
file <- file.path(data_repo, config$file_land)
data_obs <- get_obs_covars(file)

## join MIS with observation covariates ----
data_join <- left_join(data_mis, data_obs,
                       by = join_by(county_code))

## filter missing states ----
data_join2 <- data_join |>
  filter(!st_name %in% c("CALIFORNIA", "ALABAMA", "ARIZONA", "ARKANSAS"))

targets::tar_assert_true(!any(is.na(data_join2$c_road_den)))
targets::tar_assert_true(!any(is.na(data_join2$c_rugged)))
targets::tar_assert_true(!any(is.na(data_join2$c_canopy)))

## join with farm bill properties ----
data_farm_bill <- read_csv(file.path(data_repo, "All_FB_Agreements_long_2024-05-30.csv"))
farm_bill_properties <- data_farm_bill |>
  rename(alws_agrprop_id = propertyID) |>
  select(-agreement_name, -property_name) |>
  mutate(farm_bill = 1)

data_final <- left_join(data_join2, farm_bill_properties) |>
  mutate(property = as.numeric(as.factor(propertyID)),
         county = as.numeric(as.factor(county_code))) |>
  filter(farm_bill == 1)

print_info <- function(df){
  fit_properties <- length(unique(df$propertyID))
  nfb <- df |>
    filter(farm_bill == 1) |>
    pull(propertyID) |>
    unique() |>
    length()
  message("=======================================")
  message("Total properties in data: ", fit_properties)
  message("Total Farm Bill properties in data: ", nfb)
  message("Total counties in data: ", length(unique(df$county_code)))
  message("=======================================")
}

print_info(data_final)

params_check <- config$params_check
out_dir <- config$out_dir
dest_mcmc <- file.path(out_dir, paste0(1, "_mcmc"))
dest_posterior <- file.path(out_dir, paste0(1, "_posterior"))
monitors_add <- "N"
custom_samplers <- NULL


finished <- prep_and_run_mcmc(
  informed = FALSE,
  post_path = NULL,
  data_repo = data_repo,
  dest_mcmc = dest_mcmc,
  dest_posterior = dest_posterior,
  df = data_final,
  monitors_add = monitors_add,
  custom_samplers = custom_samplers)

source("R/check_mcmc.R")









