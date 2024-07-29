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
         county = as.numeric(as.factor(county_code)))

data_final |>
  filter(farm_bill == 1) |>
  pull(propertyID) |>
  unique() |>
  length()

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
files_in_out_dir <- list.files(out_dir)
first_fit <- length(files_in_out_dir) == 0

monitors_add <- "N"
custom_samplers <- NULL

if(first_fit){ # run first fit

  informed <- FALSE

  data_for_nimble <- subset_data_for_development(data_final) |>
    mutate(primary_period = primary_period - min(primary_period) + 1)

  print_info(data_for_nimble)
  data_for_nimble |>
    group_by(method) |>
    summarise(n = n(),
              take = sum(take)) |> print()

  first_fit_properties <- unique(data_for_nimble$propertyID)
  all_properties <- unique(data_final$propertyID)
  not_fit_properties <- setdiff(all_properties, first_fit_properties)

  np <- length(first_fit_properties)
  dest_mcmc <- file.path(out_dir, paste0(1, "_mcmc"))
  dest_posterior <- file.path(out_dir, paste0(1, "_posterior"))

  fit_successfully <- tibble(
    property = first_fit_properties,
    round = 1,
    fit = TRUE,
    n_properties = np,
    dir = dest_posterior
  ) |>
    bind_rows(
      tibble(
        property = not_fit_properties,
        round = 0,
        fit = NA,
        n_properties = NA,
        dir = NA
      )
    )

  write_rds(fit_successfully, file.path(data_repo, "iterativeFitting.rds"))

  finished <- prep_and_run_mcmc(
    informed,
    NULL,
    data_repo,
    dest_mcmc,
    dest_posterior,
    data_for_nimble,
    monitors_add,
    custom_samplers)

  source("R/check_mcmc.R")

} else { # run iterative fitting

  informed <- TRUE

  fit_successfully <- read_rds(file.path(data_repo, "iterativeFitting.rds"))
  good_fits <- fit_successfully |> filter(fit)

  last_good_fit <- good_fits |>
    filter(round == max(round)) |>
    pull(dir) |>
    unique() |>
    basename()

  last_fit <- read_rds(file.path(out_dir, last_good_fit, "modelData.rds"))
  n_total_properties <- length(unique(data_final$propertyID))
  n_properties_to_fit <- n_total_properties - length(unique(last_fit$propertyID))

  start <- good_fits |>
    pull(round) |>
    max() +
    1

  for(i in start:n_properties_to_fit){

    # need to add a check so that if a property failed, we don't try to read the posterior from that fit
    # need to keep track of posterior file paths and order of properties
    # look back to last fit group and use
    fit_successfully <- read_rds(file.path(data_repo, "iterativeFitting.rds"))
    good_fits <- fit_successfully |> filter(fit)

    post_path <- good_fits |>
      filter(round == max(round)) |>
      pull(dir) |>
      unique()
    last_good_fit <- basename(post_path)

    data_last_fit <- read_rds(file.path(out_dir, last_good_fit, "modelData.rds"))
    data_for_nimble <- get_next_property(data_final, data_last_fit) |>
      mutate(primary_period = primary_period - min(primary_period) + 1)

    message("\n==========================================================")
    print_info(data_for_nimble)
    message("==========================================================\n")

    np <- length(data_for_nimble$propertyID)
    dest_mcmc <- file.path(out_dir, paste0(i, "_mcmc"))
    dest_posterior <- file.path(out_dir, paste0(i, "_posterior"))

    finished <- prep_and_run_mcmc(
      informed,
      post_path,
      data_repo,
      dest_mcmc,
      dest_posterior,
      data_for_nimble,
      monitors_add,
      custom_samplers)

    prop <- last(data_for_nimble$propertyID)

    fit_successfully <- fit_successfully |>
      mutate(round = if_else(property == prop, i, round),
             fit = if_else(property == prop, finished, fit),
             n_properties = if_else(property == prop, np, n_properties),
             dir = if_else(property == prop, dest_posterior, dir))

    write_rds(fit_successfully, file.path(data_repo, "iterativeFitting.rds"))

    if(finished) source("R/check_mcmc.R")

  }

}










