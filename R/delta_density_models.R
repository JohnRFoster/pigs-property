library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(mgcv)

source("R/functions_delta_density_models.R")

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])

cutoff_date <- ymd("2023-12-31")

config_name <- "hpc_dev"
# config_name <- "default"
config <- config::get(config = config_name)

source("R/functions_data.R")

# ===================================================
# Data ingest ----
# ===================================================

data_repo <- config$data_repo

center_scale <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

## observation covariates ----
file <- file.path(data_repo, config$file_land)
data_obs <- get_obs_covars(file) |>
  mutate(c_prop.pub.land = center_scale(prop.pub.land)) |>
  select(county_code, contains("c_"))


## join with farm bill properties ----
data_farm_bill <- read_csv(file.path(data_repo, "All_FB_Agreements_long_2024-05-30.csv"))
farm_bill_properties <- data_farm_bill |>
  rename(alws_agrprop_id = propertyID) |>
  select(-agreement_name, -property_name) |>
  mutate(farm_bill = 1)

if_dir <- "11_posterior"
top_dir <- config$out_dir

posterior_path <- file.path(top_dir, if_dir, "modelData.rds")
data <- read_rds(posterior_path) |>
  mutate(method = if_else(method == "FIREARMS", "Firearms", method),
         method = if_else(method == "FIXED WING", "Fixed wing", method),
         method = if_else(method == "HELICOPTER", "Helicopter", method),
         method = if_else(method == "SNARE", "Snare", method),
         method = if_else(method == "TRAPS", "Trap", method))

posterior_path <- file.path(top_dir, if_dir, "densitySummaries.rds")
density <- read_rds(posterior_path)

n_return <- data |>
  select(propertyID, primary_period, method) |>
  distinct() |>
  pivot_wider(names_from = method,
              values_from = method) |>
  unite(method,
        -c(propertyID, primary_period),
        sep = ", ",
        na.rm = TRUE) |>
  group_by(propertyID, method) |>
  mutate(return_interval = c(0, diff(primary_period))) |>
  ungroup() |>
  rename(methods_used = method)

data_grouped <- data |>
  group_by(propertyID, agrp_prp_id, start_dates, end_dates, st_name, cnty_name, farm_bill,
           alws_agrprop_id, property, primary_period, property_area_km2, county_code) |>
  summarise(total_take = sum(take),
            total_effort_per = sum(effort_per),
            n_events = n()) |>
  ungroup() |>
  mutate(take_density = total_take / property_area_km2,
         farm_bill = if_else(is.na(farm_bill), 0, farm_bill)) |>
  left_join(density) |>
  left_join(n_return) |>
  mutate(year = year(end_dates))

data_mis <- data_grouped |>
  filter(end_dates <= cutoff_date) |>
  mutate(abundance_estimate = round(`0.5` * property_area_km2)) |>
  rename(density_estimate = `0.5`) |>
  select(propertyID, end_dates, year, st_name, cnty_name, county_code,
         property_area_km2, total_take, take_density, density_estimate, abundance_estimate,
         n_events, methods_used)

first_flag <- -1e4

change_df <- data_mis |>
  group_by(propertyID, year, st_name, county_code, property_area_km2) |>
  summarise(med_density = median(density_estimate),
            sum_take = sum(total_take),
            avg_take = mean(total_take),
            sum_events = sum(n_events),
            avg_events = mean(n_events)) |>
  arrange(propertyID, year) |>
  group_by(propertyID) |>
  mutate(delta_density = c(first_flag, diff(med_density)),
         delta_take = c(first_flag, diff(sum_take)),
         delta_events = c(first_flag, diff(sum_events)),
         cum_take = cumsum(sum_take)) |>
  ungroup() |>
  filter(delta_density != first_flag)

assertthat::assert_that(all(change_df != first_flag))

data <- change_df |>
  left_join(data_obs) |>
  mutate(
    # delta_density = center_scale(delta_density),
    # sum_take = center_scale(sum_take),
    # cum_take = center_scale(cum_take),
    # avg_take = center_scale(avg_take),
    # sum_events = center_scale(sum_events),
    # avg_events = center_scale(avg_events),
    # delta_take = center_scale(delta_take),
    # delta_events = center_scale(delta_events),
    # property_area_km2 = center_scale(property_area_km2),
    st_name = factor(st_name),
    year = factor(year),
    county_code = factor(county_code),
    state_year = factor(paste(st_name, year)),
    county_year = factor(paste(county_code, year)))

data |>
  select(contains("take"), contains("events"), contains("c_"), contains("delta")) |>
  select(-delta_density) |>
  as.matrix() |>
  cor()

# re-coding state and county because, for example, 2020 in TX is not the same as 2020 in MO
# same goes for counties

glimpse(data)

method <- "REML"
family <- "gaussian"

# Groups that could matter
# state
# county
# year
# state-year (see above)
# county-year (see above)

# create models by group/smoother and variable criteria

# G model - A single common smoother for all observations; only has a Global smoother

# GS model - A global smoother plus group-level smoothers that have the *same* wiggliness
#   - Global smoother with individual effects that have a *Shared* penalty

# GI model - A global smoother plus group-level smoothers with *differing* wiggliness
#   - Global smoother with individual effects that have *Individual* penalties
#   - This is useful if different groups differ substantially in how wiggly they are.

# S model - Group-specific smoothers without a global smoother
#   - all smoothers having the *same* wiggliness

# I model - Group-specific smoothers with different wiggliness
#   - all smoothers having the *different* wiggliness

# T models - yearly totals

# M models - yearly means

# D models - yearly change


all_props <- unique(data$propertyID)

draw <- sample.int(length(all_props), 50)
test_props <- all_props[draw]

if(config_name == "default"){
  model_data <- data |>
    filter(propertyID %in% test_props)
} else {
  model_data <- data
}


k_sum_take <- 30
k_state_year <- length(unique(model_data$state_year))
k_county_year <- length(unique(model_data$county_year))

dest <- config$out_delta

if(task_id == 1) m0(model_data, method, family, file.path(dest, "m0.rds"))
if(task_id == 2) m1(model_data, method, family, file.path(dest, "m1.rds"))

if(task_id == 3) G_T(model_data, method, family, k_sum_take, k_state_year, file.path(dest, "G_T.rds"))
if(task_id == 4) G_M(model_data, method, family, k_sum_take, k_state_year, file.path(dest, "G_M.rds"))
if(task_id == 5) G_D(model_data, method, family, k_sum_take, k_state_year, file.path(dest, "G_D.rds"))

if(task_id == 6) GS_T(model_data, method, family, file.path(dest, "GS_T.rds"))
if(task_id == 7) GS_M(model_data, method, family, file.path(dest, "GS_M.rds"))
if(task_id == 8) GS_D(model_data, method, family, file.path(dest, "GS_D.rds"))

if(task_id == 9)  GI_T(model_data, method, family, k_state_year, file.path(dest, "GI_T.rds"))
if(task_id == 10) GI_M(model_data, method, family, k_state_year, file.path(dest, "GI_M.rds"))
if(task_id == 11) GI_D(model_data, method, family, k_state_year, file.path(dest, "GI_D.rds"))

if(task_id == 12) S_T(model_data, method, family, k_state_year, file.path(dest, "S_T.rds"))
if(task_id == 13) S_M(model_data, method, family, k_state_year, file.path(dest, "S_M.rds"))
if(task_id == 14) S_D(model_data, method, family, k_state_year, file.path(dest, "S_D.rds"))

if(task_id == 15) I_T(model_data, method, family, k_state_year, file.path(dest, "I_T.rds"))
if(task_id == 16) I_M(model_data, method, family, k_state_year, file.path(dest, "I_M.rds"))
if(task_id == 17) I_D(model_data, method, family, k_state_year, file.path(dest, "I_D.rds"))


# modGS1 <- gam(delta_density ~
#                te(avg_take, avg_events, bs = c("tp", "tp")) +  # thin plate regression spline
#                t2(avg_take, avg_events, st_name, bs = c("tp", "tp", "re")) +
#                s(state_year, bs = "re", k = k_state_year), # random effect
#              method = method,
#              family = family,
#              data = data)
#
# summary(modGS1)

# m2 <- gam(delta_density ~
#             s(sum_take) +                                    #S1
#             s(sum_take, year, bs = "sz") +                   #S2
#             s(sum_take, st_name, bs = "fs", by = year),      #S3
#           data = data,
#           method = method)

# summary(m2)
# plot(m2, pages = 1)
# gam.check(m2)

# S1 is the average smooth effect of total yearly take, over the entire data set.
# This would capture, for example, the general tendency for the change in yearly density
# to increase, reach a peak, and then decline again as take increases.

# S2 is a difference smooth, parametrised in such a way as to be orthogonal to the average
# smooth of total yearly take. This would describe the smooth differences between the "average"
# total take effect and the smooth effects for each year.

# S3 is a random factor smooth term, with a smooth of total yearly take for each state,
# which are allowed to differ by year. Modelling this more like a random effect because
# one could think of many states from which the change in density is drawn, and there are
# likely quite a few states, but they all have roughly the same wiggliness, but different
# shapes — I wouldn't expect the smooth effect of total take to be really wiggly for TX
# and really smooth OK for example.










