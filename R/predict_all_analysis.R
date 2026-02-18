library(dplyr)
library(tidyr)
library(readr)
library(nimble)

set.seed(326)

source("R/functions_data.R")
source("R/functions_prep_nimble.R")
source("R/nimble_removal_model.R")
source("R/functions_nimble.R")

config_name <- "hpc_dev"
config <- config::get(config = config_name)
out_dir <- config$out_dir
post_dir <- "1_posterior"
read_dir <- file.path(out_dir, post_dir)

message("Reading posterior predictive samples from: ", fname)
fname <- "allAundanceSamples_2025-10-01.rds"
rds_name <- file.path(read_dir, fname)
pp_samples <- read_rds(rds_name)

colnames(pp_samples) <- paste0("N[", seq_len(ncol(pp_samples)), "]")

message("Collating posterior predictive samples")
abundance_sample <- pp_samples |>
	as_tibble() |>
	pivot_longer(cols = everything(), names_to = "node", values_to = "value") |>
	mutate(n_id = as.numeric(stringr::str_extract(node, "(?<=\\[)\\d*(?=\\])")))

head(abundance_sample)

message("Getting data for all properties")
# get data for all properties
file <- "insitu/2025-10-01/MIS.Effort.Take.All.Methods.Daily.Events.csv"
interval <- config$interval
dev <- config$dev
data_repo <- config$data_repo
data_mis <- get_data(file.path(data_repo, file), interval, data_repo)

## observation covariates ----
file <- file.path(data_repo, config$file_land)
data_obs <- get_obs_covars(file)

## join MIS with observation covariates ----
data_join <- left_join(data_mis, data_obs, by = join_by(county_code))

## filter missing states ----
data_join2 <- data_join |>
	filter(!st_name %in% c("CALIFORNIA", "ALABAMA", "ARIZONA", "ARKANSAS"))

targets::tar_assert_true(!any(is.na(data_join2$c_road_den)))
targets::tar_assert_true(!any(is.na(data_join2$c_rugged)))
targets::tar_assert_true(!any(is.na(data_join2$c_canopy)))

## join with farm bill properties ----
data_farm_bill <- read_csv(file.path(
	data_repo,
	"All_FB_Agreements_long_2024-05-30.csv"
))
farm_bill_properties <- data_farm_bill |>
	rename(alws_agrprop_id = propertyID) |>
	select(-agreement_name, -property_name) |>
	mutate(farm_bill = 1)

data_final <- left_join(data_join2, farm_bill_properties) |>
	mutate(
		property = as.numeric(as.factor(propertyID)),
		county = as.numeric(as.factor(county_code))
	)

all_pp <- create_all_primary_periods(data_final) |>
	select(-timestep)

all_property <- data_final |>
	select(agrp_prp_id, property, property_area_km2) |>
	distinct() |>
	right_join(all_pp)

message("Summarise predictive samples by property")
property_match <- left_join(abundance_sample, all_property) |>
	mutate(density = value / property_area_km2) |>
	group_by(
		node,
		n_id,
		agrp_prp_id,
		property,
		primary_period,
		property_area_km2
	) |>
	summarise(
		mean = mean(density),
		variance = var(density),
		`0.025` = quantile(density, 0.025),
		`0.05` = quantile(density, 0.05),
		`0.25` = quantile(density, 0.25),
		`0.5` = quantile(density, 0.5),
		`0.75` = quantile(density, 0.75),
		`0.95` = quantile(density, 0.95),
		`0.975` = quantile(density, 0.975)
	) |>
	ungroup()

write_rds(
	property_match,
	file.path(read_dir, "allDensitySummaries_2025-10-01.rds")
)
