#---------
#
# Fit all MIS data using posterior samples from the MCMC fit
# to all properties. Simulate from the nimble model to get
# abundance/density at all properties
#
#---------

library(dplyr)
library(tidyr)
library(readr)
library(nimble)

set.seed(838)

source("R/functions_data.R")
source("R/functions_prep_nimble.R")
source("R/nimble_removal_model.R")
source("R/functions_nimble.R")

config_name <- "hpc_dev"
#config_name <- "default"
config <- config::get(config = config_name)

# get posterior samples from MCMC fit to initial properties
out_dir <- config$out_dir
post_dir <- "1_posterior"
read_dir <- file.path(out_dir, post_dir)

fname <- "posteriorSamples.rds"
posterior_samples <- read_rds(file.path(read_dir, fname))

# get data for all properties
data_store <- "C:/Users/John.Foster/OneDrive - USDA/Desktop/fosteR/data-store/"
file <- "insitu/2025-10-01/MIS.Effort.Take.All.Methods.Daily.Events.csv"
interval <- config$interval
dev <- config$dev
data_repo <- config$data_repo
data_mis <- get_data(file.path(data_store, file), interval, data_repo)

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

constants <- nimble_constants(
	data_final,
	4,
	config$data_repo,
	FALSE,
	NULL
)

data <- nimble_data(data_final)
inits <- nimble_inits(constants, data)


Rmodel <- nimbleModel(
	code = modelCode,
	constants = constants,
	data = data,
	inits = inits,
	calculate = TRUE
)

Cmodel <- compileNimble(Rmodel)

## Ensure we have the nodes needed to simulate new datasets
dataNodes <- Rmodel$getNodeNames(dataOnly = TRUE)
parentNodes <- Rmodel$getParents(dataNodes, stochOnly = TRUE)

## Ensure we have both data nodes and deterministic intermediates (e.g., lifted nodes)
simNodes <- Rmodel$getDependencies(parentNodes, self = FALSE)

# default MCMC configuration
# mcmcConf <- configureMCMC(Rmodel, useConjugacy = TRUE)

# control_rw <- list(
# 	adaptInterval = 100,
# 	adaptFactorExponent = 0.6
# )

# mcmcConf$removeSamplers("beta1")
# mcmcConf$addSampler(
# 	target = "beta1",
# 	type = "RW_block",
# 	control = control_rw
# )

# mcmcConf$removeSamplers("log_nu")
# mcmcConf$addSampler(
# 	target = "log_nu",
# 	type = "slice"
# )

# mcmcConf$addMonitors("N")

# Rmcmc <- buildMCMC(mcmcConf)
# Cmcmc <- compileNimble(Rmcmc)

nodes <- colnames(posterior_samples)

n_samp <- nrow(posterior_samples)
n_N <- max(constants$nH_p)
pp_samples <- matrix(NA, n_samp, n_N)

sty <- if_else(config_name == "default", 3, 1)
pb <- txtProgressBar(max = n_samp, style = sty)

for (i in 1:n_samp) {
	for (j in seq_along(nodes)) {
		Cmodel[[nodes[j]]] <- posterior_samples[i, nodes[j]]
	}
	Cmodel$simulate(simNodes, includeData = TRUE)
	pp_samples[i, ] <- Cmodel[["N"]]
	setTxtProgressBar(pb, i)
}
close(pb)

write_rds(
	pp_samples,
	file.path(read_dir, "allAundanceSamples_2025-10-01.rds")
)

abundance_sample <- pp_samples |>
	as_tibble() |>
	pivot_longer(cols = everything(), names_to = "node", values_to = "value") |>
	mutate(n_id = as.numeric(stringr::str_extract(node, "(?<=\\[)\\d*(?=\\])")))

all_pp <- create_all_primary_periods(data_final) |>
	select(-timestep)

all_property <- data_final |>
	select(agrp_prp_id, property, property_area_km2) |>
	distinct() |>
	right_join(all_pp)

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
