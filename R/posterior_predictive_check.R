library(dplyr)
library(tidyr)
library(readr)
library(nimble)

set.seed(1256)

source("R/functions_prep_nimble.R")
source("R/nimble_removal_model.R")
source("R/functions_nimble.R")

config_name <- "hpc_dev"
config <- config::get(config = config_name)

out_dir <- config$out_dir
post_dir <- "1_posterior"
read_dir <- file.path(out_dir, post_dir)
fname <- "modelData.rds"

model_data <- read_rds(file.path(read_dir, fname))

fname <- "posteriorSamples.rds"
posterior_parameters <- read_rds(file.path(read_dir, fname))

fname <- "stateSamples.rds"
posterior_states <- read_rds(file.path(read_dir, fname))

posterior_samples <- bind_cols(
	posterior_parameters,
	posterior_states
)

constants <- nimble_constants(
	df = model_data,
	interval = config$interval
)

data <- nimble_data(model_data)
inits <- nimble_inits(constants = constants, data = data)

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
mcmcConf <- configureMCMC(Rmodel, useConjugacy = TRUE)

control_rw <- list(
	adaptInterval = 100,
	adaptFactorExponent = 0.6
)

mcmcConf$removeSamplers("beta1")
mcmcConf$addSampler(
	target = "beta1",
	type = "RW_block",
	control = control_rw
)

mcmcConf$removeSamplers("log_nu")
mcmcConf$addSampler(
	target = "log_nu",
	type = "slice"
)

mcmcConf$addMonitors("N")

Rmcmc <- buildMCMC(mcmcConf)
Cmcmc <- compileNimble(Rmcmc)

# don't need these samples but run to get chains moving
samples <- runMCMC(Cmcmc, niter = 50000, nburnin = 5000)

nodes <- colnames(posterior_samples)

pp_samples <- matrix(NA, n_samp, nrow(data_for_nimble))

pb <- txtProgressBar(max = n_samp, style = 1)
for (i in 1:n_samp) {
	for (j in seq_along(nodes)) {
		Cmodel[[nodes[j]]] <- posterior_samples[i, nodes[j]]
	}
	Cmodel$simulate(simNodes, includeData = TRUE)
	pp_samples[i, ] <- Cmodel[["y"]]
	setTxtProgressBar(pb, i)
}
close(pb)

write_rds(
	pp_samples,
	file.path(read_dir, "posterior_predictive_samples.rds")
)
