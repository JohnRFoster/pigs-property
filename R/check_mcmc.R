#---------
#
# Workflow for checking property-level Bayes model fits of the MIS data
# - combines mcmc chunks
# - checks for convergence
#
#---------

library(dplyr)
library(tidyr)
library(readr)
library(parallel)
library(coda)

source("R/functions_nimble.R")


config_name <- "default"
config <- config::get(config = config_name)

out_dir <- "/lustrefs/ceah/feral-swine/property-fits"
np_dir <- paste0("dev", config$np)
dest <- file.path(out_dir, np_dir)

mcmc_dirs <- list.files(dest)
param_file_name <- "paramSamples.rds"

mcmc_rds <- file.path(dest, mcmc_dirs[1], param_file_name)
rds <- read_rds(mcmc_rds)
mcmc <- rds$params
n_chains <- length(mcmc)
store_mcmc <- list()

for(j in seq_len(n_chains)){
  store_mcmc[[j]] <- as.matrix(mcmc[[j]])
}


for(i in 2:5){
  mcmc_rds <- file.path(dest, mcmc_dirs[i], param_file_name)
  rds <- read_rds(mcmc_rds)
  mcmc <- rds$params
  for(j in seq_len(n_chains)){
    store_mcmc[[j]] <- rbind(store_mcmc[[j]], as.matrix(mcmc[[j]]))
  }
}

params_mcmc_list <- as.mcmc.list(lapply(store_mcmc, as.mcmc))
diagnostic <- continue_mcmc(params_mcmc_list, effective_size = 5000, max_psrf = 15)

np_dir <- paste0("dev", config$np, "combined")
dest <- file.path(out_dir, np_dir)
if(!dir.exists(dest)) dir.create(dest, recursive = TRUE, showWarnings = FALSE)

message("Creating traceplots...")
pdf(file = file.path(dest, "mcmcTimeseries%03d.pdf"), onefile = FALSE)
plot(params_mcmc_list)
dev.off()
message("  done")
