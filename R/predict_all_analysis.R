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

fname <- file.path(read_dir, "allAundanceSamples_2025-10-01.rds")

message("Reading posterior predictive samples from: ", fname)
pp_samples <- read_rds(fname)

message("Collating posterior predictive samples")
abundance_sample <- pp_samples |>
	as_tibble() |>
	pivot_longer(cols = everything(), names_to = "node", values_to = "value") |>
	mutate(n_id = as.numeric(stringr::str_extract(node, "(?<=\\[)\\d*(?=\\])")))

head(abundance_sample)
