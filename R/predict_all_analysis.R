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

samples <- read_rds(fname)

print(glimpse(samples))
