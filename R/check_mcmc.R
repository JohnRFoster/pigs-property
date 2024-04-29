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
library(ggplot2)

source("R/functions_nimble.R")

config_name <- "hpc_dev"
config <- config::get(config = config_name)

out_dir <- "/lustrefs/ceah/feral-swine/property-fits"
np <- 22
np_dir <- file.path("dev", paste0(np, "_properties"))

# where mcmc chunks are stored
dest <- file.path(out_dir, np_dir)
mcmc_dirs <- list.files(dest)
mcmc_dirs <- setdiff(mcmc_dirs, "modelData.rds")
param_file_name <- "paramSamples.rds"

# use the first mcmc chunk to initialize storage for each chain
mcmc_rds <- file.path(dest, mcmc_dirs[1], param_file_name)
rds <- read_rds(mcmc_rds)
mcmc <- rds$params
n_chains <- length(mcmc)
store_mcmc <- list()

# store each chain, will append below
for(j in seq_len(n_chains)){
  store_mcmc[[j]] <- as.matrix(mcmc[[j]])
}

# read each mcmc chunk, store each chain from the chunk as a matrix
pb <- txtProgressBar(min = 2, max = length(mcmc_dirs), style = 1)
for(i in 2:length(mcmc_dirs)){
  mcmc_rds <- file.path(dest, mcmc_dirs[i], param_file_name)
  rds <- read_rds(mcmc_rds)
  mcmc <- rds$params
  for(j in seq_len(n_chains)){
    store_mcmc[[j]] <- rbind(store_mcmc[[j]], as.matrix(mcmc[[j]]))
  }
  setTxtProgressBar(pb, i)
}
close(pb)

# need to create an mcmc list object to check for convergence
params_mcmc_list <- as.mcmc.list(lapply(store_mcmc, as.mcmc))

# calculate psrf (convergence stat) and effective sample size
diagnostic <- continue_mcmc(params_mcmc_list, effective_size = 5000, max_psrf = 15)

# X11 is not on the HPC, so need to use ggplot2 to create traceplots
# which means we need to create a tibble of samples and identify chains and iterations
posterior <- tibble()
message("Create posterior tibble...")
pb <- txtProgressBar(min = 1, max = n_chains, style = 1)
for(i in seq_len(n_chains)){
  chain_i <- as.matrix(store_mcmc[[i]]) |>
    as_tibble() |>
    mutate(chain = i)
  posterior <- bind_rows(posterior, chain_i)

  setTxtProgressBar(pb, i)
}
close(pb)
message("  done")

# get a random sample of the posterior and save
draws <- sample.int(nrow(posterior), 5000, replace = TRUE)
posterior_samples <- posterior |>
  slice(draws) |>
  mutate(np = config$np)

# create dir for posterior samples and traceplots
np_dir <- paste0(np_dir, "_combined")
dest <- file.path(out_dir, np_dir)
if(!dir.exists(dest)) dir.create(dest, recursive = TRUE, showWarnings = FALSE)

write_rds(posterior_samples, file.path(dest, "posteriorSamples.rds"))

# function to create traceplots from thinned posterior
trace_plot <- function(post, nodes_2_plot, thin = 5000){
  df <- post |>
    select(chain, all_of(nodes_2_plot)) |>
    group_by(chain) |>
    mutate(iteration = 1:n()) |>
    ungroup()

  posterior_mat <- df |>
    select(-chain, -iteration) |>
    as.matrix()

  print(apply(posterior_mat, 2, quantile, c(0.025, 0.5, 0.975)))

  total_iterations <- max(df$iteration)
  thin_interval <- floor(seq(1, total_iterations, length.out = thin))

  gg <- df |>
    filter(iteration %in% thin_interval) |>
    pivot_longer(cols = -c(iteration, chain),
                 names_to = "node") |>
    mutate(chain = as.character(chain)) |>
    ggplot() +
    aes(x = iteration, y = value, color = chain) +
    geom_line() +
    facet_wrap(~ node, scales = "free_y") +
    labs(x = "Iteration",
         y = "Value") +
    theme_bw()
  return(gg)
}

nodes <- setdiff(colnames(posterior), "chain")
n_plots_per_page <- 4   # want to put 4 nodes on a single plot
idx <- rep(seq(1, ceiling(length(nodes) / n_plots_per_page), by = 1),
           each = n_plots_per_page)[1:length(nodes)]

plots <- tibble(
  nodes = nodes,
  idx = idx
)

message("Creating traceplots...")
# options(bitmapType = 'Xlib')
pb <- txtProgressBar(min = 1, max = max(plots$idx), style = 1)
for(i in seq_along(unique(plots$idx))){
  n2p <- plots |>
    filter(idx == i) |>
    pull(nodes)

  gg <- trace_plot(posterior, n2p)

  filename <- file.path(dest, paste0("mcmcTimeseries_", sprintf("%03d", i), ".pdf"))
  ggsave(filename, gg)

  setTxtProgressBar(pb, i)
}
close(pb)

message("  done")
