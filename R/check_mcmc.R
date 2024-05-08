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
source("R/functions_prep_nimble.R")

config_name <- "hpc_dev"
config <- config::get(config = config_name)

out_dir <- "/lustrefs/ceah/feral-swine/property-fits"
np <- 77
# np_dir <- file.path("dev", paste0(np, "_properties"))
np_dir <- "dev/161_properties_prior22"
# where mcmc chunks are stored
dest <- file.path(out_dir, np_dir)

# get samples
mcmc_list <- collate_mcmc_chunks(dest)
params_mcmc_list <- mcmc_list$params
states_mcmc_list <- mcmc_list$states

mcmc_matrix <- as.matrix(states_mcmc_list)
draws <- sample.int(nrow(mcmc_matrix), 5000, replace = TRUE)
abundance_sample <- mcmc_matrix[draws,] |>
  as_tibble() |>
  pivot_longer(cols = everything(),
               names_to = "node",
               values_to = "value") |>
  mutate(n_id = as.numeric(stringr::str_extract(node, "(?<=\\[)\\d*(?=\\])")))

data_path <- file.path(dest, "modelData.rds")
model_data <- read_rds(data_path)

all_pp <- create_all_primary_periods(data) |>
  select(-timestep)

message("\n\nall pp")
print(glimpse(all_pp))

property_info <- model_data |>
  select(agrp_prp_id, property, n_id, property_area_km2) |>
  left_join(all_pp)

message("\n\nproperty info")
print(glimpse(property_info))

property_match <- left_join(abundance_sample, property_info) |>
  mutate(density = value / property_area_km2)

message("\n\nproperty match")
print(glimpse(property_match))

stop("TESTING END")


quants = c(0.025, 0.05, seq(0.1, 0.9, by = 0.1), 0.95, 0.975)
mcmc_quants <- t(apply(mcmc_matrix, 2, quantile, quants))
node_names <- rownames(mcmc_quants)

abundance_quants <- mcmc_quants |>
  as_tibble() |>
  mutate(node = node_names)

mcmc_mean <- t(apply(mcmc_matrix, 2, mean)) |>
  as_tibble() |>
  pivot_longer(cols = everything(),
               names_to = "node",
               values_to = "mean")
mcmc_var <- t(apply(mcmc_matrix, 2, var)) |>
  as_tibble() |>
  pivot_longer(cols = everything(),
               names_to = "node",
               values_to = "variance")



print(mcmc_quants)
print(mcmc_mean)
print(mcmc_var)





total_iter <- nrow(params_mcmc_list[[1]])
n_chains <- length(params_mcmc_list)
burnin <- round(total_iter / 2)
params_burnin <- window(params_mcmc_list, start = burnin)

# calculate psrf (convergence stat) and effective sample size
diagnostic <- continue_mcmc(params_mcmc_list, effective_size = 5000, max_psrf = 15)

# X11 is not on the HPC, so need to use ggplot2 to create traceplots
# which means we need to create a tibble of samples and identify chains and iterations
posterior <- tibble()
message("Create posterior tibble...")
pb <- txtProgressBar(min = 1, max = n_chains, style = 1)
for(i in seq_len(n_chains)){
  chain_i <- as.matrix(params_mcmc_list[[i]]) |>
    as_tibble() |>
    mutate(chain = i)
  posterior <- bind_rows(posterior, chain_i)

  setTxtProgressBar(pb, i)
}
close(pb)
message("  done")

# get a random sample of the posterior and save
posterior_burnin <- params_burnin |>
  as.matrix()

draws <- sample.int(nrow(posterior_burnin), 5000, replace = TRUE)
posterior_samples <- posterior_burnin |>
  as_tibble() |>
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
