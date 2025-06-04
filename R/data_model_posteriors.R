library(dplyr)
library(tidyr)
library(readr)

source("R/functions_data.R")

config_name <- "hpc_dev"
config <- config::get(config = config_name)

data_repo <- config$data_repo

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])

message("Task: ", task_id)

## MIS data ----
file <- file.path(data_repo, config$file_mis)
interval <- config$interval
dev <- config$dev

# top_dir <- "C:/Users/John.Foster/Downloads/iterativeFits2"

if_dir <- "11_posterior"
posterior_path <- file.path(config$out_dir, if_dir, "posteriorSamples.rds")
params <- read_rds(posterior_path)

glimpse(params)

posterior_path <- file.path(config$out_dir, if_dir, "modelData.rds")
data <- read_rds(posterior_path) |>
  mutate(method = if_else(method == "FIREARMS", "Sharpshooting", method),
         method = if_else(method == "FIXED WING", "Fixed wing", method),
         method = if_else(method == "HELICOPTER", "Helicopter", method),
         method = if_else(method == "SNARE", "Snare", method),
         method = if_else(method == "TRAPS", "Trap", method))

glimpse(data)

effort_table <- data |>
  group_by(method) |>
  summarise(
    min_units = min(trap_count),
    max_units = max(trap_count),
    min_effort = min(effort),
    max_effort = max(effort),
    min_area = min(property_area_km2),
    max_area = max(property_area_km2)) |>
  ungroup()

glimpse(effort_table)

get_grid <- function(m){

  tmp <- effort_table |> filter(method == m)

  glimpse(tmp)

  if(m == "Snare"){
    b <- 10
  } else if(m == "Trap"){
    b <- 5
  } else {
    b <- 1
  }

  tc <- seq(tmp$min_units, tmp$max_units, by = b)

  expand_grid(
    method = m,
    trap_count = tc,
    effort = seq(tmp$min_effort, tmp$max_effort, length.out = 10),
    area = seq(2, 150, by = 16),
    c_road_den = -3:3,
    c_rugged = -2:2,
    c_canopy = -2:2
  )

}

catch_model <- function(params, effort_df){

  method = effort_df$method

  if(method == "Sharpshooting") id <- 1
  if(method == "Fixed wing") id <- 2
  if(method == "Helicopter") id <- 3
  if(method == "Snare") id <- 4
  if(method == "Trap") id <- 5

  log_rho <- params |>
    pull(paste0("log_rho[", id, "]"))

  beta1 <- params |>
    pull(paste0("beta1[", id, "]"))

  beta_p_nodes <- paste0("beta_p[", id, ", ", 1:3, "]")
  beta_p <- params |>
    select(all_of(beta_p_nodes))

  beta_p_1 <- beta_p |> pull(1)
  beta_p_2 <- beta_p |> pull(2)
  beta_p_3 <- beta_p |> pull(3)

  n_ens <- nrow(params)

  trap_count <- effort_df$trap_count
  effort <- effort_df$effort
  area <- effort_df$area
  log_e <- log(effort / trap_count)

  c_road_den <- effort_df$c_road_den
  c_rugged <- effort_df$c_rugged
  c_canopy <- effort_df$c_canopy


  if(method %in% c("Trap", "Snare")){

    log_gamma <- params |>
      pull(paste0("log_gamma[", id - 3, "]"))

    p_unique <- params |>
      pull(paste0("p_mu[", id - 3, "]")) |>
      boot::inv.logit()

    lpa <- log(pi) +
      (2 * (log_rho + log_e -
              log(exp(log_gamma) + exp(log_e)))) +
      log(1 + (p_unique * (trap_count - 1)))

  } else {
    lpa <- log_rho + log_e
  }

  log_prop_searched <- pmin(0, lpa - log(area))

  log_theta <- log(
    nimble::ilogit(beta1 +
                     c_road_den * beta_p_1 +
                     c_rugged * beta_p_2 +
                     c_canopy * beta_p_3
    )
  ) +
    log_prop_searched

  effort_df |>
    mutate(theta = median(exp(log_theta)),
           area_searched = median(exp(lpa)),
           prop_prop_searched = median(exp(log_prop_searched)))

}


tasks <- c(
  "Sharpshooting",
  "Fixed wing",
  "Helicopter",
  "Snare",
  "Trap"
)

method_task <- tasks[task_id]

hyper_grid <- get_grid(method_task)

all_out <- tibble()
pb <- txtProgressBar(min = 1, max = nrow(hyper_grid), style = 1)
for(i in seq_len(nrow(hyper_grid))){

  tmp <- catch_model(params, hyper_grid[i, ])
  all_out <- bind_rows(all_out, tmp)

  setTxtProgressBar(pb, i)

}
close(pb)

f <- paste0("data_model_", method_task, ".csv")
write_csv(all_out, file.path(config$out_eom, f))

glimpse(all_out)
