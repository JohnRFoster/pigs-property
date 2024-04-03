## functions called throughout the target pipeline

# ==========================================
# Helper functions ----
# ==========================================

# need all timesteps whether there are observations or not
create_all_primary_periods <- function(df){
  pp_min_max <- df |>
    select(property, primary_period) |>
    distinct() |>
    group_by(property) |>
    filter(primary_period == min(primary_period) |
             primary_period == max(primary_period)) |>
    ungroup()

  properties <- unique(pp_min_max$property)

  all_pp <- tibble()
  message("Include all primary periods")
  pb <- txtProgressBar(max = length(properties), style = 1)
  for(i in seq_along(properties)){
    pid <- pp_min_max |> filter(property == properties[i])
    p_min <- min(pid$primary_period)
    p_max <- max(pid$primary_period)
    pp <- tibble(
      property = properties[i],
      primary_period = p_min:p_max,
      timestep = 1:length(p_min:p_max)
    )
    all_pp <- bind_rows(all_pp, pp)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  all_pp
}

# need to know the total number of timesteps in each property (sampled or not) for indexing
n_timesteps <- function(df){
  df |>
    group_by(property) |>
    filter(timestep == max(timestep)) |>
    pull(timestep)
}

# the actual primary period sequence for each property in wide format
create_wide_seq <- function(df){
  df |>
    pivot_wider(names_from = timestep,
                values_from = primary_period) |>
    select(-property) |>
    as.matrix()
}

# index (as a matrix) for tracking abundance in long format (converting wide to long)
# used in the process model
N_lookup_table <- function(df){
  df |>
    mutate(n_id = 1:n()) |>
    select(-primary_period) |>
    pivot_wider(names_from = timestep,
                values_from = n_id) |>
    select(-property) |>
    as.matrix()
}

# calculate the cumulative number of pigs taken as a primary period progresses
removed_in_pp_cumsum <- function(df){
  df |>
    group_by(property, primary_period) |>
    mutate(ysum = cumsum(take) - take) |>
    ungroup() |>
    pull(ysum)
}

# the total number of pigs taken in a primary period across all methods
# including periods without removals effort (equal to 0)
# wide format
total_take <- function(df_take, df_pp){
  sum_take <- df_take |>
    group_by(property, primary_period) |>
    summarise(sum_take = sum(take)) |>
    ungroup()

  left_join(df_pp, sum_take) |>
    mutate(sum_take = if_else(is.na(sum_take), 0, sum_take)) |>
    select(-primary_period) |>
    pivot_wider(names_from = timestep,
                values_from = sum_take) |>
    select(-property) |>
    as.matrix()
}

# index (long format) for which primary periods have removal effort
# and therefore included in the data model
N_lookup_data <- function(df_take, df_pp){
  tH <- df_take |>
    select(property, primary_period)

  df_pp |>
    select(property, primary_period) |>
    mutate(n_id = 1:n()) |>
    right_join(tH) |>
    pull(n_id)
}

# need start and end indicies for data model
# used for estimiating p
create_start_end <- function(df_take, df_pp){

  start <- end <- numeric(nrow(df_take))

  df <- left_join(df_take, df_pp)

  message("Creating start/end indicies")
  pb <- txtProgressBar(max = nrow(df), style = 1)
  for (i in seq_len(nrow(df))) {
    if (df$order[i] > 1) {
      idx <- which(df$county == df$county[i] &
                     df$property == df$property[i] &
                     df$timestep == df$timestep[i] &
                     df$order < df$order[i])
      start[i] <- idx[1]
      end[i] <- idx[length(idx)]
      assertthat::are_equal(idx, start[i]:end[i])
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  tibble(start = start, end = end)
}

# informed hyper parameters for beta distribution on global pig survival
create_surv_prior <- function(interval, data_repo){

  data <- read_csv(file.path(data_repo, "insitu/Vital_Rate_Data.csv"))

  data_usa <- data |>
    filter(country == "USA",
           time.period.end != "null",
           time.period.start != "null",
           !paper.ID %in% c(128, 1007, 130, 136)) |> # these papers don't have specified date ranges or are meta-analysis
    mutate(time.period.end = mdy(time.period.end),
           time.period.start = mdy(time.period.start))

  surv_data <- data_usa |>
    filter(!is.na(survival.prop)) |>
    select(
      unique.ID,
      paper.ID,
      N.hogs.in.study,
      contains("survival"),
      contains("hunting"),
      state,
      contains("time"),
      method.for.data)

  surv_mu <- surv_data |>
    mutate(weeks = as.numeric(time.period.end - time.period.start)/7,
           weeks4 = weeks / interval,
           survival.per.4week = survival.prop ^ (1/weeks4),
           logit.survival.per.4week = boot::logit(survival.per.4week)) |>
    filter(survival.per.4week > 0) |>
    mutate(scale_factor = survival.per.4week/survival.prop)

  surv_mu_summary <- surv_mu |>
    summarise(mu = mean(survival.per.4week),
              mu.logit = mean(logit.survival.per.4week))

  surv_var <- surv_data |>
    filter(survival.var.type %in% c("SD", "95% CI"))

  surv_sd <- surv_var |>
    filter(survival.var.type == "SD") |>
    mutate(sd = as.numeric(survival.var))

  surv_sd_calc <- surv_var |>
    filter(survival.var.type == "95% CI") |>
    mutate(low.CI = as.numeric(stringr::str_extract(survival.var, "[[:graph:]]*(?=\\-)")),
           high.CI = as.numeric(stringr::str_extract(survival.var, "(?<=\\-)[[:graph:]]*")),
           sd_low = (low.CI - survival.prop) / -1.96,
           sd_high = (high.CI - survival.prop) / 1.96) |>
    group_by(unique.ID) |>
    summarise(sd = max(sd_high, sd_low))

  surv_var_join <- left_join(surv_var, surv_sd_calc) |>
    filter(survival.var.type != "SD")

  scale_ids <- surv_mu |>
    select(unique.ID, scale_factor)

  surv_variance <- bind_rows(surv_var_join, surv_sd) |>
    left_join(scale_ids) |>
    mutate(variance = sd^2,
           variance.4week = variance * scale_factor^2,
           sd.4week = sqrt(variance.4week))

  surv_sd_summary <- surv_variance |>
    pull(sd.4week) |>
    mean()

  mu <- surv_mu_summary$mu
  psi <- 1 / mean(surv_variance$variance.4week)
  alpha <- mu * psi
  beta <- (1 - mu) * psi

  return(list(
    alpha = alpha,
    beta = beta
  ))
}

# matrix of landscape covariates for data model
create_X <- function(df, cols = c("c_road_den", "c_rugged", "c_canopy")){
  df |>
    select(all_of(cols)) |>
    as.matrix()
}




# ==========================================
# Create lists for nimble ----
# ==========================================

nimble_constants <- function(df, data_ls, interval, data_repo){

  all_primary_periods <- create_all_primary_periods(df)
  n_time_prop <- n_timesteps(all_primary_periods)
  all_pp <- create_wide_seq(all_primary_periods)
  nH <- N_lookup_table(all_primary_periods)
  nH_p <- N_lookup_data(df, all_primary_periods)
  N_full_unique <- nH_p |> unique()
  N_quant_unique <- setdiff(seq(1, max(N_full_unique)), N_full_unique)
  y_sum <- removed_in_pp_cumsum(df)
  rem <- total_take(df, all_primary_periods)
  X <- create_X(df)
  start_end <- create_start_end(df, all_primary_periods)
  survival_prior <- create_surv_prior(interval, data_repo)

  list(
    n_survey = nrow(df),
    n_ls = length(data_ls),
    n_property = length(unique(df$property)),
    n_first_survey = length(which(df$order == 1)),
    n_not_first_survey = length(which(df$order != 1)),
    n_method = length(unique(df$method)),
    n_time_prop = n_time_prop,
    n_Nfull = length(N_full_unique),
    n_Nquant = length(N_quant_unique),
    all_pp = all_pp,
    nH = nH,
    nH_p = nH_p,
    N_full_unique = N_full_unique,
    N_quant_unique = N_quant_unique,
    y_sum = y_sum,
    rem = rem,
    log_pi = log(pi),
    first_survey = which(df$order == 1),
    not_first_survey = which(df$order != 1),
    p_property_idx = as.numeric(as.factor(df$property)),
    p_pp_idx = df$primary_period,
    m_p = ncol(X),
    start = start_end$start,
    end = start_end$end,
    method = as.numeric(as.factor(df$method)),
    pp_len = interval * 7,
    phi_mu_a = survival_prior$alpha,
    phi_mu_b = survival_prior$beta
  )
}

nimble_data <- function(df, data_ls){

  X <- create_X(df)

  list(
    y = df$take,
    J = data_ls,
    X_p = X,
    effort_per = df$effort_per,
    log_effort_per = log(df$effort_per),
    n_trap_m1 = df$trap_count - 1,
    log_survey_area_km2 = log(df$property_area_km2)
  )
}


# ==========================================
# Create inits for nimble ----
# ==========================================


nimble_inits <- function(constants_nimble, data_nimble, buffer = 250){

  with(append(constants_nimble, data_nimble), {

    beta1 <- rnorm(n_method, 0, 0.25)
    beta_p <- matrix(rnorm(m_p*n_method, 0, 0.1), n_method, m_p)
    p_mu <- rnorm(2)
    log_gamma <- log(runif(2, 0.1, 2))
    log_rho <- log(
      c(runif(1, 0.1, 5), runif(1, 50, 150), runif(1, 50, 150), runif(1, 5, 15), runif(1, 5, 15))
    )
    psi_phi <- runif(1, 2, 4)
    phi_mu <- runif(1, 0.7, 0.8)
    mean_ls <- round(runif(1, 5, 8))

    a <- phi_mu * psi_phi
    b <- (1 - phi_mu) * psi_phi
    mean_lpy <- 1
    zeta <- mean_lpy / 365 * pp_len * mean_ls
    N <- phi <- rep(NA, max(nH, na.rm = TRUE))
    n_init <- rep(NA, n_property)
    for(i in 1:n_property){
      n_init[i] <- round(exp(log_survey_area_km2[i]) * 5) + sum(rem[i, ], na.rm = TRUE) * 2
      N[nH[i, 1]] <- n_init[i]
      for(j in 2:n_time_prop[i]){
        phi[nH[i, j-1]] <- rbeta(1, a, b)
        z <- N[nH[i, j-1]] - rem[i, j-1]
        z <- max(2, z)
        lambda <- z * zeta / 2 + z * phi[nH[i, j-1]]

        N[nH[i, j]] <- rpois(1, lambda)
      }
    }

    list(
      log_lambda_1 = log(n_init + buffer),
      beta_p = beta_p,
      beta1 = beta1,
      p_mu = p_mu,
      phi_mu = phi_mu,
      psi_phi = psi_phi,
      N = N + buffer,
      log_nu = log(mean_ls),
      log_gamma = log_gamma,
      log_rho = log_rho,
      phi = phi
    )

  })
}







