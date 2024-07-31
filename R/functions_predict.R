
catch_model <- function(params){

  require(nimble)

  constants$log_rho <- as.matrix(params[, grep("log_rho", colnames(params))])
  constants$log_gamma <- as.matrix(params[, grep("log_gamma", colnames(params))])
  constants$p_unique <- as.matrix(ilogit(params[, grep("p_mu", colnames(params))]))
  constants$beta1 <- as.matrix(params[, grep("beta1", colnames(params))])
  constants$beta_p <- as.matrix(params[, grep("beta_p", colnames(params))])
  constants$n_ens <- nrow(params)

  with_list <- append(data, constants)

  with(with_list, {

    log_potential_area <- matrix(NA, n_ens, n_survey)
    log_theta <- matrix(NA, n_ens, n_survey)
    p <- matrix(NA, n_ens, n_survey)

    for(i in 1:n_survey){

      if(method[i] <= 3){
        log_potential_area[, i] <- shooting_lpa(
          log_rho = log_rho[, method[i]],
          effort_per = effort_per[i]
        )
      } else {
        log_potential_area[, i] <- trap_snare_lpa(
          log_rho = log_rho[, method[i]],
          log_gamma = log_gamma[, method[i] - 3],
          p_unique = p_unique[, method[i] - 3],
          effort_per = effort_per[i],
          n_trap_m1 = n_trap_m1[i]
        )
      }

      beta_p_nodes <- paste0("beta_p[", method[i], ", ", 1:3, "]")

      log_theta[, i] <- log(
        ilogit(beta1[, method[i]] +
                 X_p[i, 1] * beta_p[, beta_p_nodes[1]] +
                 X_p[i, 2] * beta_p[, beta_p_nodes[2]] +
                 X_p[i, 3] * beta_p[, beta_p_nodes[3]]
        )
      ) +
        pmin(0, log_potential_area[, i] - log_survey_area_km2[i])

      # the probability an individual is captured on the first survey
      p[, first_survey] <- exp(log_theta[, first_survey])

      # the probability an individual is captured after the first survey
      for(e in 1:n_ens){
        for(i in 1:n_not_first_survey){
          p[e, not_first_survey[i]] <- exp(log_theta[e, start[not_first_survey[i]]] +
                                             sum(log(1 - exp(log_theta[e, start[not_first_survey[i]]:end[not_first_survey[i]]]))))
        }
      }
    }
    p
  })
}

process_model <- function(params, constants, data, N_post, p){

  phi_mu <- params |> pull(phi_mu)
  psi_phi <- params |> pull(psi_phi)

  constants$log_nu <- params |> pull(log_nu)
  constants$a_phi <- phi_mu * psi_phi
  constants$b_phi <- (1 - phi_mu) * psi_phi
  constants$n_ens <- nrow(params)
  constants$p <- p

  with_list <- append(data, constants)
  with_list <- append(with_list, N_post)

  with(with_list, {

    y_pred <- matrix(NA, n_ens, n_survey)
    N <- matrix(NA, n_ens, n_time_prop)

    zeta <- exp(log(pp_len) + log_nu - log(365))
    rem <- rem[1,]

    get_N_from_post <- function(dfn, t, n){
      low <- dfn |> slice(t) |> pull(low)
      high <- dfn |> slice(t) |> pull(high)
      round(seq(low, high, length.out = n))
    }

    idx <- which(nH_p == 1)
    y_pred[,idx] <- -1

    # need a lookup table to pull postrior N
    tH <- tibble(
      time = N_full_unique,
      row = N_post$observed_timestep
    )

    N_init <- get_N_from_post(N_post, 1, n_ens)
    N[,1] <- N_init

    for(t in 2:n_time_prop){

      phi <- rbeta(n_ens, a_phi, b_phi)
      Z <- N_init - rem[t-1]
      lambda <- Z * phi + Z * zeta / 2
      N[,t] <- rpois(n_ens, lambda)

      if(t %in% N_full_unique){

        idx <- which(nH_p == t)
        for(j in idx){
          ZZ <- pmax(0, N[,t] - y_sum[j])
          take <- rpois(n_ens, p[,j] * ZZ)
          y_pred[,j] <- pmin(take, N[,t])
        }

        tt <- tH |> filter(time == t) |> pull(row)
        N_init <- get_N_from_post(N_post, tt, n_ens)

      } else {
        N_init <- N[,t]
      }
    }

    N <- tibble(
      time = 1:length(rem),
      sum_take = rem,
      `5%` = apply(N, 2, quantile, 0.05),
      `50%` = apply(N, 2, quantile, 0.50),
      `95%` = apply(N, 2, quantile, 0.95)) |>
      left_join(tH) |>
      rename(observed_timestep = row)

    Y_rep <- tibble(take = y,
                `5%` = apply(y_pred, 2, quantile, 0.05),
                `50%` = apply(y_pred, 2, quantile, 0.50),
                `95%` = apply(y_pred, 2, quantile, 0.95)) |>
      mutate(time = nH_p,
             observed_timestep = as.numeric(as.factor(time)))

    colnames(y_pred) <- paste0("c", 1:length(nH_p))
    y_pred <- as_tibble(y_pred)

    t_lookup <- tibble(
      nH_p = nH_p,
      col = paste0("c", 1:length(nH_p))
    )

    sum_take <- tibble(
      time = 1:length(rem),
      sum_take = rem
    )

    Y_pp <- y_pred |>
      mutate(iter = 1:n()) |>
      pivot_longer(cols = -iter,
                   names_to = "col") |>
      left_join(t_lookup) |>
      group_by(iter, nH_p) |>
      summarise(sum_take = sum(value)) |>
      group_by(nH_p) |>
      summarise(`5%` = quantile(sum_take, 0.05),
                `50%` = quantile(sum_take, 0.50),
                `95%` = quantile(sum_take, 0.95)) |>
      ungroup() |>
      rename(time = nH_p) |>
      left_join(sum_take)

    list(N = N, Y_rep = Y_rep, Y_pp = Y_pp)

  })
}


