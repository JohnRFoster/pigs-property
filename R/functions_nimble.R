## model and helper functions for nimble
library(nimble)


# calculate the potential area (on the log scale) for a survey given method
calc_log_potential_area <- nimbleFunction(
  run = function(
    log_rho = double(1),
    log_gamma = double(1),
    p_unique = double(1),
    log_effort_per = double(0),
    effort_per = double(0),
    n_trap_m1 = double(0),
    log_pi = double(0),
    method = double(0)
  ){
    m <- method

    if(m <= 3){ # firearms, fixed wing, and helicopter
      log_potential_area <- log_rho[m] + log_effort_per
    } else { # traps and snares
      log_potential_area <- log_pi +
        (2 * (log_rho[m] + log_effort_per -
                log(exp(log_gamma[m-3]) + effort_per))) +
        log(1 + (p_unique[m-3] * n_trap_m1))
    }
    return(log_potential_area)
    returnType(double(0))
  },
  buildDerivs = TRUE
)


# figure out if we need to keep sampling
continue_mcmc <- function(mcmc, effective_size, max_psrf){
  require(coda)
  require(purrr)

  message("Checking convergence and sample size")
  psrf <- gelman.diag(mcmc, multivariate = FALSE)$psrf
  effective_samples <- effectiveSize(mcmc)

  converged <- all(psrf[, 2] < 1.1)
  enough_samples <- all(effective_samples >= effective_size)
  funky <- any(is.nan(psrf)) | max(psrf) > max_psrf

  message('Convergence [', converged, ']')
  message('Enough effective samples [', enough_samples, ']')
  message('Mixing [', !funky, ']')

  done <- converged & enough_samples
  if(done) message("MCMC complete!")

  # TODO determine burnin and effective sample size POST burnin

  if(funky){
    message("\n*** Something is wrong with the mcmc! ***")
    done <- FALSE
  }

  if(!done | funky) print(psrf)

  return(list(
    done = done,
    psrf = psrf
  ))

}

test_build <- function(code, constants, data, inits){

  message("\n\n=== Test build ===")

  Rmodel <- nimbleModel(
    code = code,
    constants = constants,
    data = data,
    inits = inits
  )

  Rmodel$initializeInfo()

  N <- Rmodel$N
  nH_p <- constants$nH_p
  n_survey <- constants$n_survey
  y_sum <- data$y_sum

  for(i in 1:n_survey){
    N_model <- N[nH_p[i]]
    n <- round(N_model - y_sum[i])
    if(n <= 0){
      print(i)
      Rmodel$N[nH_p[i]] <- N_model + abs(n)^2
    }
  }

  calc <- Rmodel$calculate()
  if(is.infinite(calc) | is.nan(calc) | is.na(calc)){
    stop(paste0("Model log probability is ", calc))
  }

  message("==================\n")
}
