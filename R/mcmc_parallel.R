


mcmc_parallel <- function(cl, model_code, model_constants, model_data, model_inits,
                          params_check, monitors_add = NULL, custom_samplers = NULL){

  require(foreach)
  require(coda)
  require(doParallel)

  single_mcmc_chain <- function(){
    require(nimble)
    source("R/functions_nimble.R")

    Rmodel <- nimbleModel(
      code = model_code,
      constants = model_constants,
      data = model_data,
      inits = init,
      calculate = FALSE
    )

    Rmodel$initializeInfo()

    N <- Rmodel$N
    nH_p <- model_constants$nH_p
    n_survey <- model_constants$n_survey
    y_sum <- model_constants$y_sum

    for(i in 1:n_survey){
      N_model <- N[nH_p[i]]
      n <- round(N_model - y_sum[i])
      if(n <= 0){
        print(i)
        n <- ifelse(n == 0, 2, n)
        Rmodel$N[nH_p[i]] <- N_model + n^2
      }
    }

    calc <- Rmodel$calculate()
    if(is.infinite(calc) | is.nan(calc) | is.na(calc)){
      stop(paste0("Model log probability is ", calc))
    }

    # default MCMC configuration
    mcmcConf <- configureMCMC(Rmodel,
                              useConjugacy = TRUE)


    if(!is.null(custom_samplers)){
      for(i in seq_len(nrow(custom_samplers))){
        node <- custom_samplers$node[i]
        type <- custom_samplers$type[i]
        mcmcConf$removeSampler(node)
        mcmcConf$addSampler(node, type)
      }
    }

    n_method <- model_constants$n_method
    m_p <- model_constants$m_p
    for(i in 1:n_method){
      node <- paste0("beta_p[", i, ", ", 1:m_p, "]")
      node <- c(paste0("beta1[", i, "]"), node)
      mcmcConf$removeSampler(node)
      mcmcConf$addSampler(node, "AF_slice")
    }

    if(!is.null(monitors_add)){
      mcmcConf$addMonitors(monitors_add)
    }

    Rmcmc <- buildMCMC(mcmcConf)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <<- compileNimble(Rmcmc)
    Cmcmc$run(niter = 1000, nburnin = 0)
    samples <- as.matrix(Cmcmc$mvSamples)
    return(samples)

  }

  continue_sampling <- function(){
    require(nimble)
    require(coda)
    Cmcmc$run(1000, reset = FALSE, resetMV = TRUE)
    samples <<- as.matrix(Cmcmc$mvSamples)
    return(samples)
  }

  subset_mcmc <- function(node){
    require(dplyr)
    s <- samples[,grep(node, colnames(samples), value = TRUE, fixed = TRUE)] |>
      as.matrix() |>
      as_tibble()
    if(ncol(s) == 1) colnames(s) <- node
    return(s)
  }

  subset_params <- function(){
    require(dplyr)
    require(purrr)
    map_dfc(lapply(params_check, subset_mcmc), as_tibble) |> as.matrix()
  }

  subset_N_observed <- function(){
    map_dfc(lapply(paste0("N[", observed_state_nodes, "]"), subset_mcmc), as_tibble) |> as.matrix()
  }

  subset_N_unobserved <- function(){
    N_unobserved <- map_dfc(lapply(paste0("N[", unobserved_state_nodes, "]"), subset_mcmc), as_tibble) |> as.matrix()
    t(apply(N_unobserved, 2, quantile, c(0.025, 0.5, 0.975)))
  }

  observed_state_nodes <- model_constants$N_full_unique
  unobserved_state_nodes <- model_constants$N_quant_unique

  export <- c(
    "single_mcmc_chain",
    "continue_sampling",
    "subset_mcmc",
    "subset_params",
    "subset_N_observed",
    "subset_N_unobserved",
    "model_code",
    "model_data",
    "model_constants",
    "custom_samplers",
    "monitors_add",
    "params_check",
    "observed_state_nodes",
    "unobserved_state_nodes"
  )

  clusterExport(cl, export, envir = environment())

  for(i in seq_along(cl)){
    set.seed(i)
    init <- model_inits[[i]]
    clusterExport(cl[i], "init", envir = environment())
  }

  # initialize model and first samples
  c <- 1
  start <- Sys.time()
  out <- clusterEvalQ(cl, single_mcmc_chain())
  message("Model compile and initial ", 1000, " iterations completed in:")
  print(Sys.time() - start)

  start2 <- Sys.time()
  out2 <- clusterEvalQ(cl, continue_sampling())
  message("Additional ", 5000, " iterations completed in:")
  print(Sys.time() - start2)

  message("\nTotal iterations ", 5000 * c, " completed in:")
  print(Sys.time() - start)

  # use mcmc on clusters to subset parameters, observed states, and unobserved states
  params <- clusterEvalQ(cl, subset_params())
  params <- as.mcmc.list(lapply(params, as.mcmc))

  N_observed <- clusterEvalQ(cl, subset_N_observed())
  N_observed <- as.mcmc.list(lapply(N_observed, as.mcmc))

  N_unobserved <- clusterEvalQ(cl, subset_N_unobserved()) |> as.matrix()

  params_mcmc_list <- as.mcmc.list(lapply(params, as.mcmc))
  diagnostic <- continue_mcmc(params_mcmc_list, effective_size = 5000, max_psrf = 15)

  write_out <- function(p, no, nu, d, c, dest){
    out_rds <- list(params = p, N_observed = no, N_unobserved = nu, diagnostic = d)
    f <- "mcmcSamples.rds"
    if(c < 10) c <- paste0("0", c)
    ff <- paste0(c, "_", f)
    write_rds(out_rds, file.path(dest, f))
  }

  if(!dir.exists(dest)) dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  write_out(params, N_observed, N_unobserved, diagnostic, c, dest)

  continue <- !diagnostic$done
  while(continue){
    c <- c + 1
    resetMV <- TRUE
    clusterExport(cl, "resetMV", envir = environment())

    start2 <- Sys.time()
    out2 <- clusterEvalQ(cl, continue_sampling())
    message("Additional ", 5000, " iterations completed in:")
    print(Sys.time() - start2)

    # use mcmc on clusters to subset parameters, observed states, and unobserved states
    params <- clusterEvalQ(cl, subset_params())
    params <- as.mcmc.list(lapply(params, as.mcmc))

    N_observed <- clusterEvalQ(cl, subset_N_observed())
    N_observed <- as.mcmc.list(lapply(N_observed, as.mcmc))

    N_unobserved <- clusterEvalQ(cl, subset_N_unobserved()) |> as.matrix()

    params_mcmc_list <- as.mcmc.list(lapply(params, as.mcmc))
    diagnostic <- continue_mcmc(params_mcmc_list, effective_size = 5000, max_psrf = 15)

    write_out(params, N_observed, N_unobserved, diagnostic, c, dest)

    continue <- !diagnostic$done
  }


}




