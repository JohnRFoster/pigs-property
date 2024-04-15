nimbleMCMCdefs <- list(

  default = function(model){
    mcmcConf <- configureMCMC(model)
    mcmcConf
  },

  # default_thetaBlock_afSlice = function(model){
  #   mcmcConf <- configureMCMC(model)
  #   for(i in 1:5){
  #     node <- paste0("beta_p[", i, ", ", 1:m_p, "]")
  #     node <- c(paste0("beta1[", i, "]"), node)
  #     mcmcConf$removeSampler(node)
  #     mcmcConf$addSampler(node, "AF_slice")
  #   }
  #   mcmcConf
  # },
  #
  # default_thetaBlock_rwBlock = function(model){
  #   mcmcConf <- configureMCMC(model)
  #   for(i in 1:5){
  #     node <- paste0("beta_p[", i, ", ", 1:m_p, "]")
  #     node <- c(paste0("beta1[", i, "]"), node)
  #     mcmcConf$removeSampler(node)
  #     mcmcConf$addSampler(node, "RW_block")
  #   }
  #   mcmcConf
  # },

  slice_scalars1 = function(model){
    mcmcConf <- configureMCMC(model)
    mcmcConf$removeSamplers(c("phi_mu", "psi_phi", "log_nu"))
    mcmcConf$addSampler(target = c("phi_mu"), type = "slice")
    mcmcConf$addSampler(target = c("psi_phi"), type = "slice")
    mcmcConf$addSampler(target = c("log_nu"), type = "slice")
    mcmcConf
  },

  slice_scalars2 = function(model){
    mcmcConf <- configureMCMC(model)
    mcmcConf$removeSamplers(c("phi_mu", "psi_phi", "log_nu"))
    mcmcConf$addSampler(target = c("phi_mu"), type = "slice")
    mcmcConf$addSampler(target = c("psi_phi"), type = "slice")
    mcmcConf$addSampler(target = c("log_nu"), type = "ess")
    mcmcConf
  },

  rwBlock_data = function(model){
    mcmcConf <- configureMCMC(model)
    mcmcConf$removeSamplers(c("log_gamma", "p_mu"))
    mcmcConf$addSampler(target = c("log_gamma"), type = "RW_block")
    mcmcConf$addSampler(target = c("p_mu"), type = "RW_block")
    mcmcConf
  },

  afSlice_data = function(model){
    mcmcConf <- configureMCMC(model)
    mcmcConf$removeSamplers(c("log_gamma", "p_mu", "log_rho"))
    mcmcConf$addSampler(target = c("log_gamma"), type = "AF_slice")
    mcmcConf$addSampler(target = c("log_rho"), type = "AF_slice")
    mcmcConf$addSampler(target = c("p_mu"), type = "AF_slice")
    mcmcConf
  },

  slice = function(model){
    mcmcConf <- configureMCMC(model)
    mcmcConf$removeSamplers(c("beta_p", "beta1", "log_gamma", "p_mu", "log_rho", "phi_mu", "psi_phi", "log_nu"))
    mcmcConf$addSampler(target = c("phi_mu"), type = "slice")
    mcmcConf$addSampler(target = c("psi_phi"), type = "slice")
    mcmcConf$addSampler(target = c("log_nu"), type = "slice")
    mcmcConf$addSampler(target = c("log_gamma"), type = "AF_slice")
    mcmcConf$addSampler(target = c("log_rho"), type = "AF_slice")
    mcmcConf$addSampler(target = c("p_mu"), type = "AF_slice")
    mcmcConf$addSampler(target = c("beta_p"), type = "AF_slice")
    mcmcConf$addSampler(target = c("beta1"), type = "AF_slice")
    mcmcConf
  }

)
