
data_list <- readRDS(./data/data_list.rds)
inits <- readRDS(./data/inits.rds)

# fit the model
start <- Sys.time()
cl <- parallel::makeCluster(18)
mout <- run.jags(
  model = "./jags_script/conditional_model_autologistic.R",
  monitor = c(
    'psi',
    'rho',
    'omega',
    'gamma',
    'sd_psi',
    'sd_rho',
    'sd_omega',
    'psi_ranef',
    'rho_ranef',
    'omega_ranef',
    "theta_psi",
    "theta_omega",
    'n_coyote',
    'n_mange',
    'psi_mu',
    'rho_mu',
    'omega_mu',
    'gamma_mu',
    'psi_among',
    'omega_among',
    'theta_psi_mu',
    'theta_psi_tau',
    'theta_omega_mu',
    'theta_omega_tau'
  ),
  data = data_list,
  n.chains = 6,
  inits = inits,
  adapt = 1000, #ideal 1000
  burnin = 25000,#former 25000, #ideal 20000, #fast run: no burning
  sample = ceiling(200000/4), #former 200000/4), #ideal ceiling(75000/12) #fast run: 300000/3 (chains)
  thin = 4, #2, #try 4 next
  module = 'glm',
  # cl = cl,
  method = 'parallel' #'rjparallel',
)
end <- Sys.time()
end - start

parallel::stopCluster(cl)

m2 <- as.mcmc.list(mout)
saveRDS(mout, paste0("./results/coyote_mcmc_autologistic.RDS"))


