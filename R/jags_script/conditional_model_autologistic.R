model{
  #priors for among city level
  for(acp in 1:ncov_psi_among){
    psi_among[acp] ~ dlogis(0,1)
  }
  for(acp1 in 1:ncov_omega_among){
    omega_among[acp1] ~ dlogis(0,1)
  }
  # City-level random effects
  #  for all model parameters
  for(pc in 1:ncov_psi){
    psi_mu[pc] ~ dlogis(0, 1)
    psi_tau[pc] ~ dgamma(1,1)
    psi_sd[pc] <- 1 / sqrt(psi_tau[pc])
    for(city in 1:ncity){
      psi[pc, city] ~ dnorm(
        psi_mu[pc],
        psi_tau[pc]
      )
    }
  }
  for(rc in 1:ncov_rho){
    rho_mu[rc] ~ dlogis(0,1)
    rho_tau[rc] ~ dgamma(1,1)
    rho_sd[rc] <- 1 / sqrt(rho_tau[rc])
    for(city in 1:ncity){
      rho[rc, city] ~ dnorm(
        rho_mu[rc],
        rho_tau[rc]
      )
    }
  }
  for(oc in 1:ncov_omega){
    omega_mu[oc] ~ dlogis(0,1)
    omega_tau[oc] ~ dgamma(1,1)
    omega_sd[oc] <- 1 / sqrt(omega_tau[oc])
    for(city in 1:ncity){
      omega[oc, city] ~ dnorm(
        omega_mu[oc],
        omega_tau[oc]
      )
    }
  }
  for(gc1 in 1:ncov_gamma){
    gamma_mu[gc1] ~ dlogis(0,1)
    gamma_tau[gc1] ~ dgamma(1,1)
    gamma_sd[gc1] <- 1 / sqrt(gamma_tau[gc1])
    for(city in 1:ncity){
      gamma[gc1, city] ~ dnorm(
        gamma_mu[gc1],
        gamma_tau[gc1]
      )
    }
  }
  # I've coded this up in sort of a lazy way
  #  in that every city has a parameter for
  #  every season. That being said, we don't
  #  have data for every season. We could 
  #  abstract this out if we feel like
  #  this is going to create some sort
  #  of computational bottleneck, but
  #  I didn't want to spend a bunch of
  #  time on this specific problem
  #  unless we need to.
  tau_psi ~ dgamma(1, 1)
  tau_rho ~ dgamma(1, 1)
  tau_omega ~ dgamma(1, 1)
  sd_psi <- 1 / sqrt(tau_psi)
  sd_rho <- 1 / sqrt(tau_rho)
  sd_omega <- 1 / sqrt(tau_omega)
  for(ns in 1:nseason){
    psi_ranef[ns] ~ dnorm(0, tau_psi)
    rho_ranef[ns] ~ dnorm(0, tau_rho)
    omega_ranef[ns] ~ dnorm(0, tau_omega)
  }
  # Autologistic term for both psi and omega
  theta_psi_mu ~ dlogis(0, 1)
  theta_psi_tau ~ dgamma(1,1)
  theta_psi_sd <- 1 / sqrt(theta_psi_tau)
  theta_omega_mu ~ dlogis(0, 1)
  theta_omega_tau ~ dgamma(1,1)
  theta_omega_sd <- 1 / sqrt(theta_omega_tau)
  for(city in 1:ncity){
    theta_psi[city] ~ dnorm(
      theta_psi_mu,
      theta_psi_tau
    )
    theta_omega[city] ~ dnorm(
      theta_omega_mu,
      theta_omega_tau
    )
  }
  # first season occupancy, mangy or otherwise
  for(first_season in 1:nfirst_season){
    logit(psi_lp[first_season]) <- inprod(
      psi[1:ncov_psi, city_vec[first_season]], 
      psi_cov[first_season,,season_vec[first_season]]
    ) + 
      inprod(
        psi_among[1:ncov_psi_among],
        psi_among_cov[first_season,]
      ) +
      psi_ranef[sample_vec[first_season]]
    z[first_season] ~ dbern(psi_lp[first_season])
  }
  # the remaining seasons occupancy, mangy or otherwise
  for(site in (nfirst_season+1):nsite){
    logit(psi_lp[site]) <- inprod(
      psi[1:ncov_psi, city_vec[site]],
      psi_cov[site,,season_vec[site]]
    ) + 
      inprod(
        psi_among[1:ncov_psi_among],
        psi_among_cov[site,]
      ) +
      psi_ranef[sample_vec[site]] + 
      theta_psi[city_vec[site]] * z[last_sample_vec[site]]
    z[site] ~ dbern(psi_lp[site])
  }
  # all seasons detection probability, mangy or otherwise
  for(site in 1:nsite){
    logit(rho_lp[site]) <- inprod(
      rho[1:ncov_rho, city_vec[site]],
      rho_cov[site,,season_vec[site]]
    ) + 
      rho_ranef[sample_vec[site]]
    y[site] ~ dbin(rho_lp[site] * z[site], J[site])
  }
  # first season occupancy, conditional mange
  for(first_season in 1:nfirst_season){
    logit(ome_mu[first_season]) <- inprod(
      omega[1:ncov_omega, city_vec[first_season]],
      omega_cov[first_season,,season_vec[first_season]])+ 
      inprod(
        omega_among[1:ncov_omega_among],
        omega_among_cov[first_season,]
      ) + 
      omega_ranef[sample_vec[first_season]]
    x[first_season] ~ dbern(ome_mu[first_season] * z[first_season])
  }
  # the remaining seasons occupancy, conditional mange
  for(site in (nfirst_season+1):nsite){ #nsite is the number of site-seasons available, 
    logit(ome_mu[site]) <- inprod(
      omega[1:ncov_omega, city_vec[site]], 
      omega_cov[site,,season_vec[site]]
    ) + 
      inprod(
        omega_among[1:ncov_omega_among],
        omega_among_cov[site,]
      ) + 
      omega_ranef[sample_vec[site]] + #city_vec[site]] + #could this need the vector with the city id values in the mange dataset
      theta_omega[city_vec[site]] * x[last_sample_vec[site]]
    x[site] ~ dbern(ome_mu[site] * z[site])
  }
  # image specific mange detection
  for(photo in 1:nphoto){ #for each row of the mange dataset
    logit(gam_mu[photo]) <- inprod(
      gamma[1:ncov_gam, city_photo_vec1[photo]], #1 to 4 covariates gamma, for each line, get gamma[1:4, cityID of that photo]
      gamma_cov[photo,]
    ) 
    q[photo] ~ dbern(gam_mu[photo] * x[siteseasvec_coy[photo]]) #instead of site_vec changed to a site_season specific as we changed x_guess and z
    #q[1036] giving trouble
  }
  # # derived quantities
  # for(ns in 1:nseason){
  #   n_coyote[ns] <- sum(z[c(st[ns,1]:st[ns,2])])
  #   n_mange[ns] <-  sum(x[c(st[ns,1]:st[ns,2])])
  # }
}
var z[nsite], psi_lp[nsite], ome_mu[nsite], x[nsite];