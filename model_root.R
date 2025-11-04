model {
  
  for (i in 1:N) {
    # Likelihood: cell length follows logistic mean curve
    Length[i] ~ dnorm(mu[i], within.prec)
    
    mu[i] <- A[root[i]] + (B[root[i]] - A[root[i]]) /
      (1 + exp((xmid_eff[i] - x01[i]) / scal[root[i]]))
    
    # Effective midpoint includes side and genotype effects
    xmid_eff[i] <- alpha[root[i]] +
      gamma1 * side[i] +
      gamma2 * geno[i] +
      gamma3 * side[i] * geno[i]
  }
  
  #------------------------------------------------------
  # Random effects by root (level-2 variation)
  #------------------------------------------------------
  for (j in 1:n.root) {
    A[j]     ~ dnorm(mu_A, prec_A)
    B[j]     ~ dnorm(mu_B, prec_B)
    scal[j]  ~ dnorm(mu_scal, between.prec)
    alpha[j] ~ dnorm(mu_alpha, between.prec)
  }
  
  #------------------------------------------------------
  # Hyperparameters informed by NLS
  #------------------------------------------------------
  mu_A     ~ dnorm(muA_hat,      pow(prior_sd_A,     -2))
  mu_B     ~ dnorm(muB_hat,      pow(prior_sd_B,     -2))
  mu_alpha ~ dnorm(muAlpha_hat,  pow(prior_sd_alpha, -2))
  mu_scal  ~ dnorm(muScal_hat,   pow(prior_sd_scal,  -2))
  
  prec_A ~ dgamma(0.001, 0.001)
  prec_B ~ dgamma(0.001, 0.001)

  #------------------------------------------------------
  # Fixed effects (global parameters)
  #------------------------------------------------------
  gamma1 ~ dnorm(0, 1.0E-2)   # side effect (Outer vs Inner)
  gamma2 ~ dnorm(0, 1.0E-2)   # genotype effect (WT vs MU)
  gamma3 ~ dnorm(0, 1.0E-2)   # interaction effect
  
  #------------------------------------------------------
  # Variance components and priors
  #------------------------------------------------------
  within.prec  <- pow(sigma2, -1)
  between.prec <- pow(tau2, -1)
  
  sigma2 <- pow(sigma, 2)
  tau2   <- pow(tau, 2)
  sigma ~ dunif(0, 10)
  tau   ~ dunif(0, 5)
}