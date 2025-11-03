model {
  for (i in 1:N) {
    # Likelihood: cell length follows logistic mean curve
    Length[i] ~ dnorm(mu[i], within.prec)
    mu[i] <- A[root[i]] + (B[root[i]] - A[root[i]]) /
      (1 + exp((xmid_eff[i] - x01[i]) / scal[root[i]]))
    
    # xmid shifts by side and genotype
    xmid_eff[i] <- alpha[root[i]] +
                    gamma1 * side[i] +
                    gamma2 * geno[i] +
                    gamma3 * side[i] * geno[i]
  }
  
  # Random effects by root (level-2 variation)
  for (j in 1:n.root) {
    A[j] ~ dunif(0, 100)
    B[j] ~ dunif(0, 100)
    scal[j] ~ dnorm(mu_scal, between.prec)
    alpha[j] ~ dnorm(mu_alpha, between.prec)
  }
  
  # Hyperparameters
  mu_scal ~ dnorm(0, 1.0E-2)
  mu_alpha ~ dnorm(0, 1.0E-2)
  
  # Fixed effects (global)
  gamma1 ~ dnorm(0, 1.0E-2)   # side effect (Outer vs Inner)
  gamma2 ~ dnorm(0, 1.0E-2)   # genotype effect (WT vs MU)
  gamma3 ~ dnorm(0, 1.0E-2)   # interaction
  
  # Precision definitions
  within.prec <- pow(sigma2, -1)
  between.prec <- pow(tau2, -1)
  
  # Priors on variance
  sigma2 <- pow(sigma, 2)
  tau2 <- pow(tau, 2)
  sigma ~ dt(0, 0.01, 1)
  tau ~ dt(0, 0.01, 1)
}