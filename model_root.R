model {
  for (i in 1:N) {
    # Likelihood: cell length follows logistic mean curve
    Length[i] ~ dnorm(mu[i], within.prec)
    
    mu[i] <- A[root[i]] + (B[root[i]] - A[root[i]]) /
      (1 + exp((xmid_eff[i] - x01[i]) / scal[root[i]]))
    
    # Effective midpoint includes side and genotype effects
    # Outer: side = 1; Mutant: geno = 1
    xmid_eff[i] <- alpha[root[i]] +
      gamma1 * side[i] +
      gamma2 * geno[i] +
      gamma3 * side[i] * geno[i]
  }
  
  # Random effects by root
  for (j in 1:n.root) {
    # Latent parameters for root-level variation
    a[j]      ~ dnorm(mu_a,  prec_a)         # log(A_j)
    d[j]      ~ dnorm(mu_d,  prec_d)         # log(B_j - A_j)
    z0[j]     ~ dnorm(mu_z0, prec_z0)        # logit(x0_j)
    eta_s[j]  ~ dnorm(mu_eta_s, prec_eta_s)  # log(scal_j)
    
    # Transformations to enforce constraints
    A[j]     <- exp(a[j])                    # A_j > 0
    B[j]     <- A[j] + exp(d[j])             # B_j > A_j
    alpha[j] <- 1 / (1 + exp(-z0[j]))        # x_{j0} in (0,1)
    scal[j]  <- exp(eta_s[j])                # scal_j > 0
  }
  
  #------------------------------------------------------
  # Hyperpriors (data-anchored, weakly informative)
  #------------------------------------------------------
  # mean of log(A_j) anchored at log(lower length)
  # precision = 1  -> sd = 1 on log scale  (~×2.7 up or down)
  mu_a  ~ dnorm(log(L_lower), 1)
  
  # mean of log(B_j - A_j) anchored at log(span)
  mu_d  ~ dnorm(log(L_span), 1)
  
  # midpoint mean around 0.5 on logit scale (logit(0.5) = 0)
  # precision = 4 -> sd = 0.5  (midpoints mostly within ~0.2–0.8)
  mu_z0 ~ dnorm(0, 4)
  
  # mean of log(scal) on normalized x axis
  # center at ~0.2 with sd = 0.707 (scale mostly within 0.05-0.8)
  mu_eta_s ~ dnorm(log(0.2), 2)
  
  # between-root variation (still fairly vague)
  prec_a      ~ dgamma(1, 1)
  prec_d      ~ dgamma(1, 1)
  prec_z0     ~ dgamma(1, 1)
  prec_eta_s  ~ dgamma(1, 1)
  
  #------------------------------------------------------
  # Global (fixed) effects on midpoint (x in [0,1])
  #------------------------------------------------------
  # precision = 4 -> sd = 0.5 on x-scale
  # so we expect shifts of order ≤ half a root length
  gamma1 ~ dnorm(0, 4)
  gamma2 ~ dnorm(0, 4)
  gamma3 ~ dnorm(0, 4)

  #------------------------------------------------------
  # Residual and hierarchical variances
  #------------------------------------------------------
  within.prec  <- pow(sigma, -2)
  sigma ~ dunif(0, 10)
}
