

#=================================================================#
#                       Beta Distribution                         #
#=================================================================#
generate_beta <- function(nsim, n, a, b, alpha, beta, sig) {
  # X <- alpha*rbeta(n*nsim,shape1=a, shape2=b, ncp=0)-beta

  sims <- lapply(seq_len(nsim), function(i) {
    data.frame(
      x = seq_len(n),
      y = alpha * rbeta(n, shape1 = a, shape2 = b, ncp = 0) + beta,
      iteration = i
    )
  })
  sims <- do.call(rbind, sims)

  muX_beta <- a / (a+b)
  varX_beta <- {
    alpha^2 * (a * b)
  } / {
    (a + b)^2 * (a + b + 1)
  }
  eps <- randn(n * nsim, 1) * sig * sqrt(var(sims$y))
  W <- sims$y + eps
  tau <- sig * sqrt(var(W))
  out <- list(dat = sims, mu = muX_beta, var = varX_beta, tau = tau, W = W)
  return(out)
}


# =================================================================#
#                    Normal Distribution                          #
# =================================================================#
truncate_normal <- function(n, mu, sig, xlo, xhi) {

  plo <- pnorm((xlo - mu) / sig)
  phi <- pnorm((xhi - mu) / sig)
  # generate uniform [0,1] random deviates
  r <- runif(n)
  # scale to [plo, phi]
  r <- plo + (phi - plo) * r

  # * z = array of truncated normal deviates, size(z)== N
  # Invert through standard normal
  z <- qnorm(r)
  # apply shift and scale
  z <- mu + z * sig

  tru_out <- list(dat = z, plo = plo, phi = phi)

  return(tru_out)
}

generate_normal <- function(nsim, n, mu, sig_m, xlo, xhi, sig) {

  # plo, phi will be used
  truncate_normal <- truncate_normal(n = n * nsim, mu, sig_m, xlo, xhi)
  sims <- lapply(seq_len(nsim), function(i) {
    truncate_normal <- truncate_normal(n = n, mu, sig_m, xlo, xhi)
    data.frame(
      x = seq_len(n),
      y = truncate_normal$dat,
      iteration = i
    )
  })

  sims <- do.call(rbind, sims)

  plo <- truncate_normal$plo
  phi <- truncate_normal$phi

  tempZ <- pnorm(phi) - pnorm(plo)

  muX_normal <- mu + (dnorm(xhi) - dnorm(xlo)) / tempZ * sig_m
  varX_normal <- sig_m^2 * (1 - (phi * dnorm(phi) - plo * dnorm(plo)) / tempZ - ((dnorm(phi) - dnorm(plo)) / tempZ)^2)
  # varX=parA^2*varNorm;
  eps <- randn(n * nsim, 1) * sig * sqrt(varX_normal)
  W <- sims$y + eps
  tau <- sig * sqrt(var(W))
  out <- list(dat = sims, mu = muX_normal, var = varX_normal, tau = tau, W = W)
  return(out)
}

#=================================================================#
#                    Exponential Distribution                     #
#=================================================================#
generate_exponential <- function(nsim, n, mu, tauExp, alpha, beta, sig) {
  espExp <- (mu - tauExp * exp(-tauExp / mu) - mu * exp(-tauExp / mu)) / (1 - exp(-tauExp / mu)) # E(Y)

  espX <- alpha * espExp + beta # E(X) where X=aY+b

  esp2Exp <- (-tauExp^2 * exp(-tauExp / mu) - 2 * tauExp * mu * exp(-tauExp / mu) - 2 * mu^2 * exp(-tauExp / mu) + 2 * mu^2) / (1 - exp(-tauExp / mu))
  varExp <- esp2Exp - espExp^2 # V(Y)
  varX_exp <- (alpha^2) * varExp

  # X <- beta+alpha*(-mu*log(1-runif(n*nsim, 0,1)*(1-exp(-tauExp/mu))))
  sims <- lapply(seq_len(nsim), function(i) {
    data.frame(
      x = seq_len(n),
      y = beta + alpha * (-mu * log(1 - runif(n, 0, 1) * (1 - exp(-tauExp / mu)))),
      iteration = i
    )
  })
  sims <- do.call(rbind, sims)
  eps <- randn(n * nsim, 1) * sig * sqrt(varX_exp)
  W <- sims$y + eps
  tau <- sig * sqrt(var(W))
  out <- list(dat = sims, mu = espX, var = varX_exp, tau = tau, W = W)
  return(out)
}

