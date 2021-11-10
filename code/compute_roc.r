#' @title compute_ROC
#' @import check 'requirement.r' 
#' @description Function \code{roc.mable} 
#' 1 - F(X = x) where \eqn{F} is the cumulative density 
#' function (CDF) of \eqn{X}
#' @param score A vector containing (diagnostic) scores assuming contagious data
#' @param class A vector containg the class
#' @title estimation_param 
#' @function area, getcdf
#' @param W matrix of size 'n0' and 'n1' for F0 and F1, respectively
#' @param m (nbTerms) maximum number of terms to be considered 
#' in the maxinum likelihood Bernstein polynomial 
tau_intialization <- function(w, sig){
    s2 <- var(w)
    tau <- sig * as.vector(sqrt(s2))
    return(tau)
}

cdf_f <- function(x, y) {
    rval <- approxfun(x = x, y = y,
        method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered"
    )
    class(rval) <- c("ecdf", "stepfun", class(rval))
    assign("nobs", n, envir = environment(rval))
    attr(rval, "call") <- sys.call()
    rval
}
area <- function(x,y){

  mydata <- as.data.frame(cbind(x,y))
  sorted <- mydata[order(x),]
  x <- sorted$x
  y <- sorted$y
  return(sum(diff(x)*(y[-1]+y[-length(y)]))/2)
}

getcdf <- function(x){
    if(!(methods::is(x, "density"))){
    stop("x is not of class \"density\" ")
  }
  density <- x
  densityxvals <- density$x
  area(densityxvals, density$y)
}

survival <- function(x, y, cutoff){
    densityxvals <- xx_0 ; densityyvals <- yy_0
    indexes <- which(densityxvals >= c_vals)
    ROCit::trapezoidarea(densityxvals[indexes], densityyvals[indexes])
}

simple_auc <- function(TPR, FPR){
    # inputs already sorted, best scores first 
        dFPR <- c(diff(FPR), 0)
        dTPR <- c(diff(TPR), 0)
        sum(TPR * dFPR) + sum(dTPR * dFPR)/2
}


logLikabBernsteinMat2 <- function(param, gridx, m, data){

  gridxnew <- sort(param[2]*gridx+param[3])
  #fw <- densityWabBernsteinMat(param = param, gridx = gridxnew, m = m, data = data)
  tempF <- -sum(log(densityWabBernsteinMat(param = param, gridx = gridxnew, m = m, data = data)))
  if (is.nan(tempF)){
    llk <- 999999
  } else {
    llk <- tempF
  }
  return(llk)
}

densitySigMat <- function(tau, gridx) {
  f <- dnorm(gridx, 0, tau) # 0 = mu, sig = sigma
  return(f)
}


densityWabBernsteinMat <- function(param, gridx, m, data){
    n <- length(data)
    u <- matrix(data, nrow = length(data), ncol = length(gridx)) - repmat(gridx,n,1)
    fu <- densitySigMat(param[1], u)
    #new_grid <- (gridx - param[3])/param[2]
    fb <- densityXBernstein(gridx, m = m, param = param, data = data) # m x gridx
    #fb <- repmat(densityXBernstein(new_grid, m = m, param = param, data = data), n,1) #nxgridx * m 
    f <- sum(fb%*%t(fu))/(param[2]*param[1])
    f <- f[f > 0]
  return(f)
}

densityXBernstein <- function(gridx, m, param, data) { 
  k <- seq(0, m, 1)
  theta <- param[4:length(param)]
  shape1 <- matrix(k+1, nrow = m+1, ncol = length(gridx))
  shape2 <- matrix(m-k+1, nrow = m+1, ncol = length(gridx))
  f <-  theta * dbeta(repmat(gridx, m + 1, 1), shape1, shape2)
  return(f)
}

# measurement error model 

controls <- list(sig.level=1.0e-2, eps = 1.0e-7, maxit = 5000, eps.em = 1.0e-7, maxit.em = 5000,
    eps.nt = 1.0e-7, maxit.nt = 1000, tini=1e-4)

eta_psi <- function(w, m, tau) {
    a <- min(w)-tau; b <- max(w)+tau; beta <- a ; alpha <- b - a
    y <- (w-a)/(b-a) # we can only estimate f as a density with support [x_1, x_n]. IF the support of f  differs from [0, 1] and we can find a finite interval 
    psi <- c()
    n <- length(y)
    psi <- c()
  for(i in seq_len(n)){
    for (p in seq_len(m + 1)) {
      k <- p - 1
      gn <- function(x) dnorm(x, 0, tau) # 0.2 will be replaced to tau 
      ff <- function(x) (b-a)*gn((b-a)*x)
      eta_mj <- function(x) {
        as.vector((m+1)*dbinom(k, m, x, 1-x)*ff(y[i]-x))
      }
      result <- adaptIntegrate(eta_mj, lowerLimit = 0, upperLimit = 1)$integral
      psi[i + k * n] <- result
    }
  }
  return(psi)
}


# EM algorithm for error* beta(i+1, m+1-i), i = 0, ..., m , n = sample size 
em_gBeta_mix <- function(w, p, m, tau, controls){
    n <- length(w)
    gbeta <- rep(0, (m+1)*n)
    p_gbeta <- rep(0, (m+1)*n)
    fp <- rep(0, n)
    a <- min(w)-tau; b <- max(w)+tau; beta <- a ; alpha <- b - a ; llik <- c()
    y <- (w-a)/(b-a) 
    psi_m <- c(); p_gbeta <- c()
    gbeta <- eta_psi(w, m, tau)
    llik[1] <- 0 
    for(i in seq_len(n)){
        fp[i] <- 0 
        for(j in seq_len(m+1)){
            p_gbeta[i+n*(j-1)] <- p[j] * gbeta[i+n*(j-1)]
            fp[i] <- fp[i] + p_gbeta[i + n * (j-1)]
        }
        llik <- llik + log(fp[i])
    }

    # iteration of EM algorithm
 
    ifelse(m > 0, del <- 10, del <- 0) # m must be positive
    it <- 1
    controls <- controls
    eps <- controls$eps
    maxit <- controls$maxit

    while(del>eps && it<maxit){
        for(j in seq_len(m+1)){
            p[j] <- 0 
            for(i in seq_len(n)){ 
                p[j] <- p[j] + p_gbeta[i + n * (j-1)] / fp[i]
            }
            p[j] <- p[j] / n
        }
        
        llik_nu <- 0
        for(i in i:n){
            fp[i] <- 0 
            for(j in seq_len(m+1)){
                p_gbeta[i+n*(j-1)] = p[j] * gbeta[i+n*(j-1)]
                fp[i] = fp[i] + p_gbeta[i+n*(j-1)]
            }
            llik_nu = llik_nu + log(fp[i])
        }
        del = abs(llik-llik_nu)
        it = it + 1
        llik = llik_nu
        message(paste0("Iteration: ", it-1, "\nDel = ",  del))
    }

    out <- list(gBeta = gbeta, 
                p_gBeta = p_gbeta, 
                fp = fp, 
                llik = llik, 
                p = p)
    return(out)
}

compute_ROC <- function(score, class, M = M, sig = sig){
  # obtain the marker in healthy and diseased group 
  options(warn = -1)

  gridx <- seq(0, 1, 0.001)
  xx_1 <- c(); xx_0 <- c()
  yy_1 <- c(); yy_1 <- c()

  tempdata <- ROCit::rankorderdata(score, class)
  D <- tempdata[, 1]
  Y <- tempdata[, 2]
  Ybar <- 1 - Y
  DY <- D[which(Y == 1)] # diseaded group 
  DYbar <- D[which(Y == 0)] # non-diseased group 

  # omit missing data
      # Missing data
  omit.1 <- is.na(DY) # diseased group 
  omit.0 <- is.na(DYbar) # non-diseased group 

  DY_n <- DY[!omit.1]
  DYbar_n <- DYbar[!omit.0]

  n1 <- length(DY_n) 
  n0 <- length(DYbar_n)

  nY <- sum(Y)
  nYbar <- sum(Ybar)

  TP <- cumsum(Y)
  FP <- cumsum(Ybar)
  TPR_emp <- TP / nY
  FPR_emp <- FP / nYbar

  # MABLE 
  m <- M[2] - M[1]

  est_par <- param_est(DY_n = DY_n, DYbar_n = DYbar_n)

  tau_DY <- tau_intialization(DY_n, est_par$sig1)
  #tau_DY <- tau_intialization(DY_n, sig)

  tau_DYbar <- tau_intialization(DYbar_n, est_par$sig0)
  #tau_DYbar <- tau_intialization(DYbar_n, sig)

  a1 <- min(DY_n) - tau_DY; b1 <- max(DY_n)+tau_DY; beta1 <- a1 ; alpha1 <- b1- a1
  a0 <- min(DYbar_n) - tau_DYbar; b0 <- max(DYbar_n)+tau_DYbar; beta0 <- a0 ; alpha0 <- b0- a0

  res_1 <- mable.decon(y = DY_n, gn = function(x) dnorm(x, mean = 0, sd = tau_DY), M = c(1, 10), interval = c(a1, b1), progress = FALSE)
  res_0 <- mable.decon(y = DYbar_n, gn = function(x) dnorm(x, mean = 0, sd = tau_DYbar), M = c(1, 10), interval  = c(a0, b0), progress = FALSE)
  param1 <- c(tau_DY, alpha1, beta1, res_1$p)
  param0 <- c(tau_DYbar, alpha0, beta0, res_0$p)

  # diseased => density
  xx_1 <- seq(a1, b1, len = 512)
  yy_1 <- dmixbeta(xx_1, res_1$p, c(a1, b1))
  yy_1c <- pmixbeta(xx_1, res_1$p, c(a1, b1))
  # not diseased => density
  xx_0 <- seq(a0, b0, len = 512)
  yy_0 <- dmixbeta(xx_0, res_0$p, c(a0, b0))
  yy_0c <- pmixbeta(xx_0, res_0$p, c(a0, b0))
  # ROC calculate
  F1emp <- stats::ecdf(DY_n) # diseased
  p <- seq(0.01, 0.99, length = 512)
  new_p <- seq(0.01, 0.99, length = 512)
  rocemp <- 1 - F1emp(quantile(DYbar_n, 1 - p, type = 1))

  F1mable <- cdf_f(xx_1, yy_1c) # diseased
  rocmable <- 1 - F1mable(quantile(xx_0, 1-yy_0c, type = 1))
  aucemp <- sum( outer(DYbar_n, DY_n, "<") ) / (n0 * n1) + sum(outer(DYbar_n, DY_n, "==")) / (2*n0*n1)
  #emp_auc <- round(sum(aucemp), 3)
  emp_auc <- round(area(new_p, rocemp), 3)
  mable_auc <- round(area(new_p, rocmable), 3)

  out <- list(diseased = DY_n, non_diseased = DYbar_n,
              n0 = n0, n1 = n1, TPR = TPR_emp, FPR = FPR_emp, 
              m0 = res_0$m, m1 = res_1$m,
              param_0 = param0, param_1 = param1,
              #bic0 = res_0$best_bic, bic1 = res_1$best_bic,
              xx_0 = xx_0, xx_1 = xx_1, yy_0 = yy_0, yy_1 = yy_1,
              ROC_mable = rocmable, ROC_emp = rocemp,
              AUC_mable_decon = mable_auc, AUC_emp = emp_auc
            )
  return(out)
}

