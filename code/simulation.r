#' @param n0 number of data to generate the data for n0 
#' @param n1 number of data to generate the data for n1
#' @param nsim number of simulation 
#' @param dist_0 name of distribution for F0 ("Exp", "Beta", "Normal")
#' @param dist_1 name of distribution for F1 ("Exp", "Beta", "Normal")
SimusBernstein(30, 30, 2, "Beta", "Beta", M = c(1, 10))

n0 <- 300; n1 <- 300; nsim = 1; dist_0 = "Beta"; dist_1 = "Beta"; M= c(1, 10); sig = 0.1
SimusBernstein <- function(n0, n1, nsim, dist_0, dist_1, M, sig ){
    setwd("/Users/youheekil/Desktop/application/simulation2/true")
    options(warn = -1) 
    # grid for the numerical evaluation of the integral
    gridx <- seq(0, 1, by = 0.001)
    if (strcmp(dist_0,'Exp') == TRUE){
            exp_ab <- matrix(c(0.5, 4.0, 4.0, 1.0, 
                            10, 20, 20, 1.0), 
                            ncol = 4, byrow = T)
            mu0 <- exp_ab[1,1]; tau0 <- exp_ab[1, 2]; alpha0 <- exp_ab[1, 3]; beta0 <- exp_ab[1, 4]
        #    mu1 <- exp_ab[2,1]; tau1 <- exp_ab[2, 2]; alpha1 <- exp_ab[2, 3]; beta1 <- exp_ab[2, 4]                
    } else if (strcmp(dist_0,'Beta') == TRUE){
          beta_par <- matrix(c(1.0, 1.0, 2, -1, # shape1, shape2, alpha, beta
                                1.0, 2.0, 2, -1, 
                                 1, 1, 7, -7,
                                2, 2, 4, -4),
                                ncol = 4, byrow = T)
    #    shape11 <- beta_par[3, 1]
    #    shape21 <- beta_par[3, 2]
    #    alpha1 <- beta_par[3, 3]
    #    beta1 <- beta_par[3, 4] 
        shape12 <- beta_par[4, 1]
        shape22 <- beta_par[4, 2]
        alpha2 <- beta_par[4, 3]
        beta2 <- beta_par[4, 4]              
    } else if (strcmp(dist_0,'Normal') == TRUE){
          normal_ab <- matrix(c(0, 1, -1.5, 1.5 ,
                                0, 1, -7, 7 ), 
                                ncol = 4, byrow = T)
    #        mu1 <- normal_ab[1,1] ; sig1 = normal_ab[1, 2]; xlo1 = normal_ab[1,3]; xhi1 = normal_ab[1,4]
            mu0 <- normal_ab[2,1] ; sig0 = normal_ab[2, 2]; xlo0 = normal_ab[2,3]; xhi0 = normal_ab[2,4]


    } else if (strcmp(dist_0,'Gamma') == TRUE){
        gamma_par <- matrix(c(0.5, 4, 
                              0.5, 1), 
                                ncol = 2, byrow = T)
        shape0 = gamma_par[2, 1]; scale0 = gamma_par[2,2]
    }

    
    if (strcmp(dist_1,'Exp') == TRUE){
            exp_ab <- matrix(c(0.5, 4.0, 4.0, 1.0, 
                            10, 20, 20, 1.0), 
                            ncol = 4, byrow = T)
        #    mu0 <- exp_ab[1,1]; tau0 <- exp_ab[1, 2]; alpha0 <- exp_ab[1, 3]; beta0 <- exp_ab[1, 4]
            mu1 <- exp_ab[2,1]; tau1 <- exp_ab[2, 2]; alpha1 <- exp_ab[2, 3]; beta1 <- exp_ab[2, 4]                
    } else if (strcmp(dist_1,'Beta') == TRUE){
          beta_par <- matrix(c(1.0, 1.0, 2, -1, # shape1, shape2, alpha, beta
                                1.0, 2.0, 2, -1,
                                 1, 1, 7, -7, # x1
                                1, 1, 7, -7), # x0
                                ncol = 4, byrow = T)
        shape11 <- beta_par[3, 1]
        shape21 <- beta_par[3, 2]
        alpha1 <- beta_par[3, 3]
        beta1 <- beta_par[3, 4] 
        #shape12 <- beta_par[4, 1]
        #shape22 <- beta_par[4, 2]
        #alpha2 <- beta_par[4, 3]
        #beta2 <- beta_par[4, 4]              
    } else if (strcmp(dist_1,'Normal') == TRUE){
          normal_ab <- matrix(c(2, 1, -7, 7,
                                0, 1, -7, 7 ),
                                ncol = 4, byrow = T)
            mu1 <- normal_ab[1,1] ; sig1 = normal_ab[1, 2]; xlo1 = normal_ab[1,3]; xhi1 = normal_ab[1,4]
        #    mu0 <- normal_ab[2,1] ; sig0 = normal_ab[2, 2]; xlo0 = normal_ab[2,3]; xhi0 = normal_ab[2,4]


    } else if (strcmp(dist_1,'Gamma') == TRUE){
        gamma_par <- matrix(c(0.5, 4, 
                              0.5, 1, 
                              2, 2), 
                                ncol = 2, byrow = T)
        #shape1 = gamma_par[1, 1]; scale1 = gamma_par[1,2]
        shape1 = gamma_par[3, 1]; scale1 = gamma_par[3,2]

    }

    # simulation parameters ==> arguments of the function
    # n = 300; % number of subjects per dataset
    # nsim = 500; % number of replicated datasets
    # dist='Beta'; % or 'Exp', or 'Normal'
    # par = [0.7,0.5,-1,0];%[4,3,2,1,-2,2];
    
    #set of values for the ratio std(eps)/std(X)
    #sigset <- c(0.1, 0.25, 0.5, 0.75)
    m <- M[2] - M[1]
  
    #for(s in 3:4){ # for each value of NSR
    #sig <- sigset[s]
    message("algorithm for sigma = ", sig, " is starting ...")

    # objects for the output
    # save for each significant

    optimal_m0 <- list();
    optimal_m1 <- list(); xx_0_sim <- list(); yy_0_sim <- list();
    xx_1_sim <- list(); yy_1_sim <- list(); optimal_param0 <- list();
    optimal_param1 <- list(); optimal_rocmable <- list()
    optimal_auc_mable <- list(); optimal_rocemp <- list(); optimal_aucemp <- list() ;
    tau_W0 <- list(); tau_W1 <- list(); alpha_W0 <- list(); beta_W0 <- list() ; 
    alpha_W1 <- list(); beta_W1 <- list()
    fileTemp1 <- paste0("result", as.character(dist_0), n0, " & ", as.character(dist_1), n1,"_", as.character(sig),".txt")
    cat("result: sig = ", sig , "\n" ,file = fileTemp1, append = FALSE)
    
    # save file (open file)
    #fileTemp1 <- paste0("result", as.character(dist_0), "&", as.character(dist_1), as.character(sig),".txt")
    #cat("result: sig = ", sig , "\n" ,file = fileTemp1, append = FALSE)

        for (ns in seq_len(nsim)){
            message(cat(paste0("Simulation : ", ns)))
            set.seed(ns + 100 * ( ns - 1 ))
    
            # data generation for F0
            if (strcmp(dist_0,'Exp') == TRUE){
                res_exp_0 <- generate_exponential(nsim = 1, n = n0, mu = mu0, tauExp = tau0, alpha = alpha0, beta = beta0, sig = sig)
                W0 <- res_exp_0$W            
            } else if (strcmp(dist_0,'Beta') == TRUE){
                res_beta_0 <- generate_beta(nsim = 1, n = n1, a = shape12, b = shape22, alpha = alpha2, beta = beta2, sig = sig)
                W0 <- res_beta_0$W
            } else if (strcmp(dist_0,'Normal')==TRUE){
                res_normal_0 <- generate_normal(nsim = 1, n = n0, mu= mu0, sig_m = sig0, xlo = xlo0, xhi= xhi0, sig = sig)
                W0 = res_normal_0$W    
            } else if (strcmp(dist_0,'Gamma')==TRUE){
                res_gamma_0 <- generate_gamma(nsim = 1, n = n0, shape = shape0, scale = scale0, sig = sig)
                W0 = res_gamma_0$W
            }

            # data generation for F1
            if (strcmp(dist_1,'Exp') == TRUE){
                res_exp_1 <- generate_exponential(nsim = 1, n = n1, mu = mu1, tauExp = tau1, alpha = alpha1, beta = beta1, sig = sig)
                W1 <- res_exp_1$W
            } else if (strcmp(dist_1,'Beta') == TRUE){
                res_beta_1 <- generate_beta(nsim = 1, n = n0, a = shape11, b = shape21, alpha = alpha1, beta = beta1, sig = sig)
                W1 <- res_beta_1$W
            } else if (strcmp(dist_1,'Normal')==TRUE){
                res_normal_1 <- generate_normal(nsim = 1, n = n1, mu= mu1, sig_m = sig1, xlo = xlo1, xhi= xhi1, sig = sig)
                W1 = res_normal_1$W            
            } else if (strcmp(dist_0,'Gamma')==TRUE){
                res_gamma_1 <- generate_gamma(nsim = 1, n = n1, shape = shape1, scale = scale1, sig = sig)
                W1 = res_gamma_1$W
            }


            # save the generated datasets
            score <- as.vector(rbind(W0, W1))
            class <- c(rep(0, length(W0)), rep(1, length(W1)))
            tau_w0 <- tau_intialization(W0, sig)
            tau_w1 <- tau_intialization(W1, sig)
            alpha_w0 <- max(W0) - min(W0); beta_w0 <- min(W0);
            alpha_w1 <- max(W1) - min(W1); beta_w1 <- min(W1);


            res <- true_compute_ROC(score = score, class = class, M = M, sig = sig)

            optimal_m0 <- append(optimal_m0, list(res$m0))
            optimal_m1 <- append(optimal_m1, list(res$m1))
            xx_0_sim <- append(xx_0_sim, list(res$xx_0))
            yy_0_sim <- append(yy_0_sim, list(res$yy_0))
            xx_1_sim <- append(xx_1_sim, list(res$xx_1))
            yy_1_sim <- append(yy_1_sim, list(res$yy_1))
            optimal_param0 <- append(optimal_param0, list(res$param_0))
            optimal_param1 <- append(optimal_param1, list(res$param_1))
            optimal_rocmable <- append(optimal_rocmable, list(res$ROC_mable))
            optimal_rocemp <- append(optimal_rocemp, list(res$ROC_emp))
            optimal_aucemp <- append(optimal_aucemp, list(res$AUC_emp))
            optimal_auc_mable <- append(optimal_auc_mable, list(res$AUC_mable_decon))
            tau_W0 <- append(tau_W0, tau_w0)
            tau_W1 <- append(tau_W1, tau_w1)
            alpha_W0 <- append(alpha_W0, alpha_w0)
            alpha_W1 <- append(alpha_W1, alpha_w1)
            beta_W0 <- append(beta_W0, beta_w0)
            beta_W1 <- append(beta_W1, beta_w1)



        }

        tau0 <- c(); tau1 <- c(); alpha0 <- c(); alpha1 <- c(); beta0 <- c(); beta1 <- c()
        for(i in seq_len(nsim)){
            t0 <- optimal_param0[[i]][1]
            tau0 <- append(tau0, t0)
            t1 <- optimal_param1[[i]][1]
            tau1 <- append(tau1, t1)
            a0 <- optimal_param0[[i]][2]
            alpha0 <- append(alpha0, a0)
            a1 <- optimal_param1[[i]][2]
            alpha1 <- append(alpha1, a1)
            b0 <- optimal_param0[[i]][3]
            beta0 <- append(beta0, b0)
            b1 <- optimal_param1[[i]][3]
            beta1 <- append(beta1, b1)
        }
        avg <- list( original_tau0 = round(c(mean(unlist(tau_W0)), sd(unlist(tau_W0))), 3), 
                     original_tau1 = round(c(mean(unlist(tau_W1)), sd(unlist(tau_W1))), 3),
                     avg_tau0 = round(c(mean(tau0), sd(tau0)),3),
                     avg_tau1 = round(c(mean(tau1), sd(tau1)),3),
                     original_alpha0 = round(c(mean(unlist(alpha_W0)), sd(unlist(alpha_W0))), 3), 
                     original_alpha1 = round(c(mean(unlist(alpha_W1)), sd(unlist(alpha_W1))), 3),
                     avg_alpha0 = round(c(mean(alpha0), sd(alpha0)),3),
                     avg_alpha1 = round(c(mean(alpha1), sd(alpha1)),3),
                     original_beta0 = round(c(mean(unlist(beta_W0)), sd(unlist(beta_W0))), 3), 
                     original_beta1 = round(c(mean(unlist(beta_W1)), sd(unlist(beta_W1))), 3),
                     avg_beta0 = round(c(mean(beta0), sd(beta0)),3),
                     avg_beta1 = round(c(mean(beta1), sd(beta1)),3))
        trapezoidal_integration = function(x, f){
   
            n = length(x)
            # integrate using the trapezoidal rule
            integral = 0.5*sum((x[2:n] - x[1:(n-1)]) * (f[2:n] + f[1:(n-1)]))
            # print the definite integral
            return(integral)
        }
        m0 = unlist(optimal_m0); m1 = unlist(optimal_m1)
        roc_mib <- (matrix(((unlist(optimal_rocmable)) - unlist(optimal_rocemp)), nrow = 101, ncol = nsim))
        roc_mse <- (matrix(((unlist(optimal_rocmable)) - unlist(optimal_rocemp))^2, nrow = 101, ncol = nsim))
        mse = sum(apply(roc_mse, 2, sum))/nsim
        x = seq(0, 1, len  = 101)
        mise <- c();mib <- c()
        for(i in 1:nsim){
            mise <- append(mise, trapezoidal_integration(x = x, f = roc_mse[, i]))
            mib <- append(mib, trapezoidal_integration(x = x, f = roc_mib[, i]))
        }
        mise <- sum(mise)/nsim
        mib <- sum(mib)/nsim
    cat("\n ============================================================ ", 
  sep = " ", file = fileTemp1, append = TRUE)
  cat("\n m0:", m0, sep = ", ", file = fileTemp1, append = TRUE)
  cat("\n m1:", m1, sep = ", ", file = fileTemp1, append = TRUE)
  cat("\n tau0:", tau0 , sep = ", ", file = fileTemp1, append = TRUE)
  cat("\n tau1:", tau1 , sep = ", ", file = fileTemp1, append = TRUE)
  cat("\n alpha0:", alpha0 , sep = ", ", file = fileTemp1, append = TRUE)
  cat("\n alpha1:", alpha1 , sep = ", ", file = fileTemp1, append = TRUE)
  cat("\n beta0:", beta0 , sep = ", ", file = fileTemp1, append = TRUE)
  cat("\n beta1:", beta1 , sep = ", ", file = fileTemp1, append = TRUE)
  cat("\n avg:", unlist(avg), sep = ", ", file = fileTemp1, append = TRUE)
  cat("\n roc_emp:", unlist(optimal_rocemp) , sep = ", ", file = fileTemp1, append = TRUE)
  cat("\n roc_mable:", unlist(optimal_rocmable) , sep = ", ", file = fileTemp1, append = TRUE)
  cat("\n auc_emp:", round(sum(unlist(optimal_aucemp))/length(unlist(optimal_aucemp)), 3), sep = ", ", file = fileTemp1, append = TRUE)
  cat("\n auc_mable:", round(sum(unlist(optimal_auc_mable))/length(unlist(optimal_aucemp)), 3), sep = ", ", file = fileTemp1, append = TRUE)
  cat("\n mib:", mib, sep = ", ", file = fileTemp1, append = TRUE)
  cat("\n mise:", mise, sep = ", ", file = fileTemp1, append = TRUE)

  cat("\n ============================================================ ",
    sep = " ", file = fileTemp1, append = TRUE)



        col = c("#000000", "#d7d7d8") 
        namef0 <- paste0("f0_", as.character(dist_0), n0, " & ", as.character(dist_1),n1, "_sig_", sig)
        png(paste0(namef0, ".png"))
        par(mfrow = c(1, 1))
        graphics::plot(density(W0),
            col = col[1], lwd = 2, lty = 1)
        
        for(i in 1:nsim){
            graphics::lines(xx_0_sim[[i]], yy_0_sim[[i]], type = "l",
                col = col[2], lwd = 2, lty = 1)
        }

        graphics::lines(density(W0), 
        col = col[1], lwd = 2, lty = 1)
        dev.off()

        namef1 <- paste0("f1_", as.character(dist_0), n0, " & ", as.character(dist_1),n1, "_sig_", sig)
        png(paste0(namef1, ".png"))
        graphics::plot(density(W1),  col = col[1],
        lwd = 2, lty = 1)
        for(i in 1:nsim){
            graphics::lines(xx_1_sim[[i]], yy_1_sim[[i]], type = "l",
                col = col[2], lwd = 2, lty = 1)
        }
        graphics::lines(density(W1),  col = col[1],
        lwd = 2, lty = 1)
        dev.off()



}
