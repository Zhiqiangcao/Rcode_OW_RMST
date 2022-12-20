###calculate true values of IPTW, OW, IPTW with symmetric trimming and 
###asymmetric trimming weights under each case for simulations in paper
###i.e., reproduce simulation results in Table 1 of paper

rm(list=ls())
library(survival) 
library(mvtnorm)
library(data.table)

#determine true value of a0
generatea0 <- function(prev, lower, upper, gamma, nsim=100){
  A0 <- seq(lower, upper, length = nsim)
  PREV <- rep(0,nsim)
  N <- 100000
  
  # covariates
  rho <- 0.5
  #########comment: to keep consistent with main function of generating covariates
  R <- (1-rho) * diag(3) + rho * matrix(1,3,3)
  xc <- rmvnorm(N, rep(0,3), R)
  # generate binary covariates X4 X5 X6
  xb <- cbind(rbinom(n=N,size=1,p=0.5),rbinom(n=N,size=1,p=0.5),rbinom(n=N,size=1,p=0.5))
  W <- cbind(1, xc, xb)
  
  # coefficients
  ac <- gamma*c(0.15, 0.3, 0.3)
  ab <- -gamma*c(0.2, 0.25, 0.25)
  
  # Monte Carlo estimates of prev
  for(i in 1:nsim){
    a0 <- A0[i]
    a <- c(a0, ac, ab)
    e <- plogis(c(W%*%a))
    PREV[i] <- mean(e)
    print(i)
  }
  return(A0[which.min(abs(PREV-prev))])
}


a0_cand <- sapply(list(1,2,3,4,5), function(g) generatea0(prev=0.5, lower=-5, 
                                                          upper=5, gamma=g))
#0.3535354 0.6565657 1.0606061 1.3636364 1.7676768

######################################
nsim = 1000
N = 1000000  #sample size
#gamma <- 5 #1,2,3,4,5
#L = 2 #2, 5, 10

L_set = c(2,5,10)
gamma_set = c(1,3,5)
for(L in L_set){
  for(gamma in gamma_set){
    a0 <- a0_cand[gamma]
    #8 different estimators
    true_eff_iptw <- matrix(0, nrow=nsim, ncol=3)
    true_eff_ow <- matrix(0, nrow=nsim, ncol=3)
    true_eff_strim1 <- true_eff_strim2 <- true_eff_strim3 <- matrix(0, nrow=nsim, ncol=3)
    true_eff_atrim1 <- true_eff_atrim2 <- true_eff_atrim3 <- matrix(0, nrow=nsim, ncol=3)
    
    for(i in 1:nsim){
      cat("iter=",i,"\n")
      set.seed(123456+i)
      # Generate covariates
      rho <- 0.5
      R <- (1-rho) * diag(3) + rho * matrix(1,3,3)
      # generate continuous covariates x1-x3
      xc <- rmvnorm(N, rep(0,3), R)
      # generate binary covariates X4-X6
      xb <- cbind(rbinom(n=N,size=1,p=0.5),rbinom(n=N,size=1,p=0.5),rbinom(n=N,size=1,p=0.5))
      xpop <- cbind(xc,xb)
      
      #for treatment and control
      samp <- function(par,x){
        exp(x %*% par)/(1+exp(x %*% par))
      }
      # selection probability for trial
      ac <- gamma*c(0.15, 0.3, 0.3)
      ab <- -gamma*c(0.2, 0.25, 0.25)
      a <- c(a0, ac, ab)
      
      
      lambda_trt <- exp(cbind(1,xpop) %*% c(-1,0.4,0.2,0.1,-0.1,-0.2,-0.3))
      lambda_crl <- exp(cbind(1,xpop) %*% c(-1.4,0.0,-0.2,-0.3,-0.5,-0.6,-0.7))
      #for survial time
      upop_trt <- runif(N); upop_crl <- runif(N)
      tpop_trt = -log(upop_trt)/(lambda_trt) #T^1 for treatment
      tpop_crl = -log(upop_crl)/(lambda_crl) #T^0 for control
      #selection Non-participation in larger study 
      gampar <- c(-6.5,-0.2,-0.1,-0.5, 0.1, 0.2,-0.1)
      p_selS <- samp(gampar,cbind(1,xpop))
      S <- rbinom(N,1,p_selS)	
      tpop_trt_s0 <- tpop_trt[S==0]
      tpop_crl_s0 <- tpop_crl[S==0]
      xpop_S0 <- xpop[S==0,]
      N_S0 <- sum(S==0)
      e <- plogis(c(cbind(1,xpop_S0)%*%a))
      z <- rbinom(N_S0, 1, e)
      
      ##iptw
      iptw1 <- mean(pmin(tpop_trt_s0,L))
      iptw0 <- mean(pmin(tpop_crl_s0,L))
      
      ##ow
      h1 <- e*(1-e)
      ow1 <- mean(h1*pmin(tpop_trt_s0,L))/mean(h1)
      ow0 <- mean(h1*pmin(tpop_crl_s0,L))/mean(h1)
      
      ## Symmetric trimming
      q <- 0.05
      h2 <- ((e >= q) & (e <= (1-q)))
      strim1_1 <- mean(h2*pmin(tpop_trt_s0,L))/mean(h2)
      strim0_1 <- mean(h2*pmin(tpop_crl_s0,L))/mean(h2)
      
      q <- 0.1
      h2 <- ((e >= q) & (e <= (1-q)))
      strim1_2 <- mean(h2*pmin(tpop_trt_s0,L))/mean(h2)
      strim0_2 <- mean(h2*pmin(tpop_crl_s0,L))/mean(h2)
      
      q <- 0.15
      h2 <- ((e >= q) & (e <= (1-q)))
      strim1_3 <- mean(h2*pmin(tpop_trt_s0,L))/mean(h2)
      strim0_3 <- mean(h2*pmin(tpop_crl_s0,L))/mean(h2)
      
      ## Asymmetric trimming
      ps0 <- e[z==0]
      ps1 <- e[z==1]
      lps <- max(min(ps0), min(ps1))
      ups <- min(max(ps0), max(ps1))
      
      p <- 0
      keep <- rep(NA, N_S0)
      alpha0 <- as.numeric(quantile(ps0, 1-p))
      alpha1 <- as.numeric(quantile(ps1, p))
      keep[z==0] <- ((ps0 >= alpha1) & (ps0 <= alpha0) & (ps0 >= lps) & (ps0 <= ups))
      keep[z==1] <- ((ps1 >= alpha1) & (ps1 <= alpha0) & (ps1 >= lps) & (ps1 <= ups))
      h3 <- keep
      atrim1_1 <- mean(h3*pmin(tpop_trt_s0,L))/mean(h3)
      atrim0_1 <- mean(h3*pmin(tpop_crl_s0,L))/mean(h3)
      
      p <- 0.01
      keep <- rep(NA, N_S0)
      alpha0 <- as.numeric(quantile(ps0, 1-p))
      alpha1 <- as.numeric(quantile(ps1, p))
      keep[z==0] <- ((ps0 >= alpha1) & (ps0 <= alpha0) & (ps0 >= lps) & (ps0 <= ups))
      keep[z==1] <- ((ps1 >= alpha1) & (ps1 <= alpha0) & (ps1 >= lps) & (ps1 <= ups))
      h3 <- keep
      atrim1_2 <- mean(h3*pmin(tpop_trt_s0,L))/mean(h3)
      atrim0_2 <- mean(h3*pmin(tpop_crl_s0,L))/mean(h3)
      
      p <- 0.05
      keep <- rep(NA, N_S0)
      alpha0 <- as.numeric(quantile(ps0, 1-p))
      alpha1 <- as.numeric(quantile(ps1, p))
      keep[z==0] <- ((ps0 >= alpha1) & (ps0 <= alpha0) & (ps0 >= lps) & (ps0 <= ups))
      keep[z==1] <- ((ps1 >= alpha1) & (ps1 <= alpha0) & (ps1 >= lps) & (ps1 <= ups))
      h3 <- keep
      atrim1_3 <- mean(h3*pmin(tpop_trt_s0,L))/mean(h3)
      atrim0_3 <- mean(h3*pmin(tpop_crl_s0,L))/mean(h3)
      
      ###iptw
      true_eff_iptw[i,1] = iptw1
      true_eff_iptw[i,2] = iptw0
      true_eff_iptw[i,3] = true_eff_iptw[i,1]-true_eff_iptw[i,2]
      ###ow
      true_eff_ow[i,1] = ow1
      true_eff_ow[i,2] = ow0
      true_eff_ow[i,3] = true_eff_ow[i,1]-true_eff_ow[i,2]
      ### Symmetric trimming
      true_eff_strim1[i, 1] = strim1_1
      true_eff_strim1[i, 2] = strim0_1
      true_eff_strim1[i, 3] = true_eff_strim1[i,1]-true_eff_strim1[i,2]
      
      true_eff_strim2[i, 1] = strim1_2
      true_eff_strim2[i, 2] = strim0_2
      true_eff_strim2[i, 3] = true_eff_strim2[i,1]-true_eff_strim2[i,2]
      
      true_eff_strim3[i, 1] = strim1_3
      true_eff_strim3[i, 2] = strim0_3
      true_eff_strim3[i, 3] = true_eff_strim3[i,1]-true_eff_strim3[i,2]
      
      ## Asymmetric trimming
      true_eff_atrim1[i, 1] = atrim1_1
      true_eff_atrim1[i, 2] = atrim0_1
      true_eff_atrim1[i, 3] = true_eff_atrim1[i,1]-true_eff_atrim1[i,2]
      
      true_eff_atrim2[i, 1] = atrim1_2
      true_eff_atrim2[i, 2] = atrim0_2
      true_eff_atrim2[i, 3] = true_eff_atrim2[i,1]-true_eff_atrim2[i,2]
      
      true_eff_atrim3[i, 1] = atrim1_3
      true_eff_atrim3[i, 2] = atrim0_3
      true_eff_atrim3[i, 3] = true_eff_atrim3[i,1]-true_eff_atrim3[i,2]
      
    }
    colMeans(true_eff_iptw)
    colMeans(true_eff_ow)
    colMeans(true_eff_strim1)
    colMeans(true_eff_strim2)
    colMeans(true_eff_strim3)
    colMeans(true_eff_atrim1)
    colMeans(true_eff_atrim2)
    colMeans(true_eff_atrim3)
    
    write.csv(data.table('true_eff_iptw' = colMeans(true_eff_iptw), 
                         'true_eff_ow' = colMeans(true_eff_ow),
                         'true_eff_strim1' = colMeans(true_eff_strim1), 
                         'true_eff_strim2' = colMeans(true_eff_strim2), 
                         'true_eff_strim3' = colMeans(true_eff_strim3), 
                         'true_eff_atrim1' = colMeans(true_eff_atrim1),
                         'true_eff_atrim2' = colMeans(true_eff_atrim2),
                         'true_eff_atrim3' = colMeans(true_eff_atrim3)),
              paste0('true_val_new_L', L, "_G", gamma, ".csv"))
  }
  
}


