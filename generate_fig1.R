#visualization for weak, moderate and srong overlap in Web Figure 1
#Author: Zhiqiang Cao
###Date: 2022.12.19

rm(list = ls())

library(survival) 
library(mvtnorm)
library(data.table)

#etermine true value of a0
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
    #print(i)
  }
  return(A0[which.min(abs(PREV-prev))])
}


a0_cand <- sapply(list(1,2,3,4,5), function(g) generatea0(prev=0.5, lower=-5, upper=5, gamma=g))

###gampar=1,3,5 corresponding to strong,moderate and weak overlap
visualizedata <- function(gampar,choice){
  a0 <- a0_cand[gampar]
  N = 1000000  #sample size
  set.seed(123456+500)
  rho = 0.5
  R <- (1-rho) * diag(3) + rho * matrix(1,3,3)
  #generate continuous covariates X4,X5,X6
  xc <- rmvnorm(N, rep(0,3), R)
  #generate binary covariates X4,X5,X6
  xb = cbind(rbinom(n=N,size=1,p=0.5),rbinom(n=N,size=1,p=0.5),rbinom(n=N,size=1,p=0.5))
  x = cbind(xc,xb)
  
  # selection probability for trial
  ac <- gampar*c(0.15, 0.3, 0.3)
  ab <- -gampar*c(0.2, 0.25, 0.25)
  a <- c(a0, ac, ab)
  e <- plogis(c(cbind(1,x)%*%a))
  z <- rbinom(N, 1, e)
  # histograms
  if(choice==1){
    hist(e[z==1],breaks=35,col="gray77",border="gray77",main=paste("(a) strong overlap"),
         xlim=c(0,1),ylim=c(0,3),xlab="Propensity score",freq=F,
         cex.lab = 1, cex.axis = 1 ,cex.main = 1, cex = 1)
    hist(e[z==0],breaks=35,col="NA",add=TRUE,freq=F)
  }else if(choice==2){
    hist(e[z==1],breaks=35,col="gray77",border="gray77",main=paste("(b) moderate overlap"),
         xlim=c(0,1),ylim=c(0,3),xlab="Propensity score",freq=F,
         cex.lab = 1, cex.axis = 1 ,cex.main = 1, cex = 1)
    hist(e[z==0],breaks=35,col="NA",add=TRUE,freq=F)
  }else if(choice==3){
    hist(e[z==1],breaks=35,col="gray77",border="gray77",main=paste("(c) weak overlap"),
         xlim=c(0,1),ylim=c(0,6),xlab="Propensity score",freq=F,
         cex.lab = 1, cex.axis = 1 ,cex.main = 1, cex = 1)
    hist(e[z==0],breaks=35,col="NA",add=TRUE,freq=F)
  }
}

op <- par(mfrow = c(1, 3), pty = "s")
visualizedata(gampar=1,choice=1)
visualizedata(gampar=3,choice=2)
visualizedata(gampar=4,choice=3)
par(op)
