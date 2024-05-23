###select PS, select CS, uisng 100 covariates to fit both PS and CS (94 covariates are nuisance covariates) 
###first, conduct lasso for both PS and CS models, then refit PS and CS models for RMST

rm(list=ls())
library(survival) 
library(mvtnorm)
library(glmnet) 

#input source code
source("C:/Users/82655/Dropbox/research/ow/R_code/new_simu/github_code/RMST_functions_v3.R")
setwd("C:/Users/82655/Dropbox/research/ow/estimation_results")

#determine true value of a0
generatea0 <- function(prev, lower, upper, gamma, nsim=100){
  A0 <- seq(lower, upper, length = nsim)
  PREV <- rep(0,nsim)
  N <- 100000
  # covariates
  rho <- 0.5
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
###main function
nsim <- 1000  #simulation times
N <- 500     #sample size
SurvTime = "time"
Status = "delta"
Treatment = "z"


###Main simulation
gamma_set = 1 #c(1,3,5)
L_set = c(2,5,10)
for(gamma in gamma_set){
  for(L in L_set){
    ###true value
    if (gamma == 1){
      if (L == 2) {
        true_ow_val<- c(1.503956526,1.776999254,-0.273042727)
        true_iptw_val <- c(1.498194346,1.776465598,-0.278271251)
      } else if (L == 5) {
        true_ow_val <- c(2.693529714,3.805216113,-1.111686399)
        true_iptw_val <- c(2.68834376,3.803739793,-1.115396033)
      } else if (L == 10) {
        true_ow_val <- c(3.564410608,6.129209984,-2.564799375)
        true_iptw_val <- c(3.580969651,6.127871625,-2.546901973)
      } 
    } else if (gamma == 3) {
      if (L == 2) {
        true_ow_val <- c(1.519824121,1.778325209,-0.258501088)
        true_iptw_val <- c(1.498194346,1.776465598,-0.278271251)
      } else if (L == 5) {
        true_ow_val <- c(2.708913774,3.808817305,-1.099903531)
        true_iptw_val <- c(2.68834376,3.803739793,-1.115396033)
      } else if (L == 10) {
        true_ow_val <- c(3.516710564,6.132360347,-2.615649783)
        true_iptw_val <- c(3.580969651,6.127871625,-2.546901973)
      } 
    } else if (gamma == 5) {
      if (L == 2) {
        true_ow_val <- c(1.52741699,1.778899222,-0.251482232)
        true_iptw_val <- c(1.498194346,1.776465598,-0.278271251)
      } else if (L == 5) {
        true_ow_val <- c(2.716895706,3.81034556,-1.093449854)
        true_iptw_val <-  c(2.68834376,3.803739793,-1.115396033)
      } else if (L == 10) {
        true_ow_val <- c(3.492736569,6.133647604,-2.640911036)
        true_iptw_val <- c(3.580969651,6.127871625,-2.546901973)
      } 
    }
    
    #keep estimation results for each replicates
    iptw_est = ow_est = matrix(0, ncol=3, nrow=nsim)
    iptw_se = ow_se = matrix(0, ncol=3, nrow=nsim)
    iptw_cp = ow_cp = matrix(0, ncol=3, nrow=nsim)
    
    a0 <- a0_cand[gamma]
    
    for(i in 1:nsim){
      cat(paste0('L', L, "_G", gamma,"_iter=",i),"\n")
      set.seed(123456+i)
      # Generate covariates
      rho = 0.5
      R <- (1-rho) * diag(3) + rho * matrix(1,3,3)
      #generate continuous covariates x1-x3
      xc <- rmvnorm(N, rep(0,3), R)
      #generate binary covariates X4-X6
      xb = cbind(rbinom(n=N,size=1,p=0.5),rbinom(n=N,size=1,p=0.5),rbinom(n=N,size=1,p=0.5))
      x = cbind(xc,xb)
      # selection probability for trial
      ac <- gamma*c(0.15, 0.3, 0.3)
      ab <- -gamma*c(0.2, 0.25, 0.25)
      a <- c(a0, ac, ab)
      e <- plogis(c(cbind(1,x)%*%a))
      z <- rbinom(N, 1, e)
      lambda_trt <- exp(cbind(1,x) %*% c(-1,0.4,0.2,0.1,-0.1,-0.2,-0.3))
      lambda_crl <- exp(cbind(1,x) %*% c(-1.4,0.0,-0.2,-0.3,-0.5,-0.6,-0.7))
      #for survial time
      u_trt <- runif(N); u_crl <- runif(N)
      t_trt <- -log(u_trt)/(lambda_trt) #T^1 for treatment
      t_crl <- -log(u_crl)/(lambda_crl) #T^0 for control
      #for censoring time
      lambdac_trt <- exp(cbind(1,x) %*% c(-1.6,-0.3,0.5,0.5,0.2,-0.4,-0.5))
      lambdac_crl <- exp(cbind(1,x) %*% c(-1.6,-0.3,0.5,0.5,0.2,-0.4,-0.5))
      uc_trt <- runif(N); uc_crl <- runif(N)
      c_trt <- -log(uc_trt)/(lambdac_trt) #C^1 for treatment
      c_crl <- -log(uc_crl)/(lambdac_crl) #C^0 for control
      #T = z*T^1+(1-z)*T^0
      t = z*t_trt + (1-z)*t_crl; 
      #C = z*C^1+(1-z)*C^0
      c = z*c_trt + (1-z)*c_crl; 
      #censoring indicator
      delta = 1*(t<=c)
      #observed time
      time = delta*t + (1-delta)*c
      #generate extra nuisance continuous covariates and binary covariates
      Rn <- (1-rho) * diag(47) + rho * matrix(1,47,47)
      xc_n <- rmvnorm(N, rep(0,47), Rn)
      xb_n = NULL
      for(j in 1:47){
        gen_binary = rbinom(n=N,size=1,p=0.5)
        xb_n = cbind(xb_n,gen_binary)
      }
      x = cbind(xc,xc_n,xb,xb_n)
      Covariates = paste0("x",1:100,sep = "")
      colnames(x) = Covariates
      Data = data.frame(time=time,delta=delta,z=z,x)
      ps.formula0 = as.formula(paste(Treatment,"~",paste(Covariates,collapse="+"),sep=""))
      sort_index = order(Data[,SurvTime])
      data.sort = Data[sort_index,]
      Data = data.sort
      #next, we conduct Lasso for PS
      z = as.double(Data[,Treatment])
      x = as.matrix(Data[,-(1:3)])
      cv.ps.lasso = cv.glmnet(x, z, family='binomial', alpha=1, standardize=TRUE, type.measure='auc')
      #best_ps_lambda = cv.ps.lasso$lambda.min
      #best_ps_model = glmnet(x, z, alpha = 1, family='binomial', lambda = best_ps_lambda)
      ps_coef = coef(cv.ps.lasso, s=cv.ps.lasso$lambda.min) #coef(best_ps_model)
      ps_index = which(ps_coef!=0)
      xall = model.matrix(ps.formula0,data=Data)
      if(ps_index[1]==1){  #already include intercept
        w = xall[,ps_index]  
      }else{
        w = cbind(1,xall[,ps_index]) #include intercept
      }
      colnames(w) = paste("w",1:dim(w)[2],sep="")
      #now, we conduct lasso for Cox
      y = as.matrix(Data[,c(1,2)])
      colnames(y) = c("time","status")
      cv.cox.lasso <- cv.glmnet(x, y, family = "cox", alpha = 1, type.measure = "C")
      #best_cox_lambda = cv.cox.lasso$lambda.min
      #best_cox_model <- glmnet(x, y, family = "cox",  alpha = 1, lambda = best_cox_lambda)
      #cox_coef = coef(best_cox_model)
      cox_coef = coef(cv.cox.lasso, s=cv.cox.lasso$lambda.min)
      cox_index = which(cox_coef!=0)
      xori = xall[,-1]
      len_cox_index = length(cox_index)
      if(len_cox_index == 0){ #then no covariates in Cox model
        x = matrix(0,nrow=N)  #then estimate of x will be NA in censored
      }else if(len_cox_index==1){
        x = matrix(xori[,cox_index],nrow=N)
      }else{
        x = matrix(xori[,cox_index],nrow=N,ncol=len_cox_index)
      }
      px = dim(x)[2]
      colnames(x)=paste("x",1:dim(x)[2],sep="")
      z = as.numeric(Data[,Treatment])
      time = as.numeric(Data[,SurvTime])
      delta = as.numeric(Data[,Status])
      # reconstruct the dataset
      data = as.data.frame(cbind(x,w,z,time,delta))
      #seperate two groups
      data.trt = subset(data,z==1)
      data.con = subset(data,z==0)
      # propensity score: update ps.formula
      ps.formula = as.formula(paste("z~",paste(colnames(w),collapse = "+"),"-1",sep="")) #since w[,1] is the intercept term
      ps.model = glm(ps.formula,data=data,family=binomial(link="logit"))  #logistic model default including intercept term
      ps = 1/(1+exp(-c(w %*% ps.model$coefficients))) 
      # update censor.formula
      censor.formula = as.formula(paste("Surv(time, I(1-delta)) ~",paste(colnames(x),collapse = "+"),sep=""))
      cen.trt.model = coxph(censor.formula, data=data.trt)
      cen.con.model = coxph(censor.formula, data=data.con)
      ###OW method
      Delta.OW = OW.RMST(L=L,w,x,z,time,delta,px,ps.model,cen.trt.model,cen.con.model)
      ow_est[i,] = Delta.OW$point
      ow_se[i,] = sqrt(Delta.OW$variance)
      oW_lower = ow_est[i,]-qnorm(0.975)*ow_se[i,]
      oW_upper = ow_est[i,]+qnorm(0.975)*ow_se[i,]
      ow_cp[i,] = 1*(oW_lower <= true_ow_val & true_ow_val <= oW_upper)
      
      ###IPTW method
      Delta.IPTW = IPTW.RMST(L=L,w,x,z,time,delta,px,ps.model,cen.trt.model,cen.con.model)
      iptw_est[i,] = Delta.IPTW$point
      iptw_se[i,] = sqrt(Delta.IPTW$variance)
      iptw_lower = iptw_est[i,]-qnorm(0.975)*iptw_se[i,]
      iptw_upper = iptw_est[i,]+qnorm(0.975)*iptw_se[i,]
      iptw_cp[i,] = 1*(iptw_lower <= true_iptw_val & true_iptw_val <= iptw_upper)
    }
    #average estimation
    ow_ave = colMeans(ow_est)
    iptw_ave = colMeans(iptw_est)
    
    #bias
    bias_ow = ow_ave-true_ow_val
    bias_iptw = iptw_ave-true_iptw_val
    
    #ase
    ase_ow = apply(ow_est,2,sd)
    ase_iptw = apply(iptw_est,2,sd)
    
    #average standarded error estimated from sandwich variance estimation
    esd_ow = colMeans(ow_se)
    esd_iptw = colMeans(iptw_se)
    
    #CP
    cover_ow = colMeans(ow_cp)
    cover_iptw = colMeans(iptw_cp)
    
    #summary results
    u1_est = data.frame(bias = c(bias_ow[1],bias_iptw[1]),
                        ase = c(ase_ow[1],ase_iptw[1]),
                        esd = c(esd_ow[1],esd_iptw[1]),
                        cp = c(cover_ow[1],cover_iptw[1]))
    
    u0_est = data.frame(bias = c(bias_ow[2],bias_iptw[2]),
                        ase = c(ase_ow[2],ase_iptw[2]),
                        esd = c(esd_ow[2],esd_iptw[2]),
                        cp = c(cover_ow[2],cover_iptw[2]))
    
    delta_est = data.frame(bias = c(bias_ow[3],bias_iptw[3]),
                           ase = c(ase_ow[3],ase_iptw[3]),
                           esd = c(esd_ow[3],esd_iptw[3]),
                           cp = c(cover_ow[3],cover_iptw[3]))
    res = cbind(u1_est,u0_est,delta_est)
    rownames(res) = c("ow","iptw")
    write.csv(res, paste0('lasso_sim_L_', L, '_G', gamma, '.csv'))
  }
}




