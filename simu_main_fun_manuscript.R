###additional simulation results: when PS model is in m^{(a)}(x)
###specifically, 2*e(x) is included in m^{(1)}(x) and -e(x) is included in m^{(0)}(x) 

rm(list=ls())
library(survival) 
library(mvtnorm)

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
N <- 1000     #sample size

#Identity treatment, survival outcome, and censoring indicator
SurvTime = "time"
Status = "delta"
Treatment = "z"
#Identity column names of the pre-treatment covariates
Covariates = c("x1","x2","x3","x4","x5","x6")
#ps.formula and censor.formula
ps.formula = as.formula(paste(Treatment,"~",paste(Covariates,collapse="+"),sep=""))
censor.formula = as.formula(paste("Surv(",SurvTime,",I(1-",Status,"))","~",
                                  paste(Covariates,collapse="+"),sep=""))  #update, including intercept


###Main simulation
gamma_set = 1 #c(1,3,5)
L_set = c(2,5,10)
for(gamma in gamma_set){
  for(L in L_set){
    ###true value
    if (gamma == 1){
      if (L == 2) {
        true_ow_val<- c(1.026147949,1.854129826,-0.827981877)
        true_iptw_val <- c(1.024226465,1.8529407,-0.828714235)
        true_strim_val1 <- c(1.024226876,1.852941474,-0.828714598)
        true_strim_val2 <- c(1.024265025,1.85300063,-0.828735606)
        true_strim_val3 <- c(1.024579663,1.85336542,-0.828785757)
        true_atrim_val1 <- c(1.024227261,1.852944634,-0.828717373)
        true_atrim_val2 <- c(1.02577719,1.854852491,-0.829075301)
        true_atrim_val3 <- c(1.029709501,1.857269385,-0.827559884)
        true_trun_val1 <- true_trun_val2 <- true_trun_val3 <- c(1.024226465,1.8529407,-0.828714235)
      } else if (L == 5) {
        true_ow_val <- c(1.496617516,4.183468197,-2.686850681)
        true_iptw_val <- c(1.51817085,4.178919555,-2.660748705)
        true_strim_val1 <- c(1.51816674,4.178922271,-2.660755531)
        true_strim_val2 <- c(1.517663618,4.179141304,-2.661477686)
        true_strim_val3 <- c(1.513369589,4.180517851,-2.667148263)
        true_atrim_val1 <- c(1.518144255,4.178934257,-2.660790002)
        true_atrim_val2 <- c(1.489172895,4.186302894,-2.697129999)
        true_atrim_val3 <- c(1.435790661,4.195748711,-2.75995805)
        true_trun_val1 <- true_trun_val2 <- true_trun_val3 <- c(1.51817085,4.178919555,-2.660748705)
      } else if (L == 10) {
        true_ow_val <- c(1.741371978,7.206604566,-5.465232588)
        true_iptw_val <- c(1.796519987,7.198107742,-5.401587755)
        true_strim_val1 <- c(1.796506271,7.198112446,-5.401606175)
        true_strim_val2 <- c(1.794940276,7.198511063,-5.403570788)
        true_strim_val3 <- c(1.782313578,7.201059462,-5.418745884)
        true_atrim_val1 <- c(1.796437866,7.198135337,-5.401697472)
        true_atrim_val2 <- c(1.718429096,7.212335524,-5.493906428)
        true_atrim_val3 <- c(1.591041687,7.230610509,-5.639568822)
        true_trun_val1 <- true_trun_val2 <- true_trun_val3 <- c(1.796519987,7.198107742,-5.401587755)
      } 
    } else if (gamma == 3) {
      if (L == 2) {
        true_ow_val <- c(1.025190601,1.852328368,-0.827137767)
        true_iptw_val <- c(1.015788454,1.844750714,-0.82896226)
        true_strim_val1 <- c(1.02009592,1.849500344,-0.829404424)
        true_strim_val2 <- c(1.024297757,1.852586966,-0.82828921)
        true_strim_val3 <- c(1.028153232,1.854801566,-0.826648335)
        true_atrim_val1 <- c(1.015820826,1.844862803,-0.829041977)
        true_atrim_val2 <- c(1.021742736,1.852049832,-0.830307096)
        true_atrim_val3 <- c(1.029885995,1.857131368,-0.827245373)
        true_trun_val1 <- true_trun_val2 <- true_trun_val3 <- c(1.015788454,1.844750714,-0.82896226)
      } else if (L == 5) {
        true_ow_val <- c(1.51728134,4.175968455,-2.658687115)
        true_iptw_val <- c(1.636072016,4.146984042,-2.510912026)
        true_strim_val1 <- c(1.574381108,4.164932687,-2.590551579)
        true_strim_val2 <- c(1.521261587,4.176815939,-2.655554352)
        true_strim_val3 <- c(1.477460061,4.185453819,-2.707993758)
        true_atrim_val1 <- c(1.635201052,4.147398806,-2.512197754)
        true_atrim_val2 <- c(1.531431409,4.174964119,-2.64353271)
        true_atrim_val3 <- c(1.427743451,4.194950066,-2.767206616)
        true_trun_val1 <- true_trun_val2 <- true_trun_val3 <-  c(1.636072016,4.146984042,-2.510912026)
      } else if (L == 10) {
        true_ow_val <- c(1.786670443,7.190987254,-5.404316811)
        true_iptw_val <- c(2.098935467,7.136478835,-5.037543368)
        true_strim_val1 <- c(1.923080093,7.169869587,-5.246789495)
        true_strim_val2 <- c(1.785885274,7.192348916,-5.406463642)
        true_strim_val3 <- c(1.680126981,7.208870403,-5.528743422)
        true_atrim_val1 <- c(2.09621628,7.137248146,-5.041031866)
        true_atrim_val2 <- c(1.814692921,7.18984192,-5.375149)
        true_atrim_val3 <- c(1.572418889,7.228720016,-5.656301127)
        true_trun_val1 <- true_trun_val2 <- true_trun_val3 <- c(2.098935467,7.136478835,-5.037543368)
      } 
    } else if (gamma == 5) {
      if (L == 2) {
        true_ow_val <- c(1.026468112,1.852818422,-0.826350311)
        true_iptw_val <- c(1.011756677,1.840821039,-0.829064362)
        true_strim_val1 <- c(1.022597397,1.851040929,-0.828443532)
        true_strim_val2 <- c(1.027402471,1.854014664,-0.826612193)
        true_strim_val3 <- c(1.031119087,1.855938582,-0.824819495)
        true_atrim_val1 <- c(1.012131341,1.841454666,-0.829323325)
        true_atrim_val2 <- c(1.023019912,1.852702238,-0.829682326)
        true_atrim_val3 <- c(1.03319468,1.858310942,-0.825116262)
        true_trun_val1 <- true_trun_val2 <- true_trun_val3 <- c(1.011756677,1.840821039,-0.829064362)
      } else if (L == 5) {
        true_ow_val <-  c(1.504912627,4.177605691,-2.672693064)
        true_iptw_val <- c(1.693409459,4.131601384,-2.438191925)
        true_strim_val1 <- c(1.544690571,4.170498483,-2.625807912)
        true_strim_val2 <- c(1.4885466,4.182135832,-2.693589232)
        true_strim_val3 <- c(1.448582155,4.189756267,-2.741174112)
        true_atrim_val1 <- c(1.687670115,4.133939722,-2.446269607)
        true_atrim_val2 <- c(1.514432164,4.177259178,-2.662827014)
        true_atrim_val3 <- c(1.395845438,4.199576217,-2.803730778)
        true_trun_val1 <- true_trun_val2 <- true_trun_val3 <- c(1.693409459,4.131601384,-2.438191925)
      } else if (L == 10) {
        true_ow_val <- c(1.751763093,7.193449771,-5.441686678)
        true_iptw_val <- c(2.244921694,7.106666794,-4.8617451)
        true_strim_val1 <- c(1.838784064,7.179521095,-5.340737031)
        true_strim_val2 <- c(1.703633267,7.201894178,-5.498260912)
        true_strim_val3 <- c(1.612641988,7.2167153,-5.604073312)
        true_atrim_val1 <- c(2.226819178,7.110924551,-4.884105373)
        true_atrim_val2 <- c(1.768692921,7.19367623,-5.424983309)
        true_atrim_val3 <- c(1.504094407,7.237519614,-5.733425208)
        true_trun_val1 <- true_trun_val2 <- true_trun_val3 <- c(2.244921694,7.106666794,-4.8617451)
      } 
    }
    
    #keep estimation results for each replicates
    iptw_est = ow_est = strim_est1 = strim_est2 = strim_est3 = atrim_est1 = atrim_est2 = atrim_est3 = trun_est1 = trun_est2 = trun_est3 = matrix(0, ncol=3, nrow=nsim)
    iptw_se = ow_se = strim_se1 = strim_se2 = strim_se3 = atrim_se1 = atrim_se2 = atrim_se3 = trun_se1 = trun_se2 = trun_se3 = matrix(0, ncol=3, nrow=nsim)
    iptw_cp = ow_cp = strim_cp1 = strim_cp2 = strim_cp3 = atrim_cp1 = atrim_cp2 = atrim_cp3 = trun_cp1 = trun_cp2 = trun_cp3 = matrix(0, ncol=3, nrow=nsim)
    
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
      colnames(x) = Covariates
      # selection probability for trial
      ac <- gamma*c(0.15, 0.3, 0.3)
      ab <- -gamma*c(0.2, 0.25, 0.25)
      a <- c(a0, ac, ab)
      e <- plogis(c(cbind(1,x)%*%a))
      z <- rbinom(N, 1, e)
      lambda_trt <- exp(as.numeric(cbind(1,x) %*% c(-1,0.4,0.2,0.1,-0.1,-0.2,-0.3))+2*e)
      lambda_crl <- exp(as.numeric(cbind(1,x) %*% c(-1.4,0.0,-0.2,-0.3,-0.5,-0.6,-0.7))-e)
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
      Data = data.frame(time=time,delta=delta,z=z,x)
      
      ###calculate 8 estimators
      #OW
      Delta.OW = RMST.Effect.BW(Data,L,Treatment,SurvTime,Status,ps.formula,
                                censor.formula,Method="OW")
      ow_est[i,] = Delta.OW$Estimate
      ow_se[i,] = Delta.OW$SE
      oW_lower = Delta.OW$CI.lower
      oW_upper = Delta.OW$CI.upper
      ow_cp[i,] = 1*(oW_lower <= true_ow_val & true_ow_val <= oW_upper)
      
      #IPTW
      Delta.IPTW = RMST.Effect.BW(Data,L,Treatment,SurvTime,Status,ps.formula,
                                  censor.formula,Method="IPTW")
      iptw_est[i,] = Delta.IPTW$Estimate
      iptw_se[i,] = Delta.IPTW$SE
      iptw_lower = Delta.IPTW$CI.lower
      iptw_upper = Delta.IPTW$CI.upper
      iptw_cp[i,] = 1*(iptw_lower <= true_iptw_val & true_iptw_val <= iptw_upper)
      
      #symmetric trimming - 0.05
      Delta.SW1 = RMST.Effect.BW(Data,L,Treatment,SurvTime,Status,ps.formula,
                                 censor.formula,Method="Symmetric",alpha=0.05)
      strim_est1[i,] = Delta.SW1$Estimate
      strim_se1[i,] = Delta.SW1$SE
      strim_lower1 = Delta.SW1$CI.lower
      strim_upper1 = Delta.SW1$CI.upper
      strim_cp1[i,] = 1*(strim_lower1 <= true_strim_val1 & true_strim_val1 <= strim_upper1)
      
      #symmetric trimming - 0.10
      Delta.SW2 = RMST.Effect.BW(Data,L,Treatment,SurvTime,Status,ps.formula,
                                 censor.formula,Method="Symmetric",alpha=0.1)
      strim_est2[i,] = Delta.SW2$Estimate
      strim_se2[i,] = Delta.SW2$SE
      strim_lower2 = Delta.SW2$CI.lower
      strim_upper2 = Delta.SW2$CI.upper
      strim_cp2[i,] = 1*(strim_lower2 <= true_strim_val2 & true_strim_val2 <= strim_upper2)
      
      #symmetric trimming - 0.15
      Delta.SW3 = RMST.Effect.BW(Data,L,Treatment,SurvTime,Status,ps.formula,
                                 censor.formula,Method="Symmetric",alpha=0.15)
      strim_est3[i,] = Delta.SW3$Estimate
      strim_se3[i,] = Delta.SW3$SE
      strim_lower3 = Delta.SW3$CI.lower
      strim_upper3 = Delta.SW3$CI.upper
      strim_cp3[i,] = 1*(strim_lower3 <= true_strim_val3 & true_strim_val3 <= strim_upper3)
      
      #asymmetric trimming - q = 0
      Delta.AW1 = RMST.Effect.BW(Data,L,Treatment,SurvTime,Status,ps.formula,
                                 censor.formula,Method="Asymmetric",q=0) 
      atrim_est1[i,] = Delta.AW1$Estimate
      atrim_se1[i,] = Delta.AW1$SE
      atrim_lower1 = Delta.AW1$CI.lower
      atrim_upper1 = Delta.AW1$CI.upper
      atrim_cp1[i,] = 1*(atrim_lower1 <= true_atrim_val1 & true_atrim_val1 <= atrim_upper1)
      
      #asymmetric trimming - q = 0.01
      Delta.AW2 = RMST.Effect.BW(Data,L,Treatment,SurvTime,Status,ps.formula,
                                 censor.formula,Method="Asymmetric",q=0.01) 
      atrim_est2[i,] = Delta.AW2$Estimate
      atrim_se2[i,] = Delta.AW2$SE
      atrim_lower2 = Delta.AW2$CI.lower
      atrim_upper2 = Delta.AW2$CI.upper
      atrim_cp2[i,] = 1*(atrim_lower2 <= true_atrim_val2 & true_atrim_val2 <= atrim_upper2)
      
      #asymmetric trimming - q = 0.05
      Delta.AW3 = RMST.Effect.BW(Data,L,Treatment,SurvTime,Status,ps.formula,
                                 censor.formula,Method="Asymmetric",q=0.05) 
      atrim_est3[i,] = Delta.AW3$Estimate
      atrim_se3[i,] = Delta.AW3$SE
      atrim_lower3 = Delta.AW3$CI.lower
      atrim_upper3 = Delta.AW3$CI.upper
      atrim_cp3[i,] = 1*(atrim_lower3 <= true_atrim_val3 & true_atrim_val3 <= atrim_upper3)
      
      #IPTW trnncated - alpha = 0.025
      Delta.trun1 = RMST.Effect.BW(Data,L,Treatment,SurvTime,Status,ps.formula,
                                   censor.formula,Method="Truncated",alpha=0.025)
      trun_est1[i,] = Delta.trun1$Estimate
      trun_se1[i,] = Delta.trun1$SE
      trun_lower1 = Delta.trun1$CI.lower
      trun_upper1 = Delta.trun1$CI.upper
      trun_cp1[i,] = 1*(trun_lower1 <= true_trun_val1 & true_trun_val1 <= trun_upper1)
      
      #IPTW trnncated - alpha = 0.05
      Delta.trun2 = RMST.Effect.BW(Data,L,Treatment,SurvTime,Status,ps.formula,
                                   censor.formula,Method="Truncated",alpha=0.05)
      trun_est2[i,] = Delta.trun2$Estimate
      trun_se2[i,] = Delta.trun2$SE
      trun_lower2 = Delta.trun2$CI.lower
      trun_upper2 = Delta.trun2$CI.upper
      trun_cp2[i,] = 1*(trun_lower2 <= true_trun_val2 & true_trun_val2 <= trun_upper2)
      
      #IPTW trnncated - alpha = 0.1
      Delta.trun3 = RMST.Effect.BW(Data,L,Treatment,SurvTime,Status,ps.formula,
                                   censor.formula,Method="Truncated",alpha=0.1)
      trun_est3[i,] = Delta.trun3$Estimate
      trun_se3[i,] = Delta.trun3$SE
      trun_lower3 = Delta.trun3$CI.lower
      trun_upper3 = Delta.trun3$CI.upper
      trun_cp3[i,] = 1*(trun_lower3 <= true_trun_val3 & true_trun_val3 <= trun_upper3)
    }
    
    #average estimation
    ow_ave = colMeans(ow_est)
    iptw_ave = colMeans(iptw_est)
    strim1_ave = colMeans(strim_est1)
    strim2_ave = colMeans(strim_est2)
    strim3_ave = colMeans(strim_est3)
    atrim1_ave = colMeans(atrim_est1)
    atrim2_ave = colMeans(atrim_est2)
    atrim3_ave = colMeans(atrim_est3)
    trun1_ave = colMeans(trun_est1)
    trun2_ave = colMeans(trun_est2)
    trun3_ave = colMeans(trun_est3)
    
    #bias
    bias_ow = ow_ave-true_ow_val
    bias_iptw = iptw_ave-true_iptw_val
    bias_strim1 = strim1_ave-true_strim_val1
    bias_strim2 = strim2_ave-true_strim_val2
    bias_strim3 = strim3_ave-true_strim_val3
    bias_atrim1 = atrim1_ave-true_atrim_val1
    bias_atrim2 = atrim2_ave-true_atrim_val2
    bias_atrim3 = atrim3_ave-true_atrim_val3
    bias_trun1 = trun1_ave-true_trun_val1
    bias_trun2 = trun2_ave-true_trun_val2
    bias_trun3 = trun3_ave-true_trun_val3
    
    #ase
    ase_ow = apply(ow_est,2,sd)
    ase_iptw = apply(iptw_est,2,sd)
    ase_strim1 = apply(strim_est1,2,sd)
    ase_strim2 = apply(strim_est2,2,sd)
    ase_strim3 = apply(strim_est3,2,sd)
    ase_atrim1 = apply(atrim_est1,2,sd)
    ase_atrim2 = apply(atrim_est2,2,sd)
    ase_atrim3 = apply(atrim_est3,2,sd)
    ase_trun1 = apply(trun_est1,2,sd)
    ase_trun2 = apply(trun_est2,2,sd)
    ase_trun3 = apply(trun_est3,2,sd)
    
    #average standarded error estimated from sandwich variance estimation
    esd_ow = colMeans(ow_se)
    esd_iptw = colMeans(iptw_se)
    esd_strim1 = colMeans(strim_se1)
    esd_strim2 = colMeans(strim_se2)
    esd_strim3 = colMeans(strim_se3)
    esd_atrim1 = colMeans(atrim_se1)
    esd_atrim2 = colMeans(atrim_se2)
    esd_atrim3 = colMeans(atrim_se3)
    esd_trun1 = colMeans(trun_se1)
    esd_trun2 = colMeans(trun_se2)
    esd_trun3 = colMeans(trun_se3)
    
    #CP
    cover_ow = colMeans(ow_cp)
    cover_iptw = colMeans(iptw_cp)
    cover_strim1 = colMeans(strim_cp1)
    cover_strim2 = colMeans(strim_cp2)
    cover_strim3 = colMeans(strim_cp3)
    cover_atrim1 = colMeans(atrim_cp1)
    cover_atrim2 = colMeans(atrim_cp2)
    cover_atrim3 = colMeans(atrim_cp3)
    cover_trun1 = colMeans(trun_cp1)
    cover_trun2 = colMeans(trun_cp2)
    cover_trun3 = colMeans(trun_cp3)
    
    #summary results
    u1_est = data.frame(bias = c(bias_ow[1],bias_iptw[1],bias_strim1[1],bias_strim2[1],
                                 bias_strim3[1],bias_atrim1[1],bias_atrim2[1],bias_atrim3[1],
                                 bias_trun1[1],bias_trun2[1],bias_trun3[1]),
                        ase = c(ase_ow[1],ase_iptw[1],ase_strim1[1],ase_strim2[1],
                                ase_strim3[1],ase_atrim1[1],ase_atrim2[1],ase_atrim3[1],
                                ase_trun1[1],ase_trun2[1],ase_trun3[1]),
                        esd = c(esd_ow[1],esd_iptw[1],esd_strim1[1],esd_strim2[1],
                                esd_strim3[1],esd_atrim1[1],esd_atrim2[1],esd_atrim3[1],
                                esd_trun1[1],esd_trun2[1],esd_trun3[1]),
                        cp = c(cover_ow[1],cover_iptw[1],cover_strim1[1],cover_strim2[1],
                               cover_strim3[1],cover_atrim1[1],cover_atrim2[1],cover_atrim3[1],
                               cover_trun1[1],cover_trun2[1],cover_trun3[1]))
    
    u0_est = data.frame(bias = c(bias_ow[2],bias_iptw[2],bias_strim1[2],bias_strim2[2],
                                 bias_strim3[2],bias_atrim1[2],bias_atrim2[2],bias_atrim3[2],
                                 bias_trun1[2],bias_trun2[2],bias_trun3[2]),
                        ase = c(ase_ow[2],ase_iptw[2],ase_strim1[2],ase_strim2[2],
                                ase_strim3[2],ase_atrim1[2],ase_atrim2[2],ase_atrim3[2],
                                ase_trun1[2],ase_trun2[2],ase_trun3[2]),
                        esd = c(esd_ow[2],esd_iptw[2],esd_strim1[2],esd_strim2[2],
                                esd_strim3[2],esd_atrim1[2],esd_atrim2[2],esd_atrim3[2],
                                esd_trun1[2],esd_trun2[2],esd_trun3[2]),
                        cp = c(cover_ow[2],cover_iptw[2],cover_strim1[2],cover_strim2[2],
                               cover_strim3[2],cover_atrim1[2],cover_atrim2[2],cover_atrim3[2],
                               cover_trun1[2],cover_trun2[2],cover_trun3[2]))
    
    delta_est = data.frame(bias = c(bias_ow[3],bias_iptw[3],bias_strim1[3],bias_strim2[3],
                                    bias_strim3[3],bias_atrim1[3],bias_atrim2[3],bias_atrim3[3],
                                    bias_trun1[3],bias_trun2[3],bias_trun3[3]),
                           ase = c(ase_ow[3],ase_iptw[3],ase_strim1[3],ase_strim2[3],
                                   ase_strim3[3],ase_atrim1[3],ase_atrim2[3],ase_atrim3[3],
                                   ase_trun1[3],ase_trun2[3],ase_trun3[3]),
                           esd = c(esd_ow[3],esd_iptw[3],esd_strim1[3],esd_strim2[3],
                                   esd_strim3[3],esd_atrim1[3],esd_atrim2[3],esd_atrim3[3],
                                   esd_trun1[3],esd_trun2[3],esd_trun3[3]),
                           cp = c(cover_ow[3],cover_iptw[3],cover_strim1[3],cover_strim2[3],
                                  cover_strim3[3],cover_atrim1[3],cover_atrim2[3],cover_atrim3[3],
                                  cover_trun1[3],cover_trun2[3],cover_trun3[3]))
    res = cbind(u1_est,u0_est,delta_est)
    rownames(res) = c("ow","iptw","strim1","strim2","strim3","astrim1","astrim2","astrim3",
                      "trun1","trun2","trun3")
    write.csv(res, paste0('sim_L_', L, '_G', gamma, '.csv'))
  }
}




