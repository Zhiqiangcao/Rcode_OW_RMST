###additional simulation results: when PS model is in m^{(a)}(x)
###specifically, 3*e(x) is included in m^{(1)}(x) and -e(x) is included in m^{(0)}(x) 

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
        true_ow_val<- c(0.794177915,1.854129826,-1.059951911)
        true_iptw_val <- c(0.799390562,1.8529407,-1.053550139)
        true_strim_val1 <- c(0.799389821,1.852941474,-1.053551654)
        true_strim_val2 <- c(0.799288051,1.85300063,-1.053712579)
        true_strim_val3 <- c(0.798360344,1.85336542,-1.055005076)
        true_atrim_val1 <- c(0.799384808,1.852944634,-1.053559826)
        true_atrim_val2 <- c(0.79240399,1.854852491,-1.062448501)
        true_atrim_val3 <- c(0.778459445,1.857269385,-1.078809941)
        true_trun_val1 <- true_trun_val2 <- true_trun_val3 <- c(0.799390562,1.8529407,-1.053550139)
      } else if (L == 5) {
        true_ow_val <- c(1.07731187,4.183468197,-3.106156328)
        true_iptw_val <- c(1.110376222,4.178919555,-3.068543333)
        true_strim_val1 <- c(1.110369986,4.178922271,-3.068552285)
        true_strim_val2 <- c(1.109609163,4.179141304,-3.069532141)
        true_strim_val3 <- c(1.103064036,4.180517851,-3.077453815)
        true_atrim_val1 <- c(1.110337565,4.178934257,-3.068596691)
        true_atrim_val2 <- c(1.066653821,4.186302894,-3.119649073)
        true_atrim_val3 <- c(0.985757709,4.195748711,-3.209991002)
        true_trun_val1 <- true_trun_val2 <- true_trun_val3 <- c(1.110376222,4.178919555,-3.068543333)
      } else if (L == 10) {
        true_ow_val <- c(1.211243725,7.206604566,-5.995360841)
        true_iptw_val <- c(1.273190646,7.198107742,-5.924917096)
        true_strim_val1 <- c(1.273174256,7.198112446,-5.92493819)
        true_strim_val2 <- c(1.271303649,7.198511063,-5.927207414)
        true_strim_val3 <- c(1.256310244,7.201059462,-5.944749218)
        true_atrim_val1 <- c(1.273093646,7.198135337,-5.925041691)
        true_atrim_val2 <- c(1.182848773,7.212335524,-6.029486751)
        true_atrim_val3 <- c(1.044070389,7.230610509,-6.18654012)
        true_trun_val1 <- true_trun_val2 <- true_trun_val3 <- c(1.273190646,7.198107742,-5.924917096)
      } 
    } else if (gamma == 3) {
      if (L == 2) {
        true_ow_val <- c(0.806918004,1.852328368,-1.045410363)
        true_iptw_val <- c(0.836755494,1.844750714,-1.007995219)
        true_strim_val1 <- c(0.822986458,1.849500344,-1.026513886)
        true_strim_val2 <- c(0.809411401,1.852586966,-1.043175566)
        true_strim_val3 <- c(0.797239925,1.854801566,-1.057561641)
        true_atrim_val1 <- c(0.836579582,1.844862803,-1.008283221)
        true_atrim_val2 <- c(0.811347695,1.852049832,-1.040702138)
        true_atrim_val3 <- c(0.781015427,1.857131368,-1.076115941)
        true_trun_val1 <- true_trun_val2 <- true_trun_val3 <- c(0.836755494,1.844750714,-1.007995219)
      } else if (L == 5) {
        true_ow_val <- c(1.149588489,4.175968455,-3.026379966)
        true_iptw_val <- c(1.343004924,4.146984042,-2.803979118)
        true_strim_val1 <- c(1.247483688,4.164932687,-2.917448999)
        true_strim_val2 <- c(1.159926334,4.176815939,-3.016889605)
        true_strim_val3 <- c(1.085063261,4.185453819,-3.100390559)
        true_atrim_val1 <- c(1.341773499,4.147398806,-2.805625307)
        true_atrim_val2 <- c(1.179058745,4.174964119,-2.995905373)
        true_atrim_val3 <- c(1.001972935,4.194950066,-3.192977131)
        true_trun_val1 <- true_trun_val2 <- true_trun_val3 <- c(1.343004924,4.146984042,-2.803979118)
      } else if (L == 10) {
        true_ow_val <- c(1.336809346,7.190987254,-5.854177909)
        true_iptw_val <- c(1.728923122,7.136478835,-5.407555713)
        true_strim_val1 <- c(1.511004012,7.169869587,-5.658865575)
        true_strim_val2 <- c(1.335563848,7.192348916,-5.856785068)
        true_strim_val3 <- c(1.200032974,7.208870403,-6.008837429)
        true_atrim_val1 <- c(1.725741684,7.137248146,-5.411506462)
        true_atrim_val2 <- c(1.373979145,7.18984192,-5.815862776)
        true_atrim_val3 <- c(1.066356062,7.228720016,-6.162363954)
        true_trun_val1 <- true_trun_val2 <- true_trun_val3 <- c(1.728923122,7.136478835,-5.407555713)
      } 
    } else if (gamma == 5) {
      if (L == 2) {
        true_ow_val <- c(0.806956972,1.852818422,-1.04586145)
        true_iptw_val <- c(0.854926773,1.840821039,-0.985894266)
        true_strim_val1 <- c(0.819915391,1.851040929,-1.031125539)
        true_strim_val2 <- c(0.803804842,1.854014664,-1.050209821)
        true_strim_val3 <- c(0.791402188,1.855938582,-1.064536394)
        true_atrim_val1 <- c(0.853869056,1.841454666,-0.98758561)
        true_atrim_val2 <- c(0.81048018,1.852702238,-1.042222058)
        true_atrim_val3 <- c(0.772339934,1.858310942,-1.085971007)
        true_trun_val1 <- true_trun_val2 <- true_trun_val3 <- c(0.854926773,1.840821039,-0.985894266)
      } else if (L == 5) {
        true_ow_val <-  c(1.146025362,4.177605691,-3.031580329)
        true_iptw_val <- c(1.456918931,4.131601384,-2.674682453)
        true_strim_val1 <- c(1.220930574,4.170498483,-2.949567909)
        true_strim_val2 <- c(1.121743774,4.182135832,-3.060392058)
        true_strim_val3 <- c(1.048379848,4.189756267,-3.141376419)
        true_atrim_val1 <- c(1.448820401,4.133939722,-2.685119321)
        true_atrim_val2 <- c(1.17086571,4.177259178,-3.006393468)
        true_atrim_val3 <- c(0.953548879,4.199576217,-3.246027337)
        true_trun_val1 <- true_trun_val2 <- true_trun_val3 <- c(1.456918931,4.131601384,-2.674682453)
      } else if (L == 10) {
        true_ow_val <- c(1.322850417,7.193449771,-5.870599354)
        true_iptw_val <- c(1.95686259,7.106666794,-5.149804205)
        true_strim_val1 <- c(1.442027825,7.179521095,-5.73749327)
        true_strim_val2 <- c(1.259409012,7.201894178,-5.942485166)
        true_strim_val3 <- c(1.136104204,7.2167153,-6.080611096)
        true_atrim_val1 <- c(1.935801934,7.110924551,-5.175122617)
        true_atrim_val2 <- c(1.349535427,7.19367623,-5.844140804)
        true_atrim_val3 <- c(0.993986017,7.237519614,-6.243533598)
        true_trun_val1 <- true_trun_val2 <- true_trun_val3 <- c(1.95686259,7.106666794,-5.149804205)
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
      #set.seed(100000+i)
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
      lambda_trt <- exp(as.numeric(cbind(1,x) %*% c(-1,0.4,0.2,0.1,-0.1,-0.2,-0.3))+3*e)
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
      
      #IPTW truncated - alpha = 0.025
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




