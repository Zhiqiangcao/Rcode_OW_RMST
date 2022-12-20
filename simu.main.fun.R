###Reproduce simulation results of Tables 2-4 in paper 

rm(list=ls())
library(survival) 
library(mvtnorm)

#input source code
source("https://github.com/Zhiqiangcao/Rcode_OW_RMST/blob/main/RMST.functions.R")

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
censor.formula = as.formula(paste("Surv(",SurvTime,",I(1-",Status,"))","~ -1 + ",
                                  paste(Covariates,collapse="+"),sep=""))


###Main simulation
L_set = c(2,5,10)
gamma_set = c(1,3,5)
for(gamma in gamma_set){
  for(L in L_set){
    ###true value
    if (gamma == 1){
      if (L == 2) {
        true_iptw_val <- c(1.498194346,1.776465598,-0.278271251)
        true_ow_val<- c(1.503956526,1.776999254,-0.273042727)
        true_strim_val1 <- c(1.498196719,1.776466132,-0.278269413)
        true_strim_val2 <- c(1.498424552,1.776500247,-0.278075695)
        true_strim_val3 <- c(1.500031246,1.77668737,-0.276656125)
        true_atrim_val1 <- c(1.498205631,1.776468197,-0.278262566)
        true_atrim_val2 <- c(1.506414806,1.777409914,-0.270995108)
        true_atrim_val3 <- c(1.517891774,1.778439765,-0.260547992)
      } else if (L == 5) {
        true_iptw_val <- c(2.68834376,3.803739793,-1.115396033)
        true_ow_val <- c(2.693529714,3.805216113,-1.111686399)
        true_strim_val1 <- c(2.688345354,3.803741268,-1.115395914)
        true_strim_val2 <- c(2.688499696,3.803840495,-1.115340799)
        true_strim_val3 <- c(2.689695901,3.804373678,-1.114677777)
        true_atrim_val1 <- c(2.688349093,3.803748013,-1.11539892)
        true_atrim_val2 <- c(2.694271545,3.806528564,-1.112257019)
        true_atrim_val3 <- c(2.704467941,3.809467047,-1.104999106)
      } else if (L == 10) {
        true_iptw_val <- c(3.580969651,6.127871625,-2.546901973)
        true_ow_val <- c(3.564410608,6.129209984,-2.564799375)
        true_strim_val1 <- c(3.580964825,6.127873074,-2.546908249)
        true_strim_val2 <- c(3.580427372,6.127973884,-2.547546511)
        true_strim_val3 <- c(3.576302521,6.128484069,-2.552181548)
        true_atrim_val1 <- c(3.580937304,6.127882,-2.546944697)
        true_atrim_val2 <- c(3.555099884,6.131077118,-2.575977234)
        true_atrim_val3 <- c(3.516096433,6.134298449,-2.618202016)
      } 
    } else if (gamma == 3) {
      if (L == 2) {
        true_iptw_val <- c(1.498194346,1.776465598,-0.278271251)
        true_ow_val <- c(1.519824121,1.778325209,-0.258501088)
        true_strim_val1 <- c(1.513917095,1.777914328,-0.263997234)
        true_strim_val2 <- c(1.52174541,1.778499349,-0.256753939)
        true_strim_val3 <- c(1.526495526,1.778837162,-0.252341635)
        true_atrim_val1 <- c(1.498571888,1.776525057,-0.277953169)
        true_atrim_val2 <- c(1.519474667,1.778489147,-0.259014479)
        true_atrim_val3 <- c(1.529569611,1.779283215,-0.249713604)
      } else if (L == 5) {
        true_iptw_val <- c(2.68834376,3.803739793,-1.115396033)
        true_ow_val <- c(2.708913774,3.808817305,-1.099903531)
        true_strim_val1 <- c(2.702334097,3.807757913,-1.105423816)
        true_strim_val2 <- c(2.710466511,3.809293652,-1.09882714)
        true_strim_val3 <- c(2.715656448,3.81017585,-1.094519402)
        true_atrim_val1 <- c(2.688565935,3.803918419,-1.115352484)
        true_atrim_val2 <- c(2.706678715,3.809507676,-1.102828961)
        true_atrim_val3 <- c(2.717253173,3.81165812,-1.094404947)
      } else if (L == 10) {
        true_iptw_val <- c(3.580969651,6.127871625,-2.546901973)
        true_ow_val <- c(3.516710564,6.132360347,-2.615649783)
        true_strim_val1 <- c(3.536064442,6.131522959,-2.595458517)
        true_strim_val2 <- c(3.511535547,6.13276335,-2.621227802)
        true_strim_val3 <- c(3.496115636,6.133464037,-2.637348401)
        true_atrim_val1 <- c(3.579916,6.128079189,-2.548163189)
        true_atrim_val2 <- c(3.51328607,6.133944299,-2.620658229)
        true_atrim_val3 <- c(3.47837664,6.136021353,-2.657644713)
      } 
    } else if (gamma == 5) {
      if (L == 2) {
        true_iptw_val <- c(1.498194346,1.776465598,-0.278271251)
        true_ow_val <- c(1.52741699,1.778899222,-0.251482232)
        true_strim_val1 <- c(1.526182752,1.778820189,-0.252637437)
        true_strim_val2 <- c(1.530300707,1.779126143,-0.248825436)
        true_strim_val3 <- c(1.532413825,1.779282119,-0.246868294)
        true_atrim_val1 <- c(1.500702288,1.776757939,-0.276055651)
        true_atrim_val2 <- c(1.5278267,1.779077074,-0.251250374)
        true_atrim_val3 <- c(1.533809156,1.779501344,-0.245692188)
      } else if (L == 5) {
        true_iptw_val <-  c(2.68834376,3.803739793,-1.115396033)
        true_ow_val <- c(2.716895706,3.81034556,-1.093449854)
        true_strim_val1 <- c(2.715320014,3.810127794,-1.09480778)
        true_strim_val2 <- c(2.719969997,3.810937164,-1.090967167)
        true_strim_val3 <- c(2.722348138,3.811350763,-1.089002625)
        true_atrim_val1 <- c(2.690236597,3.804575459,-1.114338862)
        true_atrim_val2 <- c(2.715920785,3.811017628,-1.095096843)
        true_atrim_val3 <- c(2.722611959,3.812154499,-1.08954254)
      } else if (L == 10) {
        true_iptw_val <- c(3.580969651,6.127871625,-2.546901973)
        true_ow_val <- c(3.492736569,6.133647604,-2.640911036)
        true_strim_val1 <- c(3.497175208,6.133430802,-2.636255594)
        true_strim_val2 <- c(3.483648365,6.134136132,-2.650487767)
        true_strim_val3 <- c(3.476463733,6.134542136,-2.658078404)
        true_atrim_val1 <- c(3.574444164,6.128681929,-2.554237765)
        true_atrim_val2 <- c(3.48678927,6.135039375,-2.648250105)
        true_atrim_val3 <- c(3.466237274,6.136069457,-2.669832184)
      } 
    }
    
    #keep estimation results for each replicates
    iptw_est = ow_est = strim_est1 = strim_est2 = strim_est3 = atrim_est1 = atrim_est2 = atrim_est3 = matrix(0, ncol=3, nrow=nsim)
    iptw_se = ow_se = strim_se1 = strim_se2 = strim_se3 = atrim_se1 = atrim_se2 = atrim_se3 = matrix(0, ncol=3, nrow=nsim)
    iptw_cp = ow_cp = strim_cp1 = strim_cp2 = strim_cp3 = atrim_cp1 = atrim_cp2 = atrim_cp3 = matrix(0, ncol=3, nrow=nsim)
    
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
    
    #bias
    bias_ow = ow_ave-true_ow_val
    bias_iptw = iptw_ave-true_iptw_val
    bias_strim1 = strim1_ave-true_strim_val1
    bias_strim2 = strim2_ave-true_strim_val2
    bias_strim3 = strim3_ave-true_strim_val3
    bias_atrim1 = atrim1_ave-true_atrim_val1
    bias_atrim2 = atrim2_ave-true_atrim_val2
    bias_atrim3 = atrim3_ave-true_atrim_val3
    
    #ase
    ase_ow = apply(ow_est,2,sd)
    ase_iptw = apply(iptw_est,2,sd)
    ase_strim1 = apply(strim_est1,2,sd)
    ase_strim2 = apply(strim_est2,2,sd)
    ase_strim3 = apply(strim_est3,2,sd)
    ase_atrim1 = apply(atrim_est1,2,sd)
    ase_atrim2 = apply(atrim_est2,2,sd)
    ase_atrim3 = apply(atrim_est3,2,sd)
    
    #average standarded error estimated from sandwich variance estimation
    esd_ow = colMeans(ow_se)
    esd_iptw = colMeans(iptw_se)
    esd_strim1 = colMeans(strim_se1)
    esd_strim2 = colMeans(strim_se2)
    esd_strim3 = colMeans(strim_se3)
    esd_atrim1 = colMeans(atrim_se1)
    esd_atrim2 = colMeans(atrim_se2)
    esd_atrim3 = colMeans(atrim_se3)
    
    #CP
    cover_ow = colMeans(ow_cp)
    cover_iptw = colMeans(iptw_cp)
    cover_strim1 = colMeans(strim_cp1)
    cover_strim2 = colMeans(strim_cp2)
    cover_strim3 = colMeans(strim_cp3)
    cover_atrim1 = colMeans(atrim_cp1)
    cover_atrim2 = colMeans(atrim_cp2)
    cover_atrim3 = colMeans(atrim_cp3)
    
    #summary results
    u1_est = data.frame(bias = c(bias_ow[1],bias_iptw[1],bias_strim1[1],bias_strim2[1],
                                 bias_strim3[1],bias_atrim1[1],bias_atrim2[1],bias_atrim3[1]),
                        ase = c(ase_ow[1],ase_iptw[1],ase_strim1[1],ase_strim2[1],
                                ase_strim3[1],ase_atrim1[1],ase_atrim2[1],ase_atrim3[1]),
                        esd = c(esd_ow[1],esd_iptw[1],esd_strim1[1],esd_strim2[1],
                                esd_strim3[1],esd_atrim1[1],esd_atrim2[1],esd_atrim3[1]),
                        cp = c(cover_ow[1],cover_iptw[1],cover_strim1[1],cover_strim2[1],
                               cover_strim3[1],cover_atrim1[1],cover_atrim2[1],cover_atrim3[1]))
    
    u0_est = data.frame(bias = c(bias_ow[2],bias_iptw[2],bias_strim1[2],bias_strim2[2],
                                 bias_strim3[2],bias_atrim1[2],bias_atrim2[2],bias_atrim3[2]),
                        ase = c(ase_ow[2],ase_iptw[2],ase_strim1[2],ase_strim2[2],
                                ase_strim3[2],ase_atrim1[2],ase_atrim2[2],ase_atrim3[2]),
                        esd = c(esd_ow[2],esd_iptw[2],esd_strim1[2],esd_strim2[2],
                                esd_strim3[2],esd_atrim1[2],esd_atrim2[2],esd_atrim3[2]),
                        cp = c(cover_ow[2],cover_iptw[2],cover_strim1[2],cover_strim2[2],
                               cover_strim3[2],cover_atrim1[2],cover_atrim2[2],cover_atrim3[2]))
    
    delta_est = data.frame(bias = c(bias_ow[3],bias_iptw[3],bias_strim1[3],bias_strim2[3],
                                    bias_strim3[3],bias_atrim1[3],bias_atrim2[3],bias_atrim3[3]),
                           ase = c(ase_ow[3],ase_iptw[3],ase_strim1[3],ase_strim2[3],
                                   ase_strim3[3],ase_atrim1[3],ase_atrim2[3],ase_atrim3[3]),
                           esd = c(esd_ow[3],esd_iptw[3],esd_strim1[3],esd_strim2[3],
                                   esd_strim3[3],esd_atrim1[3],esd_atrim2[3],esd_atrim3[3]),
                           cp = c(cover_ow[3],cover_iptw[3],cover_strim1[3],cover_strim2[3],
                                  cover_strim3[3],cover_atrim1[3],cover_atrim2[3],cover_atrim3[3]))
    res = cbind(u1_est,u0_est,delta_est)
    rownames(res) = c("ow","iptw","strim1","strim2","strim3","astrim1","astrim2","astrim3")
    write.csv(res, paste0('sim_L_', L, '_G', gamma, '.csv'))
  }
}




