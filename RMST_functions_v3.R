#This version add truncated IPW
library(survival)

# 1. Define a function to calculate E{min(T^a, L)}, i.e., \mu_a=\int_0^L(\hat S^w_{a}(t)dt)
### Input: 1. Timevec: A vector of time (say, t)
###        2. L: restricted survival time (numerical value)
###        3. sf.est: estimated survival function (corresponds to sorted observed survival time) 
###        4. n: number of observations
### Output: restricted mean survival time up to L. 
cal.rmst = function(Timevec,L,sf.est,n){
  #sorted observed survival time
  us = Timevec[order(Timevec)]
  uss = c(0,us)
  survh = c(1,sf.est)
  nii = sum(uss <= L)
  auc = 0
  #using the following summation to approximate integral for calculating RMST
  #modify when nii==(length(us)+1)
  if(nii<(n+1)){
    for(j in 1:nii){
      difft = min(uss[j+1]-uss[j],L-uss[j])
      auc = auc + difft*survh[j]
    }
  }else{ #in this case, L=max(time)
    difft = diff(uss)
    auc = sum(difft*sf.est)
  }
  return (auc)
}

# 2. Define a function to calculate the predicted censoring score of z group
### Input: 1. Coxmodel: A Cox model object
###        2. x: covariates vector for censoring model
###        3. z: treatment binary variable
###        4. delta: censoring indicator
###        5. px: dimension of x
###        6. n: number of observations
### Output: Probability of P(c_i>u_j|X_i) for i=1,...,n and j=1,...,n 
###         note: observed survival time u has been sorted, i.e., u_1<u_2<...<u_n)
CensorScore.est = function(Coxmodel,x,z,delta,px,n){
  bc_est = Coxmodel$coef
  if(px == 1){ #only one covariate
    if(is.na(bc_est)){  #then no covariate is included in Cox model
      bc_est = 0  #reset its value to 0
    }
    s2 = cumsum(z[n:1]*exp(x[n:1]*bc_est))
  }
  else {s2 = cumsum(z[n:1]*exp(x[n:1,]%*%bc_est))}
  dlc = z*(1-delta)/s2[n:1]
  dlc[(1:length(dlc))[is.nan(dlc)]] = 0
  ssc = matrix(0, n, n)
  ssc[,1] = dlc[1]*exp(as.matrix(x)%*%bc_est)
  otc = outer(c(exp(x%*%bc_est)), dlc, "*")   
  otc[,1] = rep(0, n)
  otc_cum = t(apply(otc,1,cumsum))
  ssc_1 = outer(ssc[,1],rep(1,n))
  cumhcij = ssc_1 + otc_cum
  #for patient i with covariates x[i], compute estimated censoring score at observed time u[j]      
  ssc = exp(-cumhcij)
  return(ssc)
}

# 3. Define a function to calculate the estimated survival function for group z
### Input: 1. CensorScore: censoring score matrix (CensorScore[i,j] for unit i at u_j)
###        2. bw: balancing weight
###        3. z: treatment binary variable
###        4. delta: censoring indicator
###        5. n: number of observations
### Output: estimated survival function S(u_i) for i=1,2,...,n, where u_1<u_2<...<u_n
survf.est = function(CensorScore,bw,z,delta,n){
  ssc = CensorScore
  #lower triangular matrix
  ld = lower.tri(matrix(1, n, n), diag=TRUE)
  wij = z*bw
  denom0 = (wij/ssc)*ld
  #result for denominator of \Lambda_{a}^w(t) in (2), i.e., denom1[u_i], i=1,...,n 
  #where denom1[t]=\sum_{i=1}^n[I(A_i=a)w_{ia}(\hat{\beta_a})e^{\hat{\Lambda}_{ia}^C(u)}Y_{ia}(t)
  denom1 = apply(denom0,2,sum)
  #upper triangular matrix
  ud = matrix(1, n, n)-ld   
  #for patient i at time u_j
  num = wij*delta/ssc
  denom = matrix(rep(denom1,n),n,n)
  otr = num/denom
  #note that for patient i, \Lambda_{a_i}^w(u_j)=\Lambda_{a_i}^w(u_i) for all j>=i,
  cumhij0 = otr*ld+diag(otr)*ud
  udd = diag(n)+ud
  cumhij = cumhij0*udd
  cumhij[is.na(cumhij)] = 0
  cumhz = apply(cumhij,2,sum)
  #estimated survival function under 
  sf.est = exp(-cumhz) 
  res = list(s=sf.est,cumhij=cumhij)
  return(res)
}

# 4. Define a function to calculate \psi_{ia}(t) and \phi_{ia}
### Input: 1. time: A vector of time (say, t)
###        2. z: treatment binary variable
###        3. w: covariates vector used in PS model
###        4. delta: censoring indicator
###        5. bw: balancing weight
###        6. ps: propnesity score
###        7. sfc: survival function of censoring for individual i at time j, n*n matrix
###        8. sa: estimated survival function for z=a (a=0,1)
###        9. cumhij: cumulative hazard function for individual i at time j
###        10.Ebeta: second derivative of Taylor Expansion for PS model
###        11. Sbeta: first derivative of Taylor Expansion for PS model
###        12. n: number of observations
###        13. L: restricted survival time (numerical value)
### Output: \hat \phi_{ia} for i=1,2,...,n.
cal.psi.phi = function(time,z,w,delta,bw,ps,sfc,sa,cumhij,Ebeta,Sbeta,n,L){
  ssc = sfc
  ssc[ssc<10^{-14}] = 10^{-14}
  wij = z*bw  #balancing weights for z
  #since the second term of \psi_{ia}(t) can be decomposed into two parts,
  #first, we compute part 1, i.e., 
  #\int_0^t \frac{I(A_i=a)w_{ia}(\beta_a)e^{\Lambda_{ia}^C(u)}dN_{ia}(u)}{D_{a}(u;\beta_a,\theta_a)}
  ld = lower.tri(matrix(1, n, n), diag=TRUE)   # lower triangular matrix
  ud = matrix(1, n, n)-ld   # upper triangular matrix
  num1 = z*delta*wij/diag(ssc)
  num1 = diag(num1)+num1*ud  #n*n matrix
  denom1 = (z*wij/ssc)*ld
  denom1_colmean = apply(denom1,2,mean) #using mean since it is divided by n
  denom_colmean = matrix(rep(denom1_colmean,n),byrow=T,n)
  IM1 = num1/t(denom_colmean)  #n*n matrix
  IM1[is.na(IM1)] = 0
  
  #next, we compute part 2, i.e., 
  #\int_0^t \frac{I(A_i=a)\omega_{ia}(\beta_a)e^{\Lambda_{ia}^C(u)}Y_{ia}(u)d{\Lambda}_{a}(u)}{D_{a}(u;\beta_a,\theta_a)}
  lamhm = cumhij
  num2 = z*wij*lamhm/ssc  #n*n matrix
  denom2 = denom_colmean #n*n matrix
  otr = num2/denom2
  qm = t(apply(otr, 1, cumsum)) 
  IM2 = qm*ld+diag(qm)*ud
  IM2[is.na(IM2)] = 0
  #then, the second term of \psi_{ia}(t) is part2
  part2 = IM1-IM2  #n*n matrix
  
  #next, we compute first term of B_a(t;\beta_a,\theta_a)
  bnum1 = z*delta*wij*(-1)*(1-ps)/diag(ssc) #1*n vector
  bdenom1_sum = apply(denom1,2,sum) #1*n vector
  Ba1_m = as.numeric(bnum1/bdenom1_sum)*w #n*(p+1) matrix
  Ba1_cumsum = apply(Ba1_m,2,cumsum) #n*(p+1) matrix
  Ba1_cumsum[is.na(Ba1_cumsum)] = 0
  
  #then, we compute B_a^2(t;\beta_a,\theta_a)
  bnum2s = z*delta*wij/diag(ssc) #vector, length is n
  temp0 = (z*wij/ssc)*(1-ps)  #n*n matrix
  temp1 = temp0*ld             # lower triangular matrix
  bnum2f = apply(temp1,2,function(x){
    x = as.numeric(x);
    temp2 = x*w   #n*(p+1) matrix
    res = apply(temp2,2,sum)
  })  #(p+2)*n matrix
  bnum2f = t(bnum2f)
  Ba2_m = (bnum2s/bdenom1_sum^2)*bnum2f  #n*(p+1) matrix
  Ba2_cumsum = apply(Ba2_m,2,cumsum) #n*(p+1) matrix
  Ba2_cumsum[is.na(Ba2_cumsum)] = 0
  Ba = Ba1_cumsum + Ba2_cumsum
  
  #next, we compute part 1 of \psi(t)
  part1 = matrix(0,n,n)
  EBB_inv = solve(Ebeta)
  for(i in 1:n){
    part1i = Ba%*%EBB_inv%*%matrix(Sbeta[i,],ncol = 1) #n*1
    part1[i,] = part1i 
  }
  psi = part1 + part2
  
  #next, we compute phi_{ia}
  us = time  
  phi = numeric(n)
  for(i in 1:n){
    psi_i = psi[i,] #length is n
    survhi = sa*psi_i
    phi[i] = cal.rmst(us,L,survhi,n)
  }
  return(phi)
}


# overlap weight with asymptotic variance estimator
OW.RMST = function(L,w,x,z,time,delta,px,ps.model,cen.trt.model,cen.con.model){
  # w is covariates vector for PS model, including the intercept term
  # x is covariates vector for censoring model
  # if no covariate is included in censoring model, then x is a n*1 0 matrix
  # propensity score
  ps = 1/(1+exp(-c(w %*% ps.model$coefficients))) 
  # censoring function
  n = length(z)
  Kc.trt.est = CensorScore.est(cen.trt.model,x,z,delta,px,n)   #for treatment
  Kc.con.est = CensorScore.est(cen.con.model,x,1-z,delta,px,n) #for control
  # estimated survival function in treatment group
  bw.trt = 1-ps
  sf.trt.res = survf.est(Kc.trt.est,bw.trt,z,delta,n)
  surv.trt.est = sf.trt.res$s
  cumh.trt.est = sf.trt.res$cumhij
  auc1 = cal.rmst(time,L,surv.trt.est,n)
  
  # variance estimate
  # Taylor Expansion for PS model
  Ebeta = crossprod(sqrt(ps*(1-ps)) * w) / n
  Sbeta.trt = (z-ps)*w   #logistic score for treated 
  # for treatment group
  phi.trt = cal.psi.phi(time,z,w,delta,bw.trt,ps,Kc.trt.est,surv.trt.est,cumh.trt.est,Ebeta,Sbeta.trt,n,L)
  var_mu1 = mean(phi.trt^2)/n
  
  # estimated survival function in control group
  bw.con = ps
  sf.con.res = survf.est(Kc.con.est,bw.con,1-z,delta,n)
  surv.con.est = sf.con.res$s
  cumh.con.est = sf.con.res$cumhij
  auc0 = cal.rmst(time,L,surv.con.est,n)
  # for control group
  Sbeta.con = -(z-ps)*w   #logistic score for control
  phi.con = cal.psi.phi(time,1-z,w,delta,bw.con,1-ps,Kc.con.est,surv.con.est,cumh.con.est,Ebeta,Sbeta.con,n,L)
  var_mu0 = mean(phi.con^2)/n
  var_triangle = mean((phi.trt-phi.con)^2)/n
  
  # point estimate
  p.est = c(auc1,auc0,auc1-auc0)
  var.est = c(var_mu1,var_mu0,var_triangle)
  res = list(point=p.est,variance=var.est)
  return(res)
}


# IPTW with asymptotic variance estimator
IPTW.RMST = function(L,w,x,z,time,delta,px,ps.model,cen.trt.model,cen.con.model){
  # w is covariates vector for PS model, including the intercept term
  # x is covariates vector for censoring model
  # if no covariate is included in censoring model, then x is a n*1 0 matrix
  # propensity score
  ps = 1/(1+exp(-c(w %*% ps.model$coefficients))) 
  # censoring function
  n = length(z)
  Kc.trt.est = CensorScore.est(cen.trt.model,x,z,delta,px,n)   #for treatment
  Kc.con.est = CensorScore.est(cen.con.model,x,1-z,delta,px,n) #for control
  # estimated survival function in treatment group
  bw.trt = 1/ps
  sf.trt.res = survf.est(Kc.trt.est,bw.trt,z,delta,n)
  surv.trt.est = sf.trt.res$s
  cumh.trt.est = sf.trt.res$cumhij
  auc1 = cal.rmst(time,L,surv.trt.est,n)
  
  # variance estimate
  # Taylor Expansion for PS model
  Ebeta = crossprod(sqrt(ps*(1-ps)) * w) / n
  Sbeta.trt = (z-ps)*w   #logistic score for treated 
  # for treatment group
  phi.trt = cal.psi.phi(time,z,w,delta,bw.trt,ps,Kc.trt.est,surv.trt.est,cumh.trt.est,Ebeta,Sbeta.trt,n,L)
  var_mu1 = mean(phi.trt^2)/n
  
  # estimated survival function in control group
  bw.con = 1/(1-ps)
  sf.con.res = survf.est(Kc.con.est,bw.con,1-z,delta,n)
  surv.con.est = sf.con.res$s
  cumh.con.est = sf.con.res$cumhij
  auc0 = cal.rmst(time,L,surv.con.est,n)
  # for control group
  Sbeta.con = -(z-ps)*w   #logistic score for control
  phi.con = cal.psi.phi(time,1-z,w,delta,bw.con,1-ps,Kc.con.est,surv.con.est,cumh.con.est,Ebeta,Sbeta.con,n,L)
  var_mu0 = mean(phi.con^2)/n
  var_triangle = mean((phi.trt-phi.con)^2)/n
  
  # point estimate
  p.est = c(auc1,auc0,auc1-auc0)
  var.est = c(var_mu1,var_mu0,var_triangle)
  res = list(point=p.est,variance=var.est)
  return(res)
}


# IPTW symmetric trimming with asymptotic variance estimator
IPTW.sym.RMST = function(L,w,x,z,time,delta,px,q,ps.model,ps.formula,censor.formula){
  # w is covariates vector for PS model, including the intercept term
  # x is covariates vector for censoring model
  # if no covariate is included in censoring model, then x is a n*1 0 matrix
  # original propensity score
  ps = 1/(1+exp(-c(w %*% ps.model$coefficients)))
  keep = ((ps>= q) & (ps <= (1-q)))
  ptrim = 1 - mean(keep)
  # trim the dataset
  if(px==1) {x = x[keep]; x = matrix(x, ncol=px); 
  colnames(x)=paste("x",1:px,sep="")}
  else x = x[keep,]
  w = w[keep,]
  z = z[keep]
  time = time[keep]
  delta = delta[keep]
  data = as.data.frame(cbind(x,w,z,time,delta))
  data.trt = subset(data,z==1)
  data.con = subset(data,z==0)
  #update propensity score
  ps.model = glm(ps.formula,data=data,family=binomial(link="logit"))
  ps = 1/(1+exp(-c(w %*% ps.model$coefficients))) 
  #update censoring function
  cen.trt.model = coxph(censor.formula, data=data.trt)
  cen.con.model = coxph(censor.formula, data=data.con)
  
  # censoring function
  n = length(z)
  Kc.trt.est = CensorScore.est(cen.trt.model,x,z,delta,px,n)   #for treatment
  Kc.con.est = CensorScore.est(cen.con.model,x,1-z,delta,px,n) #for control
  # estimated survival function in treatment group
  bw.trt = 1/ps
  sf.trt.res = survf.est(Kc.trt.est,bw.trt,z,delta,n)
  surv.trt.est = sf.trt.res$s
  cumh.trt.est = sf.trt.res$cumhij
  auc1 = cal.rmst(time,L,surv.trt.est,n)
  
  # variance estimate
  # Taylor Expansion for PS model
  Ebeta = crossprod(sqrt(ps*(1-ps)) * w) / n
  Sbeta.trt = (z-ps)*w   #logistic score for treated 
  # for treatment group
  phi.trt = cal.psi.phi(time,z,w,delta,bw.trt,ps,Kc.trt.est,surv.trt.est,cumh.trt.est,Ebeta,Sbeta.trt,n,L)
  var_mu1 = mean(phi.trt^2)/n  #should use original 
  
  # estimated survival function in control group
  bw.con = 1/(1-ps)
  sf.con.res = survf.est(Kc.con.est,bw.con,1-z,delta,n)
  surv.con.est = sf.con.res$s
  cumh.con.est = sf.con.res$cumhij
  auc0 = cal.rmst(time,L,surv.con.est,n)
  # for control group
  Sbeta.con = -(z-ps)*w   #logistic score for control
  phi.con = cal.psi.phi(time,1-z,w,delta,bw.con,1-ps,Kc.con.est,surv.con.est,cumh.con.est,Ebeta,Sbeta.con,n,L)
  var_mu0 = mean(phi.con^2)/n
  var_triangle = mean((phi.trt-phi.con)^2)/n
  
  # point estimate
  p.est = c(auc1,auc0,auc1-auc0)
  var.est = c(var_mu1,var_mu0,var_triangle)
  res = list(point=p.est,variance=var.est,ptrim=ptrim)
  return(res)
}


# IPTW asymmetric trimming with asymptotic variance estimator
IPTW.asym.RMST = function(L,w,x,z,time,delta,px,q,ps.model,ps.formula,censor.formula){
  # w is covariates vector for PS model, including the intercept term
  # x is covariates vector for censoring model
  # original propensity score
  ps = 1/(1+exp(-c(w %*% ps.model$coefficients)))
  ps0 <- ps[z == 0]
  ps1 <- ps[z == 1]
  lps <- max(min(ps0), min(ps1))
  ups <- min(max(ps0), max(ps1))
  
  # PS Asymmetric trimming
  keep <- rep(NA, length(z))
  alpha0 <- as.numeric(quantile(ps0, 1-q))
  alpha1 <- as.numeric(quantile(ps1, q))
  keep[z == 0] <- ((ps0 >= alpha1) & (ps0 <= alpha0) & (ps0 >= lps) & (ps0 <= ups))
  keep[z == 1] <- ((ps1 >= alpha1) & (ps1 <= alpha0) & (ps1 >= lps) & (ps1 <= ups))
  ptrim <- 1 - mean(keep)
  
  # trim the dataset
  if(px==1) {x = x[keep]; x = matrix(x, ncol=px); 
  colnames(x)=paste("x",1:px,sep="")}
  else x = x[keep,]
  w = w[keep,]
  z = z[keep]
  time = time[keep]
  delta = delta[keep]
  data = as.data.frame(cbind(x,w,z,time,delta))
  data.trt = subset(data,z==1)
  data.con = subset(data,z==0)
  #update propensity score
  ps.model = glm(ps.formula,data=data,family=binomial(link="logit"))
  ps = 1/(1+exp(-c(w %*% ps.model$coefficients))) 
  #update censoring function
  cen.trt.model = coxph(censor.formula, data=data.trt)
  cen.con.model = coxph(censor.formula, data=data.con)
  
  # censoring function
  n = length(z)
  Kc.trt.est = CensorScore.est(cen.trt.model,x,z,delta,px,n)   #for treatment
  Kc.con.est = CensorScore.est(cen.con.model,x,1-z,delta,px,n) #for control
  # estimated survival function in treatment group
  bw.trt = 1/ps
  sf.trt.res = survf.est(Kc.trt.est,bw.trt,z,delta,n)
  surv.trt.est = sf.trt.res$s
  cumh.trt.est = sf.trt.res$cumhij
  auc1 = cal.rmst(time,L,surv.trt.est,n)
  
  # variance estimate
  # Taylor Expansion for PS model
  Ebeta = crossprod(sqrt(ps*(1-ps)) * w) / n
  Sbeta.trt = (z-ps)*w   #logistic score for treated 
  # for treatment group
  phi.trt = cal.psi.phi(time,z,w,delta,bw.trt,ps,Kc.trt.est,surv.trt.est,cumh.trt.est,Ebeta,Sbeta.trt,n,L)
  var_mu1 = mean(phi.trt^2)/n  #should use original 
  
  # estimated survival function in control group
  bw.con = 1/(1-ps)
  sf.con.res = survf.est(Kc.con.est,bw.con,1-z,delta,n)
  surv.con.est = sf.con.res$s
  cumh.con.est = sf.con.res$cumhij
  auc0 = cal.rmst(time,L,surv.con.est,n)
  # for control group
  Sbeta.con = -(z-ps)*w   #logistic score for control 
  phi.con = cal.psi.phi(time,1-z,w,delta,bw.con,1-ps,Kc.con.est,surv.con.est,cumh.con.est,Ebeta,Sbeta.con,n,L)
  var_mu0 = mean(phi.con^2)/n
  var_triangle = mean((phi.trt-phi.con)^2)/n
  
  # point estimate
  p.est = c(auc1,auc0,auc1-auc0)
  var.est = c(var_mu1,var_mu0,var_triangle)
  res = list(point=p.est,variance=var.est,ptrim=ptrim)
  return(res)
}

# IPTW truncation with asymptotic variance estimator
IPTW.trun.RMST = function(L,w,x,z,time,delta,px,q,ps.model,cen.trt.model,cen.con.model){
  # w is covariates vector for PS model, including the intercept term
  # x is covariates vector for censoring model
  # if no covariate is included in censoring model, then x is a n*1 0 matrix
  # propensity score
  ps = 1/(1+exp(-c(w %*% ps.model$coefficients))) 
  # censoring function
  n = length(z)
  Kc.trt.est = CensorScore.est(cen.trt.model,x,z,delta,px,n)   #for treatment
  Kc.con.est = CensorScore.est(cen.con.model,x,1-z,delta,px,n) #for control
  #truncated ps
  bounds = c(q, 1-q) #e.g., 0.1,0.9
  tps = ifelse(ps < quantile(ps, bounds)[1], quantile(ps, bounds)[1], ifelse(ps > quantile(ps, bounds)[2], quantile(ps, bounds)[2], ps))
  # estimated survival function in treatment group
  bw.trt = 1/tps
  sf.trt.res = survf.est(Kc.trt.est,bw.trt,z,delta,n)
  surv.trt.est = sf.trt.res$s
  cumh.trt.est = sf.trt.res$cumhij
  auc1 = cal.rmst(time,L,surv.trt.est,n)
  
  # variance estimate
  # Taylor Expansion for PS model
  Ebeta = crossprod(sqrt(ps*(1-ps)) * w) / n
  Sbeta.trt = (z-ps)*w   #logistic score for treated 
  # for treatment group
  phi.trt = cal.psi.phi(time,z,w,delta,bw.trt,ps,Kc.trt.est,surv.trt.est,cumh.trt.est,Ebeta,Sbeta.trt,n,L)
  var_mu1 = mean(phi.trt^2)/n
  
  # estimated survival function in control group
  bw.con = 1/(1-tps)
  sf.con.res = survf.est(Kc.con.est,bw.con,1-z,delta,n)
  surv.con.est = sf.con.res$s
  cumh.con.est = sf.con.res$cumhij
  auc0 = cal.rmst(time,L,surv.con.est,n)
  # for control group
  Sbeta.con = -(z-ps)*w   #logistic score for control
  phi.con = cal.psi.phi(time,1-z,w,delta,bw.con,1-ps,Kc.con.est,surv.con.est,cumh.con.est,Ebeta,Sbeta.con,n,L)
  var_mu0 = mean(phi.con^2)/n
  var_triangle = mean((phi.trt-phi.con)^2)/n
  
  # point estimate
  p.est = c(auc1,auc0,auc1-auc0)
  var.est = c(var_mu1,var_mu0,var_triangle)
  res = list(point=p.est,variance=var.est)
  return(res)
}

RMST.Effect.BW = function(Data,L,Treatment,SurvTime,Status,ps.formula,
                          censor.formula,Method="IPTW",alpha=0.1,q=0.01){
  # ps.formula = as.formula(paste(Treatment,"~",paste(Covariates,collapse="+"),sep=""))
  # censor.formula = as.formula(paste("Surv(",SurvTime,",I(1-",Status,"))","~",paste(Covariates,collapse="+"),sep=""))
  # sort data by observed survival time
  data.sort = Data[order(Data[,SurvTime]),]
  Data = data.sort
  w = model.matrix(ps.formula,data=Data)
  colnames(w)=paste("W",1:dim(w)[2],sep="")
  x = model.matrix(censor.formula, data=Data)
  px = dim(x)[2]; 
  n = dim(x)[1]
 if(px == 1){  #in this case, no covariate in censoring model, since x from model.matrix is a n*1 matrix 
    x = x       #(intercept) with all values=1, then estimate of x (intercept) will be NA in Cox model (censored model) 
  }else{
    x = matrix(x[,-1],ncol=px-1) #delete the first column of x, since that is intercept and no intercept in Cox model (censord model)
  }
  px = dim(x)[2]
  colnames(x)=paste("x",1:dim(x)[2],sep="")
  z = as.numeric(Data[,Treatment])
  time = as.numeric(Data[,SurvTime])
  delta = as.numeric(Data[,Status])
  # reconstruct the dataset
  data = as.data.frame(cbind(x,w,z,time,delta))
  data.trt = subset(data,z==1)
  data.con = subset(data,z==0)
  # update ps.formula
  ps.formula = as.formula(paste("z~",paste(colnames(w),collapse = "+"),"-1",sep="")) #since w[,1] is the intercept term
  ps.model = glm(ps.formula,data=data,family=binomial(link="logit"))  #logistic model default including intercept term
  # ps = 1/(1+exp(-c(w %*% ps.model$coefficients))) 
  # update censor.formula
  censor.formula = as.formula(paste("Surv(time, I(1-delta)) ~",paste(colnames(x),collapse = "+"),sep=""))
  cen.trt.model = coxph(censor.formula, data=data.trt)
  cen.con.model = coxph(censor.formula, data=data.con)
  if (Method=="IPTW"){
    res = IPTW.RMST(L=L,w,x,z,time,delta,px,ps.model,cen.trt.model,cen.con.model)
  }else if (Method=="OW"){
    res = OW.RMST(L=L,w,x,z,time,delta,px,ps.model,cen.trt.model,cen.con.model)
  }else if (Method=="Symmetric"){
    res = IPTW.sym.RMST(L=L,w,x,z,time,delta,px,q=alpha,ps.model,ps.formula,censor.formula)
  }else if (Method=="Asymmetric"){
    res = IPTW.asym.RMST(L=L,w,x,z,time,delta,px,q=q,ps.model,ps.formula,censor.formula)
  }else if (Method=="Truncation"){
    res = IPTW.trun.RMST(L,w,x,z,time,delta,px,q=q,ps.model,cen.trt.model,cen.con.model)
  }
  est = res$point
  se.est = sqrt(res$variance)
  CI.lower = est-qnorm(0.975)*se.est
  CI.upper = est+qnorm(0.975)*se.est
  output = data.frame(Estimate=est,SE=se.est,CI.lower=CI.lower,CI.upper=CI.upper)
  row.names(output) = c("mu1","mu0","Delta")
  return(output)
}
