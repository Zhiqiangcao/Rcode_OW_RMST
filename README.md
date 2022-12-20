# Rcode_OW_RMST

The folders and the R code help to reproducing Tables 1-4 in the article "Addressing Extreme Propensity Scores in Estimating Restricted Mean Survival Time via the Overlap Weights" by Zhiqiang Cao and Fan Li (under review)

For questions or comments about the code, please contact Zhiqiang Cao zcaoae@connect.ust.hk. You will need to change the directory to use the example code in script. This folder includes the following functions:

1. cal.true.value.sim,R is to generate true values of IPTW, OW, IPTW with symmetric trimming and asymmetric trimming weights under each case for simulations in paper,
i.e., reproduce simulation results in Table 1 of paper,

2. generate.web.fig1.R is to plot web figure 1 of supplementary material of this paper.

3. RMST.functions.R are related functions to compute RMST and corresponding variance based on IPTW, OW, IPTW with symmetric trimming and asymmetric trimming weights

4. simu.main.fun.R is to reproduce simulation results of Tables 2-4 in paper.

5. surv.csv is the dataset used to demonstrate the proposed methods step-by-step in Web Appendix 3: R tutorial section of supplementary material. 

