# Rcode_OW_RMST

The folders and the R code help to reproducing Tables 1-4 in the article "1.	Zhiqiang Cao, Lama Ghazi, Claudia Mastrogiacomo, Laura Forastiere, F. Perry Wilson, Fan Li. Using overlap weights to address extreme propensity scores in estimating restricted mean counterfactual survival times. American Journal of Epidemiology, forthcoming, 2024(+)." 

For questions or comments about the code, please contact Zhiqiang Cao zcaoae@connect.ust.hk. You will need to change the directory to use the example 
code in script. This folder includes the following functions:

1. cal_true_value_sim_manuscript.R and cal_true_value_appendix are to calculate true values of IPTW, OW, IPTW with symmetric trimming and asymmetric trimming weights, as well as IPTW with truncation wieghts under each case for simulations in paper, i.e., reproduce simulation results in Table 1 of paper and Web Table 5 of Web Appnedix.

2. generate_fig1.R is to plot figure 1 of manuscript.

3. RMST_functions_v3.R are related functions to calculate RMST and corresponding variance estimation based on IPTW, OW, IPTW with symmetric trimming 
and asymmetric trimming weights, as well as IPTW with truncation weights.

4. simu_main_fun_manuscript.R is to reproduce simulation results of Tables 2-4 in paper.
   
5. iptw_ow_bootstrap_variance_small_sample.R and iptw_ow_closed_variance_small_sample.R are to reproduce simulation results of Table 5 in paper.

6. simu_main_fun_appendix.R is to reproduce simulation results of Web Tables 6-8 in Web Appnedix.

7. simu_main_fun_misspecification_appendix.R is to reproduce simulation results of Web Tables 9-11 in Web Appnedix.
8. simu_main_appendix_not_select_variable.R and simu_main_appendix_post_lasso_select_variable.R are to reproduce simulation result of Web Table 12 in Web Appnedix.

9. surv.csv is the dataset used to demonstrate the proposed methods step-by-step in Web Appendix 3: R tutorial Section of Web Appnedix.

10. rhc.csv is the dataset used to analyze the RHC data set in Web Appendix 4: Right-Heart Catheterization Study Application.


