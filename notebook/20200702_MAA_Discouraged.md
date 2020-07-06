 # 7.6. Update: 
 
 Here are the plots showing how bootstrap and permutation compare in terms of false positives and negatives. I've summed up the results in the below table: (Total For each column = 3060)
 
 |    | Covariance Bootstrap | Covariance Permutation | GxE Bootstrap | GxE Permutation |
 | ---|---|---|---|---|
 | True Positive | 2255 | 390 | 2520| 2260 | 
 | True Negative | 425 | 570 | 70 | 540|
 | False Positive | 145 | 0 | 470 | 0 | 
 | False Negative | 235 | 2100 | 0| 260 | 
 
Looks like our approach is pretty conservative given the high number of false negatives for permutation.... I see now what you mean by altering our inference... 

### Covariance Bootstrap
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.6.CovarianceBootCheck.png)

### Covariance Permutation
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.6.CovariancePErmutationcheck.png)

### GxE Bootstrap
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.6.GxEBootstrapCheck.png)

### GxE Permutation
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.6.GxEPermcheck.png)


# Previous Post from 7.1
Well, I've encountered another problem in the sims, and I'm starting to wonder how its possible to have been working on a project for over a year and still be so... wrong. Am I just going in circles? Each edit creates a new error? Am I sisyphus, doomed to my coding purgatory for all time? 

The issue now is with the confidence intervals and pvalues. As you said, confidence intervals that do not cross zero should be significant. Except we have a whole bunch that aren't significant that should be: 

### Covariance

Shown are estimates with 95% CIs from a single replicate. Each row is a unique set of parameters
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.2.CovConfint_discrep.png)

### GxE

This one has its own set of problems: in the top left, the damn GxE estimates don't overlap with the confidence intervals?! WTF. 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.2.GxEConfint_discreps.png)


### Estimates with CI's and True estimate

I already looked at whether confidence intervals overlap with true values. I found that no, they don't, but at least with covariance it is close enough it doesn't seem to be alarming. Here are the plots again, this time with "x" as the estimate (colors signify std. deviation), and red circle as the true value. 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.2.ConfintOverlap_true.png)

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.2.GxEoverlapper_true.png)

From that plot of GxE, it is clear that its just the first sets of parameters that generate such wonky CIs/estimates. Here is a subset of the parameter sets that generate those weird values:
One clear pattern is that these occur in the n_pop = 10 scenario. But why would that drive this pattern while other n_pop values don't?

```{r}
 delta_env delta_gen sample_size n_pop std_dev interaction true_GxE_emm GxE_emm GxE_emm_lwrCI GxE_emm_uprCI
0.01	-1.00	 5	10	0.5	0	0	0.05	0.06	0.08
0.50	-1.00	 5	10	0.5	0	0	0.05	0.06	0.07
1.00	-1.00	 5	10	0.5	0	0	0.04	0.04	0.06
0.01	0.01	 5	10	0.5	0	0	0.33	0.39	0.49
0.50	0.01	 5	10	0.5	0	0	0.10	0.12	0.16
1.00	0.01	 5	10	0.5	0	0	0.05	0.06	0.08
0.01	1.00	 5	10	0.5	0	0	0.05	0.06	0.08
0.50	1.00	 5	10	0.5	0	0	0.05	0.06	0.07
1.00	1.00	 5	10	0.5	0	0	0.04	0.04	0.06
0.01	-1.00	10	10	0.5	0	0	0.04	0.04	0.06
0.50	-1.00	10	10	0.5	0	0	0.04	0.04	0.05
1.00	-1.00	10	10	0.5	0	0	0.03	0.03	0.04
0.01	0.01	10	10	0.5	0	0	0.24	0.27	0.36
0.50	0.01	10	10	0.5	0	0	0.08	0.09	0.11
1.00	0.01	10	10	0.5	0	0	0.04	0.04	0.06
0.01	1.00	10	10	0.5	0	0	0.04	0.04	0.06
0.50	1.00	10	10	0.5	0	0	0.04	0.04	0.05
1.00	1.00	10	10	0.5	0	0	0.03	0.03	0.04
0.01	-1.00	 5	10	1.5	0	0	0.14	0.17	0.21
0.50	-1.00	 5	10	1.5	0	0	0.13	0.15	0.20
1.00	-1.00	 5	10	1.5	0	0	0.11	0.13	0.16
0.01	0.01	 5	10	1.5	0	0	0.33	0.39	0.49
0.50	0.01	 5	10	1.5	0	0	0.23	0.27	0.35
1.00	0.01	 5	10	1.5	0	0	0.15	0.17	0.22
0.01	1.00	 5	10	1.5	0	0	0.15	0.17	0.22
0.50	1.00	 5	10	1.5	0	0	0.13	0.16	0.20
1.00	1.00	 5	10	1.5	0	0	0.11	0.13	0.16
0.01	-1.00	10	10	1.5	0	0	0.11	0.12	0.16
0.50	-1.00	10	10	1.5	0	0	0.10	0.11	0.15
1.00	-1.00	10	10	1.5	0	0	0.08	0.09	0.12
0.01	0.01	10	10	1.5	0	0	0.24	0.27	0.36
0.50	0.01	10	10	1.5	0	0	0.17	0.19	0.26
1.00	0.01	10	10	1.5	0	0	0.11	0.12	0.16
0.01	1.00	10	10	1.5	0	0	0.11	0.12	0.16
0.50	1.00	10	10	1.5	0	0	0.10	0.11	0.15

```

# Katie's comments: (7/3/2020)

**(Molly's responses in italics)**

Hi Molly, here are some thoughts:

(1) I suspected that the bootstrap and permutation approach would give different insights into the significance, so there may not be a problem, just different insights that we have to interpret and clearly communicate in the paper. 

*Do you have a good paper I can reference to better understand the differences in inference? I have read your notes and the section in the whitlock book but I'm not 100% on the differences here*

(2) My following comments are assuming your graphs are plotting the raw data results 

*yes, these are raw*

(3) For the Cov - I suspected the permutation approach was conservative, but the bootstrap was less conservative. If your 2x2 plot at the top is based on whether P<0.025, then what is plotted here makes sense to me. The permutation is not a conservative approach. Could you make a similar 2x2 plot, but base significance on whether or not the 95% CI overlap with 0?
 
 *The above approach was 95% CI. The Covariance plot below has been remade with values based on p <0.025. (oops)*

(4) Can you send me the line of code that you use to calculate the CI? Do you use the quantile function? There are different ways to do it, which have different assumptions. Can you confirm that you use 1000 for bootstraps and permutations?

*boot_df is the dataframe containing all the bootstrapped estimates. Right now I am only using 99 bootstraps because I wanted to quickly test out the function of the new GxE means pvalue approach as well as how well setting interaction = n_pop worked. I can re-run a subset at 1000 if you would like to see what that looks like*
```{r}
# Covariance Confidence Intervals 
cov_CI = quantile(boot_df_raw$covariance, probs=c(0.025, 0.975), type=1) 
cor_CI = quantile(boot_df_raw$cor_est_boot, probs=c(0.025, 0.975), type=1) 
cov_corrected_CI = quantile(boot_df_raw$cov_corrected_boot, probs=c(0.025, 0.975), type=1) 

# GxE Confidence Intervals
GxE_orig_CI = quantile(boot_df_raw$GxE_emm_original_boot, probs=c(0.025, 0.975), type=1) 
GxE_emm_CI = quantile(boot_df_raw$GxE_emm_boot, probs = c(0.025, 0.975), type=1)
GxE_omega_CI = quantile(boot_df_raw$GxE_omega_boot, probs=c(0.025, 0.975), type=1)
GxE_eta_CI = quantile(boot_df_raw$GxE_eta_boot, probs=c(0.025,0.975), type = 1)
GxE_SSq_CI = quantile(boot_df_raw$GxE_SSq_boot, probs = c(0.025,0.975), type = 1)
```
(5) Because cov(G,E) is bounded by -1 and 1, and GxE is bounded by 0 and 1, we are bound (pun intended) to see weird behavior near the boundaries

*ha*

(6) Can you confirm that you only implement the bootstrap once, and then use the same bootstrap for CovGE and GxE CI?

*I run through the following code x number of times - meaning that the same bootstrap is used for CovGE and GxE CIs* 
```{r}
###################################
#####  Bootstrap -- Raw Data  #####
###################################

for(i in 1:n_boot){
  
  # Shuffle Data
  shuffle_dat <- bootstrap_raw(model_df) # This function resamples data with replacement within each genotype and environment
  
  # Anova model fit & GxE estimates
  m2 <- mod.GxE(shuffle_dat) # Insert shuffled raw phenotype dataframe - calculates covGE and GxE
  
  # GxE Estimates
  cov_matrix_boot <- m2[[1]]
  GxE_emm_original_boot <- m2[[2]]
  GxE_emm_boot <- m2[[3]]
  GxE_loop_output_boot <- m2[[4]] # GxE output 
  omega2_boot <- m2[[5]]
  eta2_boot <- m2[[6]]
  GxE_SSq_boot <- m2[[7]] 

  # Covariance Estimates
  cov_est_boot = cov(cov_matrix_boot$G_means,cov_matrix_boot$E_means)
  cor_est_boot = cor(cov_matrix_boot$G_means,cov_matrix_boot$E_means)
  correction_raw_boot = max(sd(cov_matrix_boot$E_means),sd(cov_matrix_boot$G_means))
  cov_corrected_boot = round(cov(cov_matrix_boot$G_means, cov_matrix_boot$E_means)/(correction_raw_boot^2),2)
  
  # Bootstrap dataframe
  boot_dat_raw <- data.frame("covariance" = cov_est_boot,
                             "cor_est_boot" = cor_est_boot,
                             "cov_corrected_boot" = cov_corrected_boot,
                             "GxE_emm_original_boot" = GxE_emm_original_boot,
                             "GxE_emm_boot" = GxE_emm_boot,
                             "GxE_omega_boot" = omega2_boot,
                             "GxE_eta_boot" = eta2_boot,
                             "GxE_SSq_boot" = GxE_SSq_boot)
  boot_df_raw <- rbind(boot_df_raw,boot_dat_raw)
}
```
(7) I would need to see a more thorough evaluation, but your 3-plot of the GxE does not seem to match your plot of the CI overlapping with the true values. Can you check that you are calculating the True positive and False positive rates based on the true values, and not on the estimated values from the sims? Maybe you could add the true values to the 2x2 plots. (edited) 

*In the below plot, I've added true values with the estimated ones. True values are shown in green, estimated values are shown in red. I've ordered according to estimated covariance. I based my false negs and pos's on estimates. I actually do not bootstrap CIs around the true value - should I?*
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.3.CovCI_ordered.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.3.GxE_confint_ordered.png)


(8) Maybe you could also organize the simulation number in a way that makes the order easier to interpret (E.g, low to high true values of Cov(G,E)) 

*see above*
