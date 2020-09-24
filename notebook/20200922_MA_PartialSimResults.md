# Simulation Results

Last week I started the final (maybe, probably not) simulation. I settled on the following set of parameters: 

I capped the total samples (calculated as sample size * number of populations * number of environments) at 256 to reduce computational time.
I also eliminated all rows in which both delta_env and delta_gen == zero. 

This produced 38,080 total sims. I ran these on the short partition, submitting 1000 sims (each job separated by 5 seconds) every 8 hours. 

Average time per simulation: 9min 23sec

```{params}
param_list <- list( 
  reps = c(10), 
  delta_env = c(0.0,0.25,0.5,.75,1.0), 
  delta_gen = c(-1,0.0,1),
  sample_size = c(2,4,8,16), 
  n_pop = c(2,4,8,16),
  env_scenario = c(1,2),  # 1 = npop = n_env; 2 = multiple pops across 2 envs
  std_dev= c(.5, 1), # Scale
  interaction = NULL) # Set as function of n_pop
```
## Covariance by GxE (Estimated Marginal Means)

There are more datapoints in the higher n_pop because we set the length of the c(interaction terms) = n_pop. Thus, there are 8x more points in the n_pop = 8 than n_pop = 2; 4x more in the n_pop = 4 than n_pop = 2. We may consider randomly sampling so there are even numbers across n_pops?    

Legend: 
*Green* = Covariance is significant (defined by 95% CIs that don't overlap with zero)
*Blue* = GxE is significant (defined by p < 0.05)
*Red* = Both GxE and Covariance significant
*Grey* = Neither is significant

**N_env == N_gen**
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_CovGxE1.png)

**N_env = 2**
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_CovGxE2.png)

Something is going on with covariance here that is restricting range. Need to look into it.

## Covariance by GxE (Omega Squared)

**N_env == N_gen**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_CovGxE_omega1.png)

**N_env = 2**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_CovGxE_Omega2.png)

## Omega Squared vs. Estimated Marginal Means
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_AnovavsEmm.png)

## Power Analyses

**N_env == N_gen**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_PowerStdScale1.png)

**N_env = 2**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_PowerStandardized2.png)

## Barplot
Shows the proportion of datapoints that are each color (refers to legend above)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_PowerBarplot.png)

## Tradeoff between GxE and Covariance
To calculate this, I filtered the data to just those values that had a significant GxE (EMM or Omega, depending) or significant covariance (based on CIs).

**Estimated Marginal Means**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_TradeOffUpdated.png)

**Omega Squared**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_Tradeoff_Omega.png)

## Covariance: Estimate vs. True
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_CovRawVsTrue_Updated.png)

## Covariance: Means vs. Raw 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_CovMeansVSRaw.png)

## Covariance: Error (CI length) from Means vs. Raw
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_CovErrormeansRawComparison.png)

## GxE: Estimate vs. True
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_GxE_RawVsTrue_updated.png)

## GxE: Means vs. Raw 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_GxEMeansvsRaw.png)

## GxE: Error (CI length) from Means vs. Raw
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_GxE_MeansVsRawError.png)


## Confusion Matrices

| --- | Covariance Confidence Intervals | Covariance Bootstrap | GxE Confidence Intervals | GxE Bootstrap |
| --- | --- | --- | --- | --- |
| False Negative | 894 | 347 | 278 | 12 |
| False Positive | 1 | 10 | 0 | 234 |
| True Negative | 133 | 124 | 281 | 47 |
| True Positive | 92 | 639 | 561 | 827 |

Legend: 
*Red* = Estimated value
*Green* = True value
Ordered according to the estimate

**Covariance Bootstrap**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_CovBootConfusion.png)

**Covariance Permutation**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_CovPermConfusion.png)

**GxE Bootstrap**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_GxEBootConfusion.png)

**GxE Permutation**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_GxEPermConfusion.png)

## GxE P-Value Sanity Checks

**GxE Raw**

Legend: 
*Purple* = True GxE value is zero
*Orange* = True GxE value is not zero

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_GxEAnovaVsEmm_updated.png)

**GxE Means**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_GxEMeans_EmmvsAnova.png)

## Parameter Coverage

**N_env = N_pop**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_Hex1.png)

**N_env = 2**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_Hex2.png)
