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

### Graph note for below: change color scheme to work in B&W, use different shapes for each point, add legend, add text to column-panels (number of populations?) and row-panels (sample size?). I'm a little concerned that the 2 pop with 2 sample size seems to have so much power, I would think that case would have no power? 

### KEL note: I've been thinking more about what to call our "GxE estimate" - I think we need to be more specific in graphs and the main text. I was thinking something like "GxE estimate (Deviation)" and in the figure legends citing the equation. We also need to decide how to talk about genotypes/populations/locations in the main text, and be consistent with the way these are graphed here. We also need to distinguish between the true measure and the estimate, for example "GxE estimate (deviation)" vs. "GxE true value (deviation)". Then we also need a term to distinguish the raw data from the means data.... "GxE estimate from means (deviation)" (?)

### KEL note: We also need to be consistent with what to call "Covariance" across the plots. Sometimes it is called CovGE. We should be consistent with the main text, for example "$Cov_{GE}$" or "Cov(G,E)". The former can be made in R, but it can be a pain. We also need to distinguish between the true measure and the estimate, for example "Cov(G,E) estimate" vs. "Cov(G,E) true value". Then we also need a term to distinguish the raw data from the means data.... Also minor but be consistent with capitalization (e.g. don't switch some graphs to "Cov(G,E) Estimate").

### KEL note: It would be helpful to start to order the graphs as the way you want to present them in the results, and indicate which you want to include in the main text vs. the Supp Mat. For example, do we want to start with the comparison of bootstrap vs. permutation with the confusion matrices?

**N_env == N_gen**
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_CovGxE1.png)

### Graph note for below: let's be sure to add the nenv=2 x npop=2 case to this plot, although it is redundant with the last plot it serves as an anchor for comparison. It's a little weird to me that the 4 population case seems to have lower power than the 2 population case in the last plot. Something weird might be going on with the calculation, since this would be a higher overall sample size. It would be good to double check the calculations for this case.

**N_env = 2**
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_CovGxE2.png)

Something is going on with covariance here that is restricting range. Need to look into it.

## Covariance by GxE (Omega Squared)

### Graph note for below: same comments as previous graphs

**N_env == N_gen**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_CovGxE_omega1.png)

**N_env = 2**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_CovGxE_Omega2.png)

### KEL Graph note for above graphs: I'm not sure if the above graphs should be in the main manuscript, as they would take up a lot of space. We might want to pick one or two comparisons to illustrate; I would prefer to see a comparison of different designs with same total sample sizes. Not a strong opinion yet, just food for thought.

## Omega Squared vs. Estimated Marginal Means

### Graph note for below: We need to come to consensus on terminology. For example if we use "GxE estimate (Deviation)", this should be for the estimate from the data sample. Then we could use "GxE (Deviation)" for the actual GxE in absence of environmental noise (and clarify this in figure legend). 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_AnovavsEmm.png)

## Power Analyses

### Graph note for below heatmaps: What amount of covariance and GxE is this calculation for (e.g. > 0.5?)? Power needs to be calculated as the power to detect a certain range of values. Can we incorporate that information into the main text of the plot? I assume these are for the raw data. Should we make equivalent plots for the means data? The confusion matrices may be enough.
### Will we be filling in the missing cells? Again for the 2_env case, I think we should add the 2 genotypes because although it is redundant, it would serve as an anchor for comparison. Overall, it's easy to see how increasing n_pops and sample size affects power, but it's hard to see what the best design is. For that reason I would lean toward these as supp graphs and the comparison of equal total sample size, for the larger SD, for the main paper.

**N_env == N_gen**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_PowerStdScale1.png)

**N_env = 2**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_PowerStandardized2.png)

## Barplot

### KEL note: I don't get much from this graph, since as you show at the bottom, the parameter space is not evenly sampled in any one case. So I think this is driven more by the way we set up the sims, but maybe you just wanted to see the proportions? I wouldn't advocate for including in the paper, but open to other thoughts.

Shows the proportion of datapoints that are each color (refers to legend above)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_PowerBarplot.png)

## Tradeoff between GxE and Covariance

### KEL note: The y-axis throws me off, maybe I would just need to read the methods and results for this. Make x-axis consisent with whatever we decide "True GxE (Deviation)" should be called and y-axis with "True $Cov_{GE}$"
### KEL note: Another way we could visualize is this to just plot "True GxE (Deviation)" vs. "True $Cov_{GE}$". Then each dataset could be a datapoint and it would allow us to see the variation more. I think a figure in this category should be in the main text.

To calculate this, I filtered the data to just those values that had a significant GxE (EMM or Omega, depending) or significant covariance (based on CIs).

**Estimated Marginal Means**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_TradeOffUpdated.png)

**Omega Squared**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_Tradeoff_Omega.png)

## Covariance: Estimate vs. True

### Note for below two graphs: Can't see the false positives from the true negatives. Use different points so works in B&W. Make axes labels consistent with what we decide for Cov(G,E). KEL likes the idea of making this a 2-panel plot for the main paper. Need to clarify what method was used to call true positives and false positives in simulation.

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_CovRawVsTrue_Updated.png)

## GxE: Estimate vs. True
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_GxE_RawVsTrue_updated.png)

## Covariance: Means vs. Raw 

### Notes for below 4 figures:  Make axes labels consistent with what we decide for GxE and Cov(G,E). Clarify if true value or estimates. KEL likes the idea of making the following 4 graphs a 4-panel plot for the Supp Mat.

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_CovMeansVSRaw.png)

## Covariance: Error (CI length) from Means vs. Raw
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_CovErrormeansRawComparison.png)


## GxE: Means vs. Raw 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_GxEMeansvsRaw.png)

## GxE: Error (CI length) from Means vs. Raw
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/923_GxE_MeansVsRawError.png)


## Confusion Matrices

### KEL thoughts on presenting confusion matrix: A table like this is confusing to present with numbers, but we should report the numbers in the Supp Mat. Also I think you have a typo, since the CI and Boostrap are the same method. Could you make an additional table, but convert the numbers to percentages so that each column adds up to 100%? In addition, we'll want to calculate the false positive (FPR) and false negative _rates_ for each method. FPR = FP / (number of Cov(G,E) that are 0, or FP + TN). FNR = FN / (number of Cov(G,E) that don't equal 0, or TP + FN). 

| --- | Covariance Confidence Intervals | Covariance Bootstrap | GxE Confidence Intervals | GxE Bootstrap |
| --- | --- | --- | --- | --- |
| False Negative | 894 | 347 | 278 | 12 |
| False Positive | 1 | 10 | 0 | 234 |
| True Negative | 133 | 124 | 281 | 47 |
| True Positive | 92 | 639 | 561 | 827 |

| --- | Covariance Confidence Intervals | Covariance Bootstrap | GxE Confidence Intervals | GxE Bootstrap |
| --- | --- | --- | --- | --- |
| False Negative | % | % | % | % |
| False Positive | % | % | % | % |
| True Negative | % | % | % | % |
| True Positive | % | % | % | % |

| --- | Covariance Confidence Intervals | Covariance Bootstrap | GxE Confidence Intervals | GxE Bootstrap |
| --- | --- | --- | --- | --- |
| False Negative Rate | % | % | % | % |
| False Positive Rate | % | % | % | % |

Legend: 
*Red* = Estimated value
*Green* = True value
Ordered according to the estimate

### Notes on below graphs: Will we include them in Supp Mat? KEL kind of likes them. If we do, let's make the error bars more transparent, add a legend, and use different points for true value vs. estimate.

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
