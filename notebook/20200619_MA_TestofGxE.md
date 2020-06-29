# Check New GxE approach and functionality of functions
To test the new GxE approach and broadly check if my new "functioned" code is working properly, I ran the following limited set of parameters on the cluster: 

```{r}
param_list <- list( 
  reps = c(5), 
  delta_env = c(0.01,0.5,1),
  delta_gen = c(-1,0.01,1),
  sample_size = c(5,10),
  n_pop = c(1), # Number of pops per environment (New approach allows # to change with greater flexibility)
  n_environments = c(2,5,10),
  std_dev= c(0.5,1.5),
  interaction = 2) # Vector LENGTH not magnitude

```
### GxE vs. Covariance

*Note* that I am using GxE with new for-loop approach based on estimated marginal means. Covariance is the new, "corrected" covariance. 
Points colored according to significance: Blue = GxE Sig, Green = CovGE sig, Red = Both sig
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/619_covGxE.png)

Some points are stacked on eachother as evidenced by: 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/619_Hex.png)

### Covariance Means vs Raw
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/619_CovMeansRaw.png)

### Covariance Raw with Error vs Raw without Error
Colored according to significance: Green = CovGE sig, red = both sig, blue = GxE sig
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/619_CovRawvsTrue.png)

**Why are there some true values that equal 1?** 

This occurs in cases where there is no interaction term (int =0) and delta_env and delta_gen equal eachother and are positive. Having no error and no interaction produces perfectly parallel lines that always have covGE = 1. 
## Phenotype Example
This top plot has Covariance of -.79 and a pvalue of 0.06. The std. deviation is 1.5. 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/619_CovSTuff.png)

This plot is what those lines would look like if std deviation = 0. Note that now the covariance is listed as 1. 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/619_CovNoerr.png)


### GxE Means vs Raw
Looks like that observation that some commit a Type II error holds up. Doesn't look too bad though. 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/619_GxE_meanvsRaw.png)

### GxE Raw with error vs Raw without error
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/619_GxENEvsE.png)

There is discrepancy at the low GxE numbers. These are colored according to significance just to check how this discrepancy affected significance. Blue = GxE is sig, Green = Cov is sig, Red = both are sig. 

Similar as covariance above, this is an artifact of the std. deviation. In this GxE case, however, delta_env and delta_gen do not necessarily equal eachother, but also have an interaction of 0. The GxE in these cases is driven by standard deviation (when error is gone, so goes the interaction). See example plots below. Top plot is phenotype plotted *with std deviation* and the bottom plot is the same data *without standard deviation*
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/619_GxEerror_phen.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/619_GxEnoerror_phen.png)

### Raw GxE pvalue vs. Anova Pvalue 
Red lines drawn at 0.05 levels. All instances for raw and means data in which the anova p-value is < 0.05 but permutation p-value is > 0.05 occur in n_pop = 2, interaction = 2, standard deviation = 1.5 (higher of 2 levels) scenarios. 

There are 55 inconsistant pvalues (out of 1080) for the raw dataset, and 100 inconsistencies for the means dataset. I think these discrepancies are fairly mild. As shown in the second plot below, the majority of GxEmag values whos pvalues don't match are low (highest value with a pvalue discrepancy for raw data is 0.15).

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/619_GxE_pvalue.png)

### Problematic subset

Here I have pulled out those 55 datapoints that have an Anova pvalue of < 0.05, but Permutation pvalue > 0.05 and colored by the magnitude of GxE (using EMMs). There are multiple points stacked on top of one another, which is why there doesn't appear to be 55 points. 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/6.29.ProblemGxEs.png)

### Means GxE pvalues vs. Anova Pvalue

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/619_GxEMeanspvalu_anova.png)
