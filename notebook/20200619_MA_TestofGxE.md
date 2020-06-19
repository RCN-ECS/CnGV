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
###GxE vs. Covariance

*Note* that I am using GxE with new for-loop approach based on estimated marginal means. Covariance is the new, "corrected" covariance. 
Points colored according to significance: Blue = GxE Sig, Green = CovGE sig, Red = Both sig
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/619_covGxE.png)

Some points are stacked on eachother as evidenced by: 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/619_Hex.png)

### Covariance Means vs Raw
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/619_CovMeansRaw.png)

### Covariance Raw with Error vs Raw without Error
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/619_CovErrvsNoerror.png)

**Why are there some true values that equal 1?** 
This occurs when there is no interaction term (int =0) and delta_env and delta_gen equal eachother and are positive. In these cases, having no error produces perfectly parallel lines, equally spaced apart that will always have covGE = 1, even though it does not reflect the actual covariance produced by the data when there is variation around the mean. 

### GxE Means vs Raw
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/619_GxE_meanvsRaw.png)

### GxE Raw with error vs Raw without error
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/619_GxENEvsE.png)

There is discrepancy at the low GxE numbers. These are colored according to significance just to check how this discrepancy affected significance. Blue = GxE is sig, Green = Cov is sig, Red = both are sig. 

Similar as covariance above, this is an artifact of the error. In this case, delta_env and delta_gen do not necessarily equal eachother, but also have an interaction of 0. The GxE in these cases is driven by standard error (when error is gone, so goes the interaction). See example plots below. Top plot is phenotype plotted *with std deviation* and the bottom plot is the same data *without standard deviation*
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/619_GxEerror_phen.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/619_GxEnoerror_phen.png)
