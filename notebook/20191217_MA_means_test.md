# Test of means and bootstrapping approach

## Update for 2/27/2020
Since writing this post, we have changed the way we generate data from using linear models to regression design matrices. As described below, we found the previously estimated means/SE data downwardly biased covariance and GxE estimates. 

I rewrote the code for the means/SE data in which I use the new data generation approach. In my new code, kept the overall approach the same (i.e., estimating covariance and GxE on the provided means without using Anovas), but I built in a small function in which I calculate the true covariance and GxE values (i.e., the values calculated on the full raw dataset with no noise). In doing so, I can compare how accurately the means/SE data are performing. 
**If CovGE and GxE estimates are accurate, they should match the TRUE values**

I ran a small simulation with the following starting parameters in which I varied only the standard deviation (but with a small interaction built in): 
```
param_list <- list(
  reps = 5,
  delta_env = c(1), # the amount the phenotype changes across 1 value of the environment (i.e., the slope). This is essentially the amount/degree of phenotypic plasticity that is the same across genotypes.
  delta_gen = c(-1), # the amount the phenotype changes from one genotype to the next. This is essitially the increase intercept from one genotype to the next.
  sample_size = c(5), 
  n_genotypes = c(3),
  n_environments = NULL,
  std_dev= c(0.01,0.5,1), # Random noise, with standard deviation of 1,
  interaction= c(0.5)) # this sd determines the amount of GxE)
```
I ran the full code which *includes standardizing the data*

Upon running this code I plotted the estimated CovGE and GxE against the true values. *If they were generated correctly, they should fall along a 1:1 line, with higher standard deviation being further from that line*
#### AND THATS WHAT I FOUND 
Looks like whatever bias that was generated in the old coding approach is reduced in the new coding approach! Wahoo! However, when standard deviation IS present, the covGE estimates do seem to be biased downwardly still. GxE seem fine.  

**Covariance:**
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/CovTest_means.png)

**GxE**
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/GxEtest_means.png)

### Why is CovGE increasing if the only thing changing is standard deviation? 
My guess was that the moderate (0.5) interaction term was leading to some amount of covariance being picked up. When I re-ran the above with interaction term of Zero, I got the expected pattern of no change in "true" covGE but just variation of standard deviation: 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/Cov_noInt.png)

### Error bars! 
After finding and fixing the bug in my bootstrap code, I ran a sim with sample sizes of 5 and 100. As expected, error gets smaller with the larger sample size. So that's a bit reassuring. Unfortunately, the error bars for the higher std. dev covGEs do not overlap with the 1:1 line, which we would hope for if there was no real bias. So it looks like for covariance, the bias of lower CovGEs with means/SE data persists. GxE all seem to overlap the 1:1 line though, suggesting those data are a bit more robust. 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/Cov_SampleSize_error.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/GxE_SampleSize_error.png)


## 12/17/2020: 
As noted [here](https://github.com/RCN-ECS/CnGV/blob/master/notebook/20191204_Categorical_Analyses.md), I recently re-made code into a more manageable unit.

In doing so, I modified the bootstrapping code for data means. (Bootstrapping approach for the raw data is the same as it was before).
In the new approach, I resampled means using the command below and assuming Central Limit Theorum.

```#Bootstrap means
# Shuffle data 
new_phen <- rnorm(nrow(cond_E), mean = cond_E$phen_data, sd =  cond_E$phen_mean_SE)      
```
To test whether this reliably estimated the same GxE and covariance estimates (and CI's) as the raw data, I ran the raw data through then `sim_means_se` function that generates means and SE from the simulated raw data.

I simulated multiple scenarios (different slopes, intercepts, errors, etc.) with each scenario identified by its "index". 

I ran both raw and means data through the `Categorical_sim` function which is set up to handle both data types. I concatenated the results from each and created plots to compare whether my means function generates CI's that are the same (or close) to raw data.

Below are the plots showing how the different approaches compare for Covariance and GxE.

![image](https://github.com/RCN-ECS/CnGV/blob/master/img/Covariance_test.png)

![image](https://github.com/RCN-ECS/CnGV/blob/master/img/GxE_test.png)

As is plainly observable, the means data are underestimating covariance and GxE magnitude, with CIs that do not overlap in many cases. 
HOWEVER -- the CIs estimated by both data types nicely match (to the naked eye anyway). This is good. 

I think the low-ball estimates are because the GxE and Covariances for the raw data are based on model fits... while the means/SE are used as/is (no model fit). I will try a couple different approaches to see if any of them more reliably match raw data.

## Update Dec. 18,2019

To see if its because of the model fits, I ran the means through some code that would generate raw data based on the means and standard error. I then ran that raw data through the same `categorical_sim`. I overlaid those points (called new_from_means) below.  The new raw data does a better job matching the predicted covariance and GxE, as it should. **However, the confidence intervals do not match the raw data**
![image](https://github.com/RCN-ECS/CnGV/blob/master/img/Cov2.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/img/GxE2.png)

This calls into question whether using means/SE is a good method. I'll continue to think.

## Update Dec. 20, 2019
Katie suggested that our method of scaling might be driving the weird patterns above. So I removed the scaling (phenotype-average/std.deviation) and re-ran all the code. 

It does seem like the current scaling strongly affected the means output. When you remove the scaling, the means/raw data all match nicely! If standardizing has such a large effect between data types, it defeats some of the purpose of standardizing in the first place, which was to be able to generalize across studies.

![cov wrong error](https://github.com/RCN-ECS/CnGV/blob/master/img/Means_raw_covariance_nocorrection.png)
![gxe wrong error](https://github.com/RCN-ECS/CnGV/blob/master/img/Gxe_meansandraw_nocorrection.png)

### How I set up Scaling

For raw data, I calculate the average and standard deviation of the full dataset. I then apply the following to each datapoint: 
```#data
# Standardize raw data
  dat_avg <- mean(input_df$phen_data) 
  dat_std <- sd(input_df$phen_data)
  input_df$phen_corrected <- ((input_df$phen_data - dat_avg)/dat_std)
```

For means/SE data, I similarly take the average and standard deviations of the means. I then apply the following to each mean datapoint: 
```#mean
# Standardize means data
  dat_avg <- mean(input_df$phen_data) 
  dat_std <- sd(input_df$phen_data)
  input_df$phen_corrected <- ((input_df$phen_data - dat_avg)/dat_std)
  
 ```
