# Test of means and bootstrapping approach

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

![image](https://github.com/RCN-ECS/CnGV/blob/master/img/GxE_test.png)

![image](https://github.com/RCN-ECS/CnGV/blob/master/img/GxE_test.png)

As is plainly observable, the means data are underestimating covariance and GxE magnitude, with CIs that do not overlap in many cases. 
I need to revisit this and see what is causing these discrepancies. Stay tuned.
