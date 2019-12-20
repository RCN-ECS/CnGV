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
