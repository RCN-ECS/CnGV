
I have gone back and forth on whether to use bootstrapping or permutation to assign significance to CovGE. So here are some plots to guide decisions. 

## False Positives 

Per Katie's instruction - False positives are subject to residual error and sampling design. 

Below I have plotted false positive rates for each metric (Covariance + bootstrap, Covariance + permutation, etc.) using the larger residual error (standard deviation = 1) according to sampling design. These are as expected. 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/12.18.FRT_FalsePos.png)

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/12.15.CG_FalsePositive.png)

# False Negatives -- Effect Size
Per Katie's instruction - False positives are subject to residual error, sampling design, and effect size. With that in mind, I've plotted Power (1-beta) as bar plots to visualize how power changes according to effect size (divided into 4 effect sizes), design, number of populations, and sample size. 

**Full Reciprocal Transplant**
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/12.18.FRT.GxE.Power.png)

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/12.18.FRT.Cov.Power.png)

**Common Garden Design**
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/12.18.CG.Power.GxE.png)

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/12.18.CG.Cov_Power.png)

# False Negatives -- Total Sample Size
As suggested by Katie, I've also plotted Power (1-beta) according total sample size. Same plots as above, just showing 128 total samples

**Full Reciprocal Transplant**
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/12.18.FRT_GxE_128samples.png)

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/12.18.FRT_Cov_128Samples.png)

**Common Garden Design**
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/12.18.CG.GxE.128Samples.png)

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/12.18.CG.Cov.128Samples.png)


**Conclusions**

I think using permutation for CovGE and GxE are both good for the paper


