# SIMULATION RESULTS! #

Here I present several plots created from the first full-scale simulation! In this simulation, I ran the following parameter set:

```{parameter set}
param_list <- list(
  reps = c(10),
  delta_env = c(0,0.25,0.5,0.75,1),
  delta_gen = c(-1,0,1),
  sample_size = c(5,10,20), 
  n_genotypes = c(2,4,8,16),
  std_dev= c(0.1,0.5),
  interaction= c(0,0.5,1)) 
  
n_bootstraps = 100
```
This generated a table with 10,800 rows, which was reduced down to 9360 after excluding rows in which both delta_env and delta_gen = 0. I submitted these 9360 in batches to Discovery and compiled the resultant data in a single .csv. 

## True covariance and True GxE
These panels show the relationship of "true" GxE and CovGE (i.e., no random noise added to phenotype). I have divided into panels according to the number of genotypes (2,4,8, or 16) and colored points according to sample_size (5,10,20).

First impressions: 
a) Number of genotypes has big impact on spread/coverage. 
b) 2 genotypes can reach +1 or -1 covariance, but multiple genotypes do not exceed covariance values of +/- 0.5
c) Tradeoff seems to be present - If range of cov = +/-1, then GxE does not exceed 1. But if cov does not exceed +/- 0.5, then GxE can extend above 1. (Im sure more can be said here, but this is just a first observation). 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_03152020/TrueCov_GxE.png)

## Covariance and GxE Estimates
Now showing the same plot layout but now using the data with random noise. These panels look similar to the above, which is a good initial indication that everything worked as it was supposed to. 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_03152020/GxE_cov_estimates.png)

## Significant Covariance and GxE Estimates
Same as above, but only showing those data that have significant covariance and gxe estimates (p<0.05). 

This plot has me thinking perhaps something went wrong: 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_03152020/SigGxE_cov.png)

Why are GxE values near 0 showing up as significant? A look at the data shows that these are estimates slightly off zero with very small confidence intervals but from a variety of different starting parameters. However, then there is the blank space as magnitude of GxE increases that are not statistically significant... why are these small estimates significant but not the larger ones?  

 delta_env |delta_gen |sample_size|n_genotypes | std_dev | interaction |true_cov |cov_estimate|cov_lwrCI|cov_uprCI|cov_pvalue|true_GxE |GxE_estimate| GxE_lwrCI|  GxE_uprCI| GxE_pvalue 
 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
  0.50 |1|5|8|0.5|0| 4.557143e-01| 0.43726556|0.414614449|0.458338031|0.000|2.664535e-15|0.025624170|0.0057347299|0.20497646 |0.04950495  
  0.75  |      -1      |    10       |    8   |  0.5    | 0| -5.477143e-01  |-0.52953110| -0.546005341 |-0.51518242| 0.00990099 |2.065015e-14  |0.007083631 |0.0012061898 |0.10145796 |0.00990099   
  1.00   |      1       |   10     |      8  |   0.5     | 0 | 5.705357e-01 |  0.55723254 | 0.544511283 | 0.57069669| 0.00000000 |2.575717e-14  |0.011259819| 0.0007435781 |0.07277978 |0.01980198
  0.00   |      1     |     20     |      8  |   0.5   |  0 | 1.031950e-15 |  0.01543271 | 0.004087529 | 0.02765174| 0.00000000 |3.552714e-15  |0.006618095 |0.0016012738 |0.10333969 |0.02970297
  0.75    |     1     |    20     |      8  |   0.5    |       0 | 5.481429e-01 |  0.53131347 | 0.520841860 | 0.54268961| 0.00000000 |8.881784e-15  |0.005441563 |0.0012904921| 0.08105676 |0.01980198
  0.25     |   -1      |     5    |      16 |    0.5   |        0| -2.507843e-01 | -0.24216729 |-0.249369251|-0.23666770 |0.00990099 |4.862777e-14  |0.005811667 |0.0022374175 |0.11026737 |0.01980198   
  
## 3.20.2020 Update
Per Katie's suggestions, I visualized boxplots of phenotype across environments for different parameter sets. 

I'm sorry this is just points - I could not for the life of me get geom_boxplot() to cooperate with me. This decently shows the patterns though... 

Plot 1 - This just shows that sample size doesn't drive much difference in spread, so I'm going to use sample size of 10 for plot 2
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_03152020/group%20%20vs%20sample%20size.png)

Plot 2 - Each columns are GxE, rows are interaction term. Covariance and GxE estimates and significance are shown. Obviously something is wrong with my permutation code. I will post permutation plots later, and I will also be going through my code to make sure there are no dumb errors. 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_03152020/deltaenv_interaction.png)


