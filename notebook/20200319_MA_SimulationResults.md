# SIMULATION RESULTS! #

## Updated for 3.30.2020
After some tinkering this morning, we came up with the plot with points colored in terms of significance: 
Red = Both GxE and covariance are significant
Blue = Just GxE is significant
Green = Just covariance is significant
Grey = Neither are significant. 

That produced this plot: 
![image]()
## Updated for 3.27.2020
Below I present the results of a sim that was wrong --- see [here](https://github.com/RCN-ECS/CnGV/blob/master/notebook/20200325_MA_SimulationDiagnostics.md) for details.

Now that its correct, I ran a new set of parameters - full spectrum but reduced representation to get a feel for how it worked 
```{new params}
# Starting list of parameters
param_list <- list( 
  reps = c(10),
  delta_env = c(0,0.5,1),
  delta_gen = c(-1,0,1),
  sample_size = c(2,4,8,16), 
  n_pop = c(2,4,8,16),
  n_environments = NULL,
  std_dev= c(0.1,0.5),
  interaction= c(0,1))

```
## True covariance and True GxE
These panels show the relationship of "true" GxE and CovGE (i.e., no random noise added to phenotype). Columns are the number of pops (2,4,8, or 16) and rows are sample size (also 2,4,8,16).
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_03152020/trutrue3.27.png)

## Covariance and GxE Estimates
Same plot layout but now using the data with random noise.
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_03152020/estest3.27.png)

## Significant Covariance and GxE Estimates
Same as above, but only showing those data that have significant covariance and gxe estimates (p<0.05). 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_03152020/sigsig3.27.png)

## Substitute Runif
New parameter set for limited test: 
```{new set}
param_list <- list( 
  reps = c(1),
  delta_env = runif(5,0,1),
  delta_gen = c(-1,0,1),
  sample_size = c(4), 
  n_pop = c(16),
  n_environments = NULL,
  std_dev= c(0.1),
  interaction= c(0.1))
```
## 16 Population Comparison of Two Approaches
Using the original approach where delta_env = c(0,0.5,1) vs. delta_env = runif(5,0,1).

Top plot only has 16 datapoints so I could get it done semi-quickly. The bottom plot is from the full simulation so there are many more datapoints, but I filtered the full dataset to only the datapoints with sample size = 4, n_pop = 16, std. dev = 0.1 and int = 0.0.
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_03152020/Runifvsfull.png)





### Results from earlier which we now know are WRONG.

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

## 3.24.2020 Update

I changed the code so that we run covariance on 2 genotypes (populations), but a correlation on 3+ genotypes. This is to account for the weirdness that happens to covariance when we standardize the data. I also checked the code to make sure there were no silly errors anywhere. I couldn't find any, so hopefully there are no hidden bugs. I then ran a quick round of sims on the extreme values (eg- low and high sample size, low and high n_pop (which used to be n_genotype)

```{values}
# Starting list of parameters
param_list <- list( 
  reps = c(10),
  delta_env = c(0,0.5,1),
  delta_gen = c(-1,0,1),
  sample_size = c(2,16), 
  n_pop = c(2,16),
  n_environments = NULL,
  std_dev= c(0.1,0.5),
  interaction= c(0,1))
```

Here is one of the results: 

This plot shows JUST significant covariance and GxE values. Columns are the number of pops, rows are sample size. 

This makes more sense in terms of the tradeoff between GxE, Covariance, population number and sample size. THe max number of signifiant values are found with 2 populations but 16 replicates. this makes sense. 

I'm still not sure why there is a cluster of data points near -1 and 1 covariance.... need to investigate that further. 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_03152020/324sim_sig.png)



