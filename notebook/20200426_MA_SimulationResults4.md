# Full simulation Results


After resolving the covariance correction [issue](https://github.com/RCN-ECS/CnGV/blob/master/notebook/20200416_MA_SimResults_Round3.md), I ran ten replicates of the full array of parameters which were: 

```{param list}
param_list <- list( 
  reps = c(10), # or more?
  delta_env = c(0.01,0.5,1),
  delta_gen = c(-1,0.01,1),
  sample_size = c(5,10,20), 
  n_pop = c(2,3,5,10,15), 
  n_environments = NULL,
  std_dev= c(0.5,1), 
  interaction = 5) # Vector LENGTH not magnitude - Interaction levels = 5 values from 0 to n_pop.
```
Note that I changed delta_env and delta_gen from 0 to 0.01. I did this to improve representation of very small changes in slope and intercept (before, zeros would result in a zero for phenotype). 

I've been playing around with plots and here's what I've come up with so far. I'd like a better way of showing the tradeoff between sample size and n_pop (see first plot below) but I can't come up with a good way of showing it. Any suggestions would be welcome.

## Covariance and GxE
Here is an updated plot showing estimates of GxE and covariance, colored according to significance. Rows are number of populations/genotypes, while columns are sample sizes. 
**Blue** = GxE is significant
**Green** = Covariance is significant
**Red** = Both are significant
**Grey** = None are significant

Just as the preliminary sims showed, we see that the number of populations/genotypes affects ability to detect GxE (higher pops = greater ability to detect GxE) while sample size affects covariance (higher sample size = greater ability to detect significant covariation). I'm still mulling over why this is (I feel like I'm on the cusp of a lightbulb moment! My lightbulbs take a while to turn on sometimes.)

Note that the magnitude of GxE increases with n_pop because we set our interaction term to increase with n_pop. 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.5.covgxe.png)

## Power Heatmaps
This is a standard heatmap showing how power changes according to standard deviation, sample size, and number of populations. 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.6.PowerPlots.png)

## Phenotype Examples
Here are plots showing examples of the reaction norms look like in co/counter gradient scenarios with low and high GxE. Its interesting how our ability to visually identify these patterns (which is how these patterns have been typically been. identified in the past) is almost impossible as the number of populations/magnitude of GxE increases. 

| Number of Populations | |
|--- | --- |
| 2| ![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/2pop.png)|
| 3| ![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/3pop.png)|
| 5| ![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5pop.png)|
| 10| ![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/10pop.png)|
| 15| ![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/15pop.png)| 

Any suggestions for plots/improvements are welcome. 
