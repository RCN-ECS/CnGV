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

Just as the preliminary sims showed, we see that the number of sample size affects ability to detect GxE (higher sample size = greater ability to detect GxE) while number of populations affects covariance (higher number of populations = greater ability to detect significant covariation). 

Note that the magnitude of GxE increases with n_pop because we set our interaction term to increase with n_pop. Data standardization  

#### GxE calculated using Estimated Marginal Means
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.7.gxe_covEMM.png)

#### GxE calculated using Omega^2
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.7cov_gxe_omega2.png)

#### HexPlot to show coverage
To ensure that all parameter space is being covered, here is a hexplot that quantifies the number of cases in hexagonal bins. 
There are a couple areas without data, but it does not seem like there are any gaps that would drive a bias or blind spot.
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.7.hexplot.png)

## Power Heatmaps
This is a standard heatmap showing how power changes according to standard deviation, sample size, and number of populations. These plots only show the power for covGE and GxE that fall between 0.4 and 0.6, since these windows encapsulate the full range of power.
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.11Gxe_Cov_heatmap.png)

## Phenotype Examples
Here are plots showing examples of the reaction norms look like in co/counter gradient scenarios with low and high GxE. Its interesting how our ability to visually identify these patterns (which is how these patterns have been typically been. identified in the past) is almost impossible as the number of populations/magnitude of GxE increases. *These all have standard deviation of 0.5*

| Number of Populations | |
|--- | --- |
| **2** | ![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/2pop.png)|
| **3** | ![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/3pop.png)|
| **5** | ![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5pop.png)|
| **10** | ![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/10pop.png)|
| **15** | ![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/15pop.png)| 

## Combines
Holding the total the same (either 50 or 100), I played around with the four different combinations of n_pop and sample_size that resulted in those 2 totals
| Sample Size | N_pop | Total |
| --- | --- | --- |
5|10|50
10|5|50
5|20|100
20|5|100

## Suggested Checks: 

#### True CovGE and GxE vs. Estimated Values
Katie suggested that I make sure my confidence intervals for covGE and GxE are overlapping with the known TRUE value (i.e., the covGE and GxE without any standard deviation. The plots below show **only those that have true values fall outside of the 95% confidence interval**. Each point is the value from a single replicate set of parameters. Red points are the true value, and the black lines are the 95% confidence intervals of the estimates. 

Across all replicates, 4392 out of 13500 (33%) have true GxE values that fall outside of the confidence intervals. In contrast - 1304/13500 (10%) true CovGE values are outside the confidence intervals. But if you look at the plots, the vast majority have confidence intervals that are very close to the true value, so I am not sure if this is acceptable or expected. There is no particular parameter driving these patterns. 

**Covariance**
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.8.covanoms.png)

**GxE**
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.8.GxEanoms.png)


#### Means vs raw

I took 2 approaches. First I checked the estimates themselves - made sure that CovGE and GxE for means and raw match. If so, they should fall along the 1:1 (red) line. Covariance matches perfectly, GxE is close but the means estimates seem to be a tad lower than raw GxE estimates. 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.8.covmeansraw.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.8.gxemeancheck.png)

I also regressed the lengths of the confidence intervals for the means against the length of the CIs for raw data. For this, I also colored the points according to significance to make sure there were no biases in significance. **Red** = both means and raw are significant, **blue** = means are sig but raw is not, **green** = raw is sig but means are not, and **grey** = neither are significant.Again, if they matched we would expect a perfect 1:1 fit. For CovGE, we see some variability - seems like most fall around the 1:1 line with perhaps a slightly greater error for means. 

For GxE, we see much greater variability - greater error for means data than raw, but apart from a few blues/greens near GxE = 0; there doesn't seem to be a major bias in significance. For CovGE, we again see more variability but a fairly even spread across the 1:1 line. It also looks like there is more significance for raw (green) but again it doesn't look like there's any systematic bias.

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.8.coverrormeans.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.8.GxEmeanserror.png)





