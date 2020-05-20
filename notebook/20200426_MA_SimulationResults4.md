# Full simulation Results

After resolving the covariance correction [issue](https://github.com/RCN-ECS/CnGV/blob/master/notebook/20200416_MA_SimResults_Round3.md), I ran **100 replicates** of the full array of parameters which were: 

```{param list}
param_list <- list( 
  reps = c(10),
  delta_env = c(0.01,0.5,1),
  delta_gen = c(-1,0.01,1),
  sample_size = c(5,10,20), 
  n_pop = c(2,3,5,10,15), 
  n_environments = NULL,
  std_dev= c(0.5,1), 
  interaction = 5) # Vector LENGTH not magnitude - Interaction levels = 5 values from 0 to n_pop.
```
Note that I changed delta_env and delta_gen from 0 to 0.01. I did this to improve representation of very small changes in slope and intercept (before, zeros would result in a zero for phenotype). 

I've been playing around with plots and here's what I've come up with so far (as of May 6, updated May 11). I'd like a better way of showing the tradeoff between sample size and n_pop (see first plot below) but I can't come up with a good way of showing it. Any suggestions would be welcome.

## Covariance and GxE
Here is an updated plot (5.6) showing estimates of GxE and covariance, colored according to significance. Columns are number of populations/genotypes, while rows are sample sizes. **THESE WERE INCORRECTLY LABELED BEFORE**

|Color|Meaning|
|---|---|
**Blue** | GxE is significant
**Green** | Covariance is significant
**Red** | Both are significant
**Grey** | None are significant

Just as the preliminary sims showed, we see that the number of populations affects ability to detect GxE (higher n_pop = greater ability to detect GxE) while sample size affects covariance (higher number of samples = greater ability to detect significant covariation). 

*Note that the magnitude of GxE increases with n_pop because we set our interaction term to increase with n_pop.* 

#### GxE calculated using Estimated Marginal Means

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.19.covandGxe.png)

#### GxE calculated using Omega^2

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.19_GxEcov_omega2.png)

#### Tradeoff with GxE and Covariance

We predicted that since as GxE increases, the relationship between genotype and environment must become increasingly independent which should reduce the amount of CovGE observed. Looks like our hypothesis is supported! *This data has been filtered to just significant GxE OR significant CovGE values.*

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.19.CovGEtradeoff.png)


## Different way to show tradeoff betweeen GxE and covariance: 
These show covGE and GxE between 0.4 and 0.6, since these windows encapsulate the full range of power.

|GxE | Cov | Both |
|---|---|---|
|![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.19.bluedat.png)|![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.19.greendat.png)|![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.19.reddat.png)|


#### 2 population claim by Conover/Schultz
Conover and Schultz claimed that 2 populations was insufficient to estimate covariance. Well, it turns out they weren't entirely wrong. Here is the power for n_pop = 2:

| Std. Deviation | Sample Size | Power | 
|---|---|---|
 0.5 | 5 | 0.10
0.5| 10 | 0.12
0.5| 20 | 0.10
1| 5 | 0.12
1| 10 | 0.12
1| 20 | 0.12

#### HexPlot to show coverage

To ensure that all parameter space is being covered, here is a hexplot that quantifies the number of cases in hexagonal bins. After running the full 100 replicates, there don't seem to be any major gaps or areas of oversampling that would drive a bias or blind spot.

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.19.hexplot.png)

## Power Heatmaps
This is a heatmap showing how power changes according to standard deviation, sample size, and number of populations. 

**Full Parameter Coverage**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.19.FulldataHeatmap.png)

**Reduced Coverage -- Only covGE and GxE between 0.4 and 0.6, since these windows encapsulate the full range of power.** 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.12.PowerHeatmap.png)

## Phenotype Examples

Examples of the reaction norms in co/counter gradient scenarios with low and high GxE. It is interesting how our ability to visually identify these patterns (which is how these patterns have been typically been identified in the past) is almost impossible as the number of populations/magnitude of GxE increases. *These all have standard deviation of 0.5*

| Number of Populations | |
|--- | --- |
| **2** | ![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/2pop.png)|
| **3** | ![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/3pop.png)|
| **5** | ![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5pop.png)|
| **10** | ![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/10pop.png)|
| **15** | ![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/15pop.png)| 

## Combinations

Holding the total the same (either 50 or 100), I played around with a couple different combinations of n_pop and sample_size that resulted in those 2 totals (see table). Here are a few plots. I think my favorite is the heatmap for the Total = 50, because it shows how the n_pop and sample size differently affect power. 

| Sample Size | N_pop | Total |
| --- | --- | --- |
5|10|50
10|5|50
10|10|100
20|5|100

#### Covariance vs. GxE

**N = 100**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.11.hundie_covgxe.png)

**N = 50**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.11.fitty_covgxe.png)

#### Heatmaps

Note the legend's power scale is a bit misleading here. I'd have to change the scale if we included it in the paper.

**N = 100**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.12.HundieHeatmap.png)

**N = 50**
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.12.FittyHeatmap.png)


## Suggested Checks: 

#### True CovGE and GxE vs. Estimated Values
Katie suggested that I make sure my confidence intervals for covGE and GxE are overlapping with the known TRUE value (i.e., the covGE and GxE without any standard deviation. **The plots below show only those that have true values fall outside of the 95% confidence interval**. Each point is the value produced by one set of parameters. Red points are the true value, and the black lines are the 95% confidence intervals of the estimates. 

Across all 10 replicates, 4392 out of 13500 (33%) have true GxE values that fall outside of the confidence intervals. In contrast - 1304/13500 (10%) of true CovGE values are outside the confidence intervals. If you look at the plots, the vast majority have confidence intervals that are *close* to the true value, so I am not sure if this is acceptable or expected. A brief visual check of the parameters suggests there is no particular parameter driving these patterns. 

**Covariance**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.8.covanoms.png)

**GxE**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.8.GxEanoms.png)


#### Means vs Raw data Check

To check whether the estimates based on means data match estimates from raw data, I took 2 approaches. First I checked the estimates themselves - made sure that CovGE and GxE for means and raw match. If so, they should fall along the 1:1 (red) line. Covariance matches perfectly, GxE is close but the means estimates seem to fall tad lower than raw GxE estimates. 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.8.covmeansraw.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.8.gxemeancheck.png)

I also regressed the *lengths of the confidence intervals* for the means against the *lengths of the CIs* for raw data. For this plot, I also colored the points according to significance to make sure there were no biases in significance. Again, if the lengths of CIs matched we would expect a perfect 1:1 fit.

Color | Meaning
|---|---|
**Red** |estimates from means and raw are significant
**Blue** | estimates from means data are significant but estimates from raw are not signficant
**Green** | raw is sig but means are not significant
**Grey** | neither are significant.  

For CovGE, we see some variability - seems like most fall around the 1:1 line with perhaps a slightly larger error for means than raw (more points cluster above the 1:1 line. It looks like there is more significance for raw data (green) but again it doesn't look like there's any systematic bias.

For GxE, we see much greater variability - larger error for means data than raw data, but apart from a few blues/greens clustering around GxE = 0, there doesn't seem to be a major bias in significance. 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.8.coverrormeans.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.8.GxEmeanserror.png)

Please let me know of any thoughts or additional checks I can make. 

# Analysis on Empirical Data: 

For this I used the raw data from my frog study and I extracted the means from Geoff's 2002 MEPS paper. I just picked one phenotype per analysis. For the raw data (molly's data), I chose the age at metamorphosis data (only using freshwater and 6ppt treatments as "freshwater" and "saltwater" for simplicity). For the means data (Geoff's data), I chose the Shell mass growth data (figure 3). 

I ran each through the framework and here's what I got: Very interesting that even though we both report CnGV -- its trending that way but not significant! 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.14.GeoffPlot.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.14.Mollyplot.png)
