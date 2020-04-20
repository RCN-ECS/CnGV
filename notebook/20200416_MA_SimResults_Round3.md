# Simulation Results Round 3

In a previous notebook [post](https://github.com/RCN-ECS/CnGV/blob/master/notebook/20200408_KEL_MultiplePop_CnGv_GxE_UPDATE2.Rmd), Katie identified that we are getting biased covariance estimates that went down as n_pop went up because scaling the data was disproportionately affecting the ratio of the standard deviation of Emean and Gmean estimates. Fewer populations caused less of a difference in standard deviation, which is why it had greater effect as n_pops increased. 

If the standard deviation of Gmeans is greater, there is a greater difference in intercept. 
If the standard deviation of Emeans is greater, there is a greater difference in slope. 

[Standardizing the covariance by the squared standard deviation](https://github.com/RCN-ECS/CnGV/blob/master/notebook/20200413%20Screen%20Shot%202020-04-09%20at%201.23.04%20PM.png) of either Gmeans or Emeans (whichever is greater) preserves the correct relationship between Gmeans and Emeans which should allow for the proper estimation of covariance.    

I added this correction in and simulated across a small set of parameters to see how the new covariance performed: 
```{param list}

param_list <- list( 
  reps = c(10),
  delta_env = c(0,1),
  delta_gen = c(-1,1),
  sample_size = c(2,5), 
  n_pop = c(5,10), 
  n_environments = NULL,
  std_dev= c(0.25),
  interaction = 3) # Vector LENGTH not magnitude - this is a change to allow for magnitude of interaction term to vary

```
First, 4 vignettes to prove that everything is working ok: Note that these covariance estimates are corrected.

| --- | Negative Covariance | Positive Covariance | 
| --- | --- | --- | 
| Low GxE| ![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/plot_row5.png) | ![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/plot_row7.png) |
| High GxE|![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/plot_row61.png) | ![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/plot_row59.png)|

Now to look at how significance lines up according to GxE and Covariance. GxE is based on estimated marginal means; Covariance is corrected. Rows are sample size; columns are n_pop.

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/CorrectedCov.png)

Using the new Correction looks really good. I know other diagnostic plots may be useful, but I am unsure which ones to create. 
**What diagnostics would be helpful?*


**Data to Choose from**: In addition to the endpoints listed below, I have the following for each set of parameters, I have a) phenotype data, b) Gmeans and Emeans, c) bootstrap data, d) permutation data. 

I generated 6 different covariance metrics for raw simulated data and again for the means data (with CIs and p-values): 
1. True covariance uncorrected
2. Covariance estimate uncorrected
3. True correlation 
4. Correlation estimated
5. True covariance corrected
6. Covariance estimate corrected

Finally I have different GxE metrics (with CIs and p-values): 
1. True GxE using EMMs
2. GxE using EMMs
3. True GxE using Omega^2
4. GxE using Omega^2



## Means vs. Raw: 

Good news - the new correction eliminates any differences in covariance estimates between means and raw data (diff = 0) and the differences between true GxE on raw vs means are tiny too (see histogram of differences below): 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/hist_gxe%20diffs.png)

**Why are there zero differences in covariance but subtle differences in GxE?**

The means and EMMs are in fact different by very small margins. When I take the unstandardized covariance, these differences are preserved. However when I standardize the covariance, these differences are lost.

| EMM - Emeans | EMM - Gmeans | raw means - Emeans | raw means - Gmeans |
|---|---|---|---|
-0.15399250  |0.58823449 |-0.1517924 |0.57983053
0.09169482 |-0.20917987|0.0903848 |-0.20619137
-0.23684704 |-0.0377928|-0.2334633 | -0.03725295 
-0.14327593  | 0.19146278|-0.1412290| 0.18872740 
0.44242065 |-0.53272451|0.4360999 |-0.52511361 

measurement | result | notes 
|---|---|---|
EMM - Covariance Unstandardized | -0.0909834 |
Raw means - covariance Unstandardized |-0.08840226 | difference of -0.00258114
EMM - Covariance Standardized | -0.5113885 | standardized by sd of Gmeans: 0.421799
Raw means - covariance Standardized |-0.5113885 | standardized by sd of Gmeans: 0.4157729 
EMM - GxE | 0.3260144 |
Raw means - GxE | 0.3213567 |difference of 0.0046577
