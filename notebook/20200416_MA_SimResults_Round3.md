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
First, 4 vignettes to prove that everything is working ok: 

| --- | Negative Covariance | Positive Covariance | 
| --- | --- | --- | 
| Low GxE| https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/plot_row5.png |https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/plot_row7.png |
| High GxE|https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/plot_row61.png | https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/plot_row59.png|
