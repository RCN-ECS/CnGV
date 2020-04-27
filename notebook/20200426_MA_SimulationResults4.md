# Full simulation Results

After resolving the covariance correction [issue](https://github.com/RCN-ECS/CnGV/blob/master/notebook/20200416_MA_SimResults_Round3.md), I ran ten replicates of the full array of parameters which were: 

```{param list}
param_list <- list( 
  reps = c(10), # or more?
  delta_env = c(0,0.5,1),
  delta_gen = c(-1,0,1),
  sample_size = c(5,10,20), 
  n_pop = c(2,3,5,10,15), 
  n_environments = NULL,
  std_dev= c(0.5,1), 
  interaction = 5) # Vector LENGTH not magnitude
```
