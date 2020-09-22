#Simulation Results

Last week I started what is *hopefully* the final simulation. I settled on the following set of parameters: 

I capped the total samples (calculated as sample size * number of populations * number of environments) at 256 to reduce computational time.
I also eliminated all rows in which both delta_env and delta_gen == zero. 

```{params}
param_list <- list( 
  reps = c(10), 
  delta_env = c(0.0,0.25,0.5,.75,1.0), 
  delta_gen = c(-1,0.0,1),
  sample_size = c(2,4,8,16), 
  n_pop = c(2,4,8,16),
  env_scenario = c(1,2),  # 1 = npop = n_env; 2 = multiple pops across 2 envs
  std_dev= c(.5, 1), # Scale
  interaction = NULL) # Set as function of n_pop
```
