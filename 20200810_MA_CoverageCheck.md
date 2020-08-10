# Parameter coverage check. 

Before going on vacation, I fiddled with the parameter generation to pick a set of parameters that gave full coverage while keeping a reasonable lid on the total number of jobs.

I ended up with the following parameters: 

```{params}
param_list <- list( 
  reps = c(100), 
  delta_env = c(0.01,0.25,0.5,0.75,1.0), # Set as function of n_pop. n
  delta_gen = c(-1,0.01,1),
  sample_size = c(3,5,10), 
  n_pop = c(2,3,5,10),
  env_scenario = c(1,2),# 1 = npop = n_env; 2 = multiple pops across 2 envs
  std_dev= c(2.5,5),
  interaction = NULL) # Set as function of n_pop
```
### Important notes on parameters: 

## Environmental Scenarios - 
I am running two scenarios: the original plan in which the number of populations = number of environments (1 pop to 1 env); and the second scenario in which there are 2 environments but multiple populations from a single environment. 
For the second scenario, I capped the max number of populations at 10, so there can only be a max of 5 pops per environment. I did this both for realism and computational speed.


## Interaction term - 
The number of interaction terms is the same as the number of populations (n_pop = 2 would produce 2 interaction terms, whereas n_pop = 10 produces 10 interaction terms). This gave good coverage for GxE, though it also means more datapoints in larger n_pop scenarios. 

