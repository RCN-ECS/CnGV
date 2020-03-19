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
