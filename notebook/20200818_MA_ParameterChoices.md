# Candidate Parameter Sets

| N_pop | N_Env | samplesize | Total Number of samples |
| --- | --- | --- | --- | 
| 2 | 2 | 16 | 64 |
| 4 | 2 | 8 | 64 |
| 4 | 4 | 4 | 64 |
| 8 | 2 | 4 | 64 |
| 16| 2 | 2 | 64 |
| 2 | 2 | 32 | 128 |
| 4 | 2 | 16 | 128 |
| 4 | 4 | 8 | 128 |
| 8 | 2 | 8 | 128 |
| 8 | 8 | 2 | 128 |
| 16 | 2 | 4 | 128 |
| 32 | 2 | 2 | 128 |
| 2 | 2 | 64 | 256 |
| 4 | 2 | 32 | 256 |
| 8 | 2 | 16 | 256 |
| 16 | 2 | 8| 256 |
| 32 | 2 | 4 | 256 |
| 4 | 4 | 16 | 256 |
| 8 | 8 | 4 | 256 |
| 4 | 4 | 16 | 256 |


So, start with the following set, and cull down to just those with total sample sizes = 64, 128, or 256

```{param}
# Starting list of parameters
param_list <- list( 
  reps = c(10), 
  delta_env = c(0.0,0.25,0.5,.75,1.0), 
  delta_gen = c(-1,0.0,1),
  sample_size = c(2,4,8,16,32), 
  n_pop = c(2,4,8,16,32),
  env_scenario = c(1,2),  # 1 = npop = n_env; 2 = multiple pops across 2 envs
  std_dev= c(.5, 1), # Scale
  interaction = NULL) # Set as function of n_pop


```
