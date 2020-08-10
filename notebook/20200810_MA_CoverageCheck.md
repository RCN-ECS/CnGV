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
## Important notes on parameters: 

### Environmental Scenarios - 
I am running two scenarios: the original plan in which the number of populations = number of environments (1 pop to 1 env); and the second scenario in which there are 2 environments but multiple populations from a single environment. 
For the second scenario, I capped the max number of populations at 10, so there can only be a max of 5 pops per environment. I did this both for realism and computational speed.


### Interaction terms - 
The number of interaction terms is the same as the number of populations (n_pop = 2 would produce 2 interaction terms, whereas n_pop = 10 produces 10 interaction terms). This gave good coverage for GxE, though it also means more datapoints in larger n_pop scenarios. 

### Seeds 
I created a different master seed for each set of parameters. The I used the master seed with the following code to generate the 3 individual seeds and 3 sets of seeds needed for various bootstrap/permutations. Thus, all places where data are randomly generated should have an assigned seed and therefore be repeatable. 
```{seed}
set.seed(seed)
sim_seeds <- round(runif(4*n_boot)*1000000) # More than enough
seed1 = sim_seeds[1] # df.foundations
seed2 = sim_seeds[2] # df.foundations
seed3 = sim_seeds[3] # mean_gxe
seed.set1 = sim_seeds[c(4:(4+n_boot))] # Bootstrap Means seeds
seed.set2 = sim_seeds[c((5+n_boot):(5+2*n_boot))] # Permutation means set 1
seed.set3 = sim_seeds[c((6+2*n_boot):(6+3*n_boot))] # Permutation means set 1
```

## Cut to the chase, how's the coverage? 

Looks good! This is just 20% of the sims, so I believe we will have good coverage of the n_pop = 2 by the end (even though it looks sparse now). 

### Scenario 1 
N_pop = N_env

**Hex Plot:** 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/8.10.hex.png)

**Cov and GxE plot** 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/8.10.CovGxE.png)

### Scenario 2 
N_env = 2; n_pop = # of total populations (divide by 2 to get number of pops per environment)

**Hex Plot:**

This one looks a bit wonky. I'll look into why we aren't getting covGE values near 1 and -1 in n_pop = 10.
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/8.10.2envHex.png)


**Cov and GxE plot**
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/8.10.2envCovGXE.png)
