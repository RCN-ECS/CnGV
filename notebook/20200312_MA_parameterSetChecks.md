# Parameter Checks

Prior to initiating a full power analysis on the cluster, I ran the set of parameters through a pared down version of the function to determine whether this set of parameters gives us the coverage we desire. 

The set of parameters is as follows: 
```{params}
param_list <- list( 
  reps = c(1),
  delta_env = c(0,0.5,1),
  delta_gen = c(-1,0,1),
  sample_size = c(5,10,20), 
  n_genotypes = c(2,4,8,16),
  n_environments = NULL,
  std_dev= c(0.1,0.5),#seq(from = 0.0, to = 2.1, by = 0.5), # Random noise
  interaction= c(0,0.5,1)) # Interaction term
```
After running the sim, I created the standard "horseshoe" plots 4x, each time coloring the points by either delta_env, delta_gen, std_dev, or interaction so we could figure out which parameter is most driving patterns and which parameters might need some filling in.

It seems like we may want a bit more resolution in our delta_envs and delta_gens.

### Delta env

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/Delta_env_trues.png)

### Delta gen

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/deltagen_trues.png)

### Standard Deviation

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/StdDev_trues.png)

### Interaction term

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/interaction_trues.png)
