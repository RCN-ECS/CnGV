#Sim Results

I ran a full round (10 replicates) for the following starting parameters. [Recall I am using a different way to set up parameters.](https://github.com/RCN-ECS/CnGV/blob/master/notebook/20201006_NewParamGeneration.md)

I have made the following plots, ordered in terms of where they will go in the paper or supplementary material. Figures 1 and 2 are handmade heuristics. They can be found [here.](https://docs.google.com/document/d/1CcoJFTX6I8zPptLAzBvITiTaTkNtdgGxNNdQlrcxlSY/edit#heading=h.68abrs8dfl95)
```{params}
param_list <- list( 
  reps = c(10), 
  delta_env = NULL, 
  delta_gen = NULL,
  sample_size = c(2,4,8,16), 
  n_pop = c(2,4,8,16),
  env_scenario = c(1,2),  # 1 = npop = n_env; 2 = multiple pops across 2 envs
  std_dev= c(.5, 1), # Scale
  interaction = NULL) 
```

## Figure 3
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Simulation_10.10.2020/PhenotypePanel.pdf)
