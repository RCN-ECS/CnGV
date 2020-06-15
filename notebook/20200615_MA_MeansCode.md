# Trouble-shooting: 

Here are 2 sets of parameters, identical except for interaction term. 

```{r}
# Expect GxE
args = c("row" = 84,"replicate" = 84, "delta_env" = 1, "delta_gen" = -1, "sample_size" = 5, "n_env" = 3, "std_dev" = 1.5, "n_pop" = 3, "interaction" = 3)

# Do not expect GxE
args = c("row" = 84,"replicate" = 84, "delta_env" = 1, "delta_gen" = -1, "sample_size" = 5, "n_env" = 3, "std_dev" = 1.5, "n_pop" = 3, "interaction" = 0)
```
You need to enter "args" but then everything you need to run the full sim is in the R files below.

Updated code: [Code](https://github.com/RCN-ECS/CnGV/blob/master/src/Cov_GxE_clusterFun.R)
Updated functions: [Func-y town](https://github.com/RCN-ECS/CnGV/blob/master/src/Cov_GxE_functions.R)

I hope you can solve the issue quickly! 
