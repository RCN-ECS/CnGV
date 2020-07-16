# Troubleshooting 

Last week, Katie asked me if there was anythign else to do before finally setting up the BIG simulation. I ran the following limited set of parameters and got the following result: 

```{params}
# Starting list of parameters
param_list <- list( 
  reps = c(5), # or more?
  delta_env = c(0.01,0.5,1),
  delta_gen = c(-1,0.01,1),
  sample_size = c(5,10),#c(5,10,20),#c(5,10,20), 
  pop_per_env = c(1),#c(2,3,4,5),#c(2,3,5,10,15), 
  n_environments = c(2,5,10),#c(2,3,5,7,10),
  std_dev= c(1,3),#c(0.5,1), 
  interaction = NULL) # Length = number of populations, number of levels = number of populations
```
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.16.WonkyResult.png)

![image](https://tenor.com/view/sad-face-rain-doctor-who-gif-11466849)

The above plot is actually less bad than previous attempts. I had used `runif()` to choose the levels of interaction, and thought originally that the error was somewhere in my parameter generation code. But after changing several things, I have found that the wonkiness remains. So it must be something else. 


