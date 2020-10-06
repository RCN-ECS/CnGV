# New Code to fill parameter space

In my previous attempts to generate starting parameters that cover the full parameter space that we desired, I had a ground up approach. I started with expand.grid for basic parameters and then added more layers as needed on top of that grid. This approach made it difficult to limit the size of the resulting dataframe since every added variable would double or triple the size. The final death knell came when I realized I needed to add another parameter (as referenced[here](https://github.com/RCN-ECS/CnGV/blob/master/notebook/20200930_MA_NewDataGeneration_2Envs.md)) to deal with the situation in env_scenario = 2. Except I could not figure out a way to do this without tripling an already way too big dataframe. 

# Flip the Script
I decided to try a completely new approach, re-write the code, and work backwards instead. I decided to choose an upper bounds on the dataframe size and work backwards to fill in the parameter space. I set the limit at 16,000 parameter sets per env_scenario (so the max number of simulations would be 32,000). 

We are ultimately interested in how different population sizes (4 levels), sample sizes (4 levels), and standard deviation (2 levels) affect covariance and GxE estimates and power. If I capped the dataframe at 16000, that would give me 10 replicates of 50 samples across those 32 levels. 

Simply put, I generated 50 values for delta_env, delta_gen, and interaction levels for each replicate. For env_scenario2, I generated 50 values for delta_env, delta_gen, interaction, and errpop levels (errpop corrects the sampling issue) for each replicate. 

Because we had gaps in that lower true_GxE parameter space, I more heavily sampled low interactions (low_int = runif(n = 40, min = 0, max = 1.99)) relative to higher interaction levels (high_int = runif(n = 10, min = 2, max = params$n_pop[i])). Note that if n_pop = 2, that all 50 would be sampled between 0 and 2. Seeds and bookkeeping rows were added at the end. 

This was repeated twice for each env_scenario across 10 replicates. After generating parameters, I filtered the parameter dataset to those sets of parameters with fewer than 500 total samples (n_env * n_gen * sample_size). This resulted in a final dataset size of 20,000. Not too shabby! 

## Errpop Side adventure: 
For the env_scenario = 2 situation. As I mention in the previous post, I set errpop = rnorm(mean = 0, sd = abs(delta_gen)). Except after running a test sim, I found that making sd in errpop = abs(delta_gen) did not correct the issue. I had 3 levels of delta_gen (-1, 0, 1). Thus the absolute value would give me a big sd for 2/3 of the cases, which again led to the majority of covariances being very low and did not extend from -1 to 1. 

I knew I wanted a normal distribution that more heavily sampled from lower values. With the new approach (detailed below), I now sample 50 delta_gens between -1 and 1 which gave a nice, relatively uniform distribution between 0 and 1 (absolute value). When I fed this vector of delta_gens into the following and took the absolute value of the distribution to bound it at 0 (errpop = **abs**(rnorm(n = 50,  mean = 0, sd = abs(delta_gen))), it gave a nice chi-squared-esque distribution.

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/10.6.ErrpopHist.png)

Yay! But the real test is whether this approach actually does a good job with covering parameter space, so I ran a **single** replicate to test. 

### Hex Plots

Env_scenario 1. is the covariance and gxe space reasonably well covered with no systematic gaps? YES! 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/10.6.Hex1.png)

Env_scenario 2. This one had me nervous. Did I finally resolve the issue and sample from -1 to 1? YES!
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/10.6.Hex2.png)

### Estimated Values

Env_scenario 1 - looks like even though the bigger n_pops start to float off zero when error is introduced, that we are still catching a good range. 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/10.6.GxEvCov1.png)

Env_scenario 2 - 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/10.6.GxEvCov2.png)
