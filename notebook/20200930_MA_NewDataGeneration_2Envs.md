# New Simulation approach

We simulated 2 scenarios. One in which the number of environments matched the number of genotypes, and one in which there were just 2 environments but multiple genotypes from each environment. 

In the results for the second scenario, my covariance estimates were bound near 0.4 and -0.4 and never reached 1 or -1 as they should.

I realized that the way we are simulating the data is confining the covariance estimates. Because I am applying delta_gen to each genotype independently, each genotype was spaced from the other genotype by a certain interval (delta_gen) even if they shared an enviroment. But if there is strong CovGE, one would expect genotypes to cluster more closely together when they share an environment. So by limiting the ability of genotypes from the same environment to cluster together, I was accidentally limiting the CovGE.

So I came up with a new way to simulate the data for the 2 environment scenarios that allows for clustering but I would like to get the official lotterhos stamp of approval to make sure I didn't make any bad decisions. 

To clarify, I only apply this new simulation approach to *JUST* the scenarios with 2 env/ >2gens. I kept the same data generation approach for the scenarios where n_env = n_gen.

## New Approach

I simulate a single average phenotype for each genotype/environment, and use that average as the anchor for genotypes that share a native environment. 

I simulate a single mean phenotype for 2 genotypes across 2 environments using the standard regression equation. I add the interaction term but no noise. The below example has interaction = 0, so there is no GxE in these example plots.

```{code1}
# True mean phenotype data using regression equation
model_df$GE_true_temp = delta_env * model_df$env + delta_gen * model_df$gen + model_df$int 
```

This produces something that looks like this: 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/930_2GE.png)

To these means, I add deviation that represents the degree that each genotype deviates from its overall mean. The bounds of that deviation are set by delta_gen. 

```{code2}
# Deviation for each genotype 
set.seed = seed3
model_df$e_pop <- rep(rnorm(n_pop, 0, sd = abs(delta_gen)), each = n_env*sample_size)
model_df$GE_true <- model_df$GE_true_temp + model_df$e_pop
```

For a 4 total genotype situation (2 genotypes per environment), this gives something like the following: 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/930_2GE_step2.png)

I then standardize these means, and as a final step, I add the individual samples by generating noise around each genotype's true means.

```{code3}
 # Standardize true means 
GE_true_means <- tapply(model_df$GE_true, model_df$GE_factor, mean)
model_df$GE_stn_true <- (model_df$GE_true - mean(GE_true_means))/sd(GE_true_means)
  
# Add random noise
set.seed = seed1
model_df$e <- rnorm(sample_size * n_pop * n_environments, 0, sd = std_dev) 
  
# Phenotype with error
model_df$phen_corrected <- model_df$GE_stn_true + model_df$e
```

To produce something like the following:
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/930_2GEstep3.png)

This particular example gives a significant but low CNGV covariance estimate of -0.28 (true_cov = -0.23) and GxE = 0.07 (true GxE = 0), which makes sense to me based on the pattern.

I have run several tests with different n_gens and delta_gens and this approach allows genotypes from the same environment to cluster, but generally keeps the genotypes from different environments apart by delta_gen. I also get better estimates of CovGE that appear (in my estimation) to visually match the patterns displayed by reaction norms.

Today I will run a replicate of the new interaction terms (finer resolution in the lower interaction terms) and 2 env_scenario on the cluster to test everything out at a larger scale. 



