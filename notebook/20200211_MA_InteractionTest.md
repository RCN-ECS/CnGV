## Quick background

The new simulation code generates data based on the anova model (as opposed to generating data based on linear model formula). The reason for the switch is to allow increased freedom in testing >2 genotypes across >2 environments.

The data generation code takes in each row from the table of starting parameters to create the data using the below code: (which is subsetted out from the larger "ring" function that simultaneously simulates data and calculates covariance functions/gxe functions, boot straps, and permutations.

The degree of interaction is set by generating an interaction number by sampling a normal distribution with a mean of 0 and standard deviation matching the starting interaction parameter. 
Thus, the higher the starting interaction parameter, the greater the variation among genotypes and environments *should* be. 

```(data generation)
data_gen <- function(param_table, n_boot){

  # Load packages
  library("emmeans")
  library("lme4")
  library("tidyverse")
  library("rlist")
  
  # Output dataframe
  output <- data.frame()
  
  for(i in 1:nrow(param_table)){
    
    # Counter
    cat(i, "\n")
    
    n_environments = param_table$n_genotypes[i]
    
    # Approximate Cov(G,E)
    cov_GE_approx = param_table$delta_env[i] * param_table$delta_gen[i]
    
    gen <- rep(1:param_table$n_genotypes[i], each = param_table$sample_size[i])
    env <- rep(1:n_environments, times = param_table$sample_size[i]) 
    
    noise <- rnorm(param_table$sample_size[i] * param_table$n_genotypes[i], 0, sd = param_table$std_dev[i]) # Random noise
    
    # Create Interactions
    int <- rnorm(param_table$n_genotypes[i] * n_environments, 0, sd = param_table$interaction[i]) # sd determines the amount of GxE
    int_df <- data.frame(expand.grid(G = 1:param_table$n_genotypes[i], E = 1:n_environments), int)
    
    # Create the model dataframe 
    model_df <- data.frame(gen, env, noise)
    model_df <- merge(model_df, int_df)
    model_df$gen_factor = factor(paste("G", model_df$G, sep = "_"))
    model_df$exp_env_factor = factor(paste("E", model_df$E, sep = "_"))
    
    # Generate phenotype data using regression equation
    phen = param_table$delta_env[i] * model_df$E + param_table$delta_gen[i] * model_df$G  + model_df$noise + model_df$int
    model_df$phen = phen
    
    # Standardize data
    dat_avg <- mean(phen) 
    dat_std <- sd(phen)
    model_df$phen_corrected <- ((phen - dat_avg)/dat_std)
    
    model_df$interaction = rep(param_table$interaction[i],length(model_df$phen))
    model_df$rep = rep(param_table$reps[i],length(model_df$phen))
    
    output = rbind(output, model_df)
  }
  return(output)
}    

```
To test this, I ran 5 levels of interaction, from 0 to 4. Here are my starting parameters: 

```(starters)
# Starting list of parameters
param_list <- list(
  reps = 100,
  delta_env = c(1), # the amount the phenotype changes across 1 value of the environment (i.e., the slope). This is essentially the amount/degree of phenotypic plasticity that is the same across genotypes.
  delta_gen = c(1), # the amount the phenotype changes from one genotype to the next. This is essitially the increase intercept from one genotype to the next.
  sample_size = c(5), 
  n_genotypes = c(2),
  n_environments = NULL,
  std_dev= c(0.0), # Random noise, with standard deviation of 1,
  interaction= c(0,1,2,3,4)) # this sd determines the amount of GxE)

```
When I do that and create box plots of the phenotypic data generated for each interaction, these are the plots I get: 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/int0.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/int1.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/int2.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/int3.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/int4.png)

