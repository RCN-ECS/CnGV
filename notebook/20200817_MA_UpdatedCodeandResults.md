Last week, we discovered that the way we were simulating and standardizing data could affect our "true" estimates. Katie designed a new way of generating simulated data shown [here](https://github.com/RCN-ECS/CnGV/blob/master/notebook/20200811_KEL_NewWayPhenData.Rmd).

I implemented this approach into the way we simulate data. Here is the new function that generates and standardizes phenotypic data (for raw data with error, raw data with no error, means data with error, means data with no error):

```{new}
df.foundations <- function(delta_env, delta_gen, sample_size, n_env, std_dev, n_pop, interaction, seed1, seed2){
  
  # Dataframe foundations
  n_environments <- n_env 
  gen <- rep(1:n_pop, each = sample_size*n_environments)
  env <- rep(1:n_environments, times = n_pop, each = sample_size) 
  
  # Interaction Term
  set.seed = seed2
  int <- rep(rnorm(n_pop * n_environments, 0, sd = interaction), each = sample_size) # interaction term - one for each GE level
  
  ############################################
  ## Phenotype data with standard deviation ##
  ############################################
  
  # Create model dataframe 
  model_df <- data.frame(gen, env, int)
  model_df$gen_factor <- factor(paste("G", model_df$gen, sep = "_"))
  model_df$exp_env_factor <- factor(paste("E", model_df$env, sep = "_"))
  model_df$GE_factor <- paste0("G",model_df$gen, "E", model_df$env)
  
  # Native environments
  model_df$nat_env_factor = rep(unique(model_df$exp_env_factor), each = (sample_size*n_pop))
  
  # True mean phenotype data using regression equation
  model_df$GE_true = delta_env * model_df$env + delta_gen * model_df$gen + model_df$int 

  # True means 
  G_true <- data.frame(G_true = tapply(model_df$GE_true, model_df$gen_factor, mean))
  G_true <- rownames_to_column(G_true, var = "gen_factor")
  E_true <- data.frame(E_true = tapply(model_df$GE_true, model_df$exp_env_factor, mean))
  E_true <- rownames_to_column(E_true, var = "exp_env_factor")
  
  model_df <- full_join(model_df, G_true, by = "gen_factor")
  model_df <- full_join(model_df, E_true, by = "exp_env_factor")
  model_df$mean_true <- mean(model_df$GE_true)
  
  # Calculate true interaction term
  model_df$int_true <- model_df$GE_true - model_df$G_true - model_df$E_true + model_df$mean_true # GEij - Gmean - Emean + Overall
  
  # Standardize true means 
  GE_true_means <- tapply(model_df$GE_true, model_df$GE_factor, mean)
  model_df$GE_stn_true <- (model_df$GE_true - mean(GE_true_means))/sd(GE_true_means)
  
  # Sanity plots
  #ggplot(model_df, aes(x = exp_env_factor, y = GE_stn_true, group = gen_factor, colour = gen_factor))+ geom_point()+geom_line()
  #ggplot(model_df, aes(x = exp_env_factor, y = GE_true, group = gen_factor, colour = gen_factor))+ geom_point()+geom_line()
  
  # Add random noise
  set.seed = seed1
  model_df$e <- rnorm(sample_size * n_pop * n_environments, 0, sd = std_dev) 
  
  model_df$phen_corrected <- model_df$GE_stn_true + model_df$e
  #ggplot(model_df, aes(x = exp_env_factor, y = phen, group = gen_factor, colour = gen_factor))+geom_smooth(se=F)
  
  # Chr -> Factor
  model_df$gen_factor <- as.factor(model_df$gen_factor)
  model_df$exp_env_factor <- as.factor(model_df$exp_env_factor)
  
  # Generate means dataframe
  GE_means <- data.frame(avg_phen_corrected = tapply(model_df$phen, model_df$GE_factor, mean))
  GE_means <- rownames_to_column(GE_means, "GE_factor")
  GE = data.frame(GE_factor = unique(model_df$GE_factor),
                  gen_factor = rep(unique(model_df$gen_factor), each = n_pop),
                  exp_env_factor = rep(unique(model_df$exp_env_factor), times = n_pop))
  mean_df <- full_join(GE, GE_means, by="GE_factor")
  
  sd_avg = model_df %>%
    group_by(gen_factor, exp_env_factor) %>%
    summarize("deviation" = mean(abs(e)))
  mean_df <- full_join(mean_df, sd_avg)
  mean_df$se <- mean_df$deviation/(sqrt(sample_size))
  mean_df$nat_env_factor <- model_df$nat_env_factor[match(mean_df$gen_factor,model_df$gen_factor)]
  
  # Chr -> Factor
  mean_df$gen_factor <- as.factor(mean_df$gen_factor)
  mean_df$exp_env_factor <- as.factor(mean_df$exp_env_factor)
  
  #ggplot(mean_df, aes(x = exp_env_factor, y = avg_phen_corrected, group = gen_factor))+geom_line(aes(colour = gen_factor))
  
  ###############################################
  ## Phenotype data with no standard deviation ##
  ###############################################
  
  # Raw data
  df.ne = data.frame(gen_factor = model_df$gen_factor, 
                     exp_env_factor = model_df$exp_env_factor, 
                     nat_env_factor = model_df$nat_env_factor,
                     phen_corrected = model_df$GE_stn_true)
  
  # Generate means dataframe
  GE_means_ne <- data.frame(avg_phen_corrected = tapply(model_df$GE_stn_true, model_df$GE_factor, mean))
  GE_means_ne <- rownames_to_column(GE_means_ne, "GE_factor")
  mean_df_ne <- full_join(GE, GE_means_ne, by="GE_factor")
  mean_df_ne$se <- 0
  
  # Add native environments
  mean_df_ne$nat_env_factor = df.ne$nat_env_factor[match(mean_df_ne$gen_factor,df.ne$gen_factor)]
  
  # Chr -> Factor
  mean_df_ne$gen_factor <- as.factor(mean_df_ne$gen_factor)
  mean_df_ne$exp_env_factor <- as.factor(mean_df_ne$exp_env_factor)
  
  return(list(model_df,mean_df,df.ne,mean_df_ne))
}
```

After setting that up and running some smaller tests, I ran 10 replicates of the following limited set of parameters: 

```{param}
param_list <- list( 
  reps = c(10), 
  delta_env = c(0.01, 0.5, 1),
  delta_gen = c(-1,0.01,1),
  sample_size = c(3,5),
  n_pop = c(2,5,10),
  env_scenario = c(1),  # 1 = npop = n_env; 2 = multiple pops across 2 envs
  std_dev= c(0.125,.25,.5, 1), # Scale (1*sd
  interaction = NULL) # set as function of n_pop
```

This yielded the following. 

**A change in Color Scheme** 
CovGE is considered significant if 95% Confidence Intervals do not intercept 0. GxE is considered significant if p-value is < 0.05. 
Red = Both are significant, green = just covariance is significant, blue = just GxE is significant.

### Covariance Vs. GxE Estimated Marginal Means
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/8.17.CovGxE.png)

### Covariance vs. GxE Omega Squared
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/8.17_CovOmega.png)

### GxE Estimated Marginal means vs. Omega Squared
Same result, Omega squared does not equal the estimated marginal means approach.
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/8.17.GxEvsOmega2.png)

### Covariance - Means vs. Raw 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/8.17.CovMeans_Raw.png)

### Covariance - True estimate (no error) vs. estimate with error
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/8.17.CovVsTrueCov.png)

### Covariance - Means 95% CI vs. Raw 95% CI
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/8.17.CovError.png)

### GxE - Means vs. Raw
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/8.17.GxEmeansvsRaw.png)

### GxE - True estimate (no error) vs. estimate with error
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/8.17.GxE_errorvs.True.png)

### GxE - Means 95% CI vs. Raw 95% CI
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/8.17.GxECIs.png)
