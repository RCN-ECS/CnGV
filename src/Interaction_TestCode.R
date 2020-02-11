
# Starting list of parameters
param_list <- list(
  reps = 100,
  delta_env = c(1,-1), # the amount the phenotype changes across 1 value of the environment (i.e., the slope). This is essentially the amount/degree of phenotypic plasticity that is the same across genotypes.
  delta_gen = c(1,-1), # the amount the phenotype changes from one genotype to the next. This is essitially the increase intercept from one genotype to the next.
  sample_size = c(5), 
  n_genotypes = c(2),
  n_environments = NULL,
  std_dev= c(0.0), # Random noise, with standard deviation of 1,
  interaction= c(0,5)) # this sd determines the amount of GxE)

# Table of parameters
table_fun <- function(param_list){
  
  # Basic parameters
  param_temp <- expand.grid("delta_env" = param_list$delta_env,
                            "delta_gen" = param_list$delta_gen,
                            "sample_size" = param_list$sample_size,
                            "n_genotypes" = param_list$n_genotypes,
                            "std_dev" = param_list$std_dev,
                            "interaction" = param_list$interaction)
  # Book keeping rows
  n_combo <- length(param_list$delta_env)*
    length(param_list$delta_gen)*
    length(param_list$sample_size)*
    length(param_list$n_genotypes)*
    length(param_list$std_dev)*
    length(param_list$interaction)
  reps <- rep(c(1:param_list$reps), each = n_combo)
  
  # Final data frame
  param_table <- data.frame("row"= seq(1:length(reps)), 
                            "reps" = reps,
                            "index" = rep(seq(1:n_combo), times = param_list$reps),
                            param_temp)
  
  return(param_table)
}
df = table_fun(param_list)
dim(df)

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
    int <- rnorm(param_table$n_genotypes[i] * n_environments,param_table$interaction[i], sd =0) # sd determines the amount of GxE
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
    model_df$row = rep(param_table$row[i],length(model_df$phen))
    
    output = rbind(output, model_df)
    
  }
  return(output)
}
temp_df = data_gen(df,0)  

########
# Plot #
########

for(i in 1:nlevels(factor(temp_df$row))){
  plot_df <- filter(temp_df, row == unique(row)[i])
  genenv = paste(plot_df$gen_factor,plot_df$exp_env_factor,sep= "_")
  
  a <- ggplot(plot_df, aes(x = exp_env_factor, y = phen, group = genenv,colour = gen_factor)) + 
    geom_boxplot() + theme_classic() + ylim(-65,85)+ 
    ggtitle(paste0("Interaction = ", unique(plot_df$interaction)))
print(a)
  }















    