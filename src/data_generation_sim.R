
# Starting parameters
init_params <- list(
  delta_env = 1, # the amount the phenotype changes across 1 value of the environment (i.e., the slope). This is essentially the amount/degree of phenotypic plasticity that is the same across genotypes.
  delta_gen = -1, # the amount the phenotype changes from one genotype to the next. This is essitially the increase intercept from one genotype to the next.
  sample_size = 5, 
  n_genotypes = 3,
  n_environments = n_genotypes,
  std_dev= 0.01, # Random noise, with standard deviation of 1,
  interaction=1) # this sd determines the amount of GxE)

data_generation <- function(init_params){
  
  # Parameters
  delta_env = init_params$delta_env
  delta_gen = init_params$delta_gen
  sample_size = init_params$sample_size
  n_genotypes = init_params$n_genotypes
  n_environments = init_params$n_environments
  std_dev = init_params$std_dev
  interaction = init_params$interaction
  
  # Output dataframe
  output <- data.frame()
  row = 0 #for indexing
  
  # Generate data
  for(a in 1:length(delta_env)){
    for(b in 1:length(delta_gen)){
      for(c in 1:length(sample_size)){
        for(d in 1:length(n_genotypes)){
          for(e in 1:length(n_environments)){
            for(f in 1:length(std_dev)){
              for(g in 1:length(interaction)){
                
                row=row+1 #index
                
                # Approximate Cov(G,E)
                cov_GE_approx = delta_env[a] * delta_gen[b]
                
                gen <- rep(1:n_genotypes[d], each = sample_size[c])
                env <- rep(1:n_environments[e],times = sample_size[c]) 
                
                noise <- rnorm(sample_size[c] * n_genotypes[d], 0, sd = std_dev[f]) # Random noise
                
                # Create Interactions
                int <- rnorm(n_genotypes[d] * n_environments[e], 0, sd = interaction[g]) # sd determines the amount of GxE
                int_df <- data.frame(expand.grid(G = 1:n_genotypes[d], E = 1:n_environments[e]), int)
                
                # Create the model dataframe 
                model_df <- data.frame(G, E, noise)
                model_df <- merge(model_df, int_df)
                
                # Generate phenotype data using regression equation
                phen = delta_env[a] * model_df$E + delta_gen[b] * model_df$G  + model_df$noise + model_df$int
                
                # Output
                temp_data <- data.frame("source" = rep("sim", (length(phen))),
                                        "index" = rep(row, (length(phen))),
                                        "GxE_approx" = rep(interaction[g], (length(phen))),
                                        "Cov_approx" = rep(cov_GE_approx, (length(phen))),
                                        "n_environments" = rep(n_environments[e], (length(phen))),
                                        "phen_n" = rep(sample_size[c], (length(phen))),
                                        "stdev" = rep(std_dev[f], (length(phen))),
                                        "nat_env_factor" = paste("E", model_df$G, sep = "_"),
                                        "gen_factor" = paste("G", model_df$G, sep = "_"),
                                        "exp_env_factor" = paste("E", model_df$E, sep = "_"),
                                        "phen_data" = phen)

                output = rbind(output, temp_data)
              }
            }
          }
        }
      }
    }
  }
  return(output)
}
