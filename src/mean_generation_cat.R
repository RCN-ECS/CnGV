mean_generation_cat <- function(genphenenv_df){
  #Input: Takes simulated, categorical, data 
  #Output: new data
  
  # Output Dataframes
  sample_data <- data.frame()
  
  # Indexing
  for(z in 1:length(unique(genphenenv_df$index))){
    subdat <- filter(genphenenv_df, index == unique(genphenenv_df$index)[z])
    
    # New dataset from means and SE
    new_data = data.frame()
    for(a in 1:length(unique(subdat$gen_factor))){
      for(b in 1:length(unique(subdat$exp_env_factor))){
        g <- unique(subdat$gen_factor)[a]
        e <- unique(subdat$exp_env_factor)[b]
        gdat <- filter(subdat, gen_factor == g)
        edat <- filter(gdat, exp_env_factor == e)

        new_data. <- data.frame("index" = rep(unique(edat$index),edat$phen_n),
                                   "phen_data" = (rnorm(edat$phen_n,edat$phen_data,edat$phen_mean_SE)),
                                   "gen_factor" = rep(g,edat$phen_n),
                                   "exp_env_factor" = rep(e,edat$phen_n))
        new_data <- rbind(new_data,new_data.)
      }
    }
    sample_data <- rbind(sample_data, new_data)
  }
  return(sample_data) 
}
