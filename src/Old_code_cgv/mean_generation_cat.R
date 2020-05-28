mean_generation_cat <- function(genphenenv_df){
  #Input: Takes simulated, categorical, data 
  #Output: new data
  require(tidyverse)
  
  # New dataset from means and SE
  tempdat = genphenenv_df %>%
    group_by(row, gen_factor, exp_env_factor) %>%
    summarise(new_mean = mean(phen_corrected),
              new_sd = sd(phen_corrected))
  tempdat$se = tempdat$new_sd/sqrt(unique(genphenenv_df$sample_size))
  
  
return(tempdat)
  
  }
test1 <- mean_generation_cat(temp_df)
