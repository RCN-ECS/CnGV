
slope_generation_cat <- function(genphenenv_df){
  #Input: Takes simulated, categorical, MEAN data 
  #Output: predicted slopes 
  
  # Output Dataframes
  sample_data <- data.frame()
  
  # Sampled dataset from means and SE
  for(a in 1:length(unique(genphenenv_df$gen_factor))){
    for(b in 1:length(unique(genphenenv_df$exp_env_cat))){
      g <- unique(genphenenv_df$gen_factor)[a]
      e <- unique(genphenenv_df$exp_env_cat)[b]
      gdat <- filter(genphenenv_df, gen_factor == g)
      edat <- filter(gdat, exp_env_cat == e)
      
      sample_data. <- data.frame("phen_data" = (rnorm(1, edat$phen_data, edat$phen_mean_SE)),
                                 "gen_factor" = g,
                                 "exp_env_cat" = e)
      sample_data <- rbind(sample_data,sample_data.)
    }
  }
  
  # Generate slopes 
  slope_temp1 <- aov(phen_data ~ gen_factor * exp_env_cat, data = sample_data)
  G1_slope <- (slope_temp1[[1]][[3]])
  G2_slope <- (slope_temp1[[1]][[3]]+slope_temp1[[1]][[4]])
  slope_temp4 <- data.frame("gen_factor" = unique(sample_data$gen_factor),
                            "slope" = c(G1_slope,G2_slope))
  
  # Predict for Gmeans/Emeans
  GE_data <- data.frame("phen_data" = predict(slope_temp1),
                        "gen_factor" = sample_data$gen_factor,
                        "exp_env_cat" = sample_data$exp_env_cat)
  
  return(list(slope_temp4,GE_data))                      
}
