
slope_generation_cat <- function(genphenenv_df){
  #Input: Takes simulated, categorical, MEAN data 
  #Output: predicted slopes 
  
  # Output Dataframes
  sample_data <- data.frame()
  
  for(z in 1:length(unique(genphenenv_df$index))){
    subdat <- filter(genphenenv_df, index == unique(genphenenv_df$index)[z])
  
    # Sampled dataset from means and SE
  for(a in 1:length(unique(subdat$gen_factor))){
    for(b in 1:length(unique(subdat$exp_env_factor))){
      g <- unique(subdat$gen_factor)[a]
      e <- unique(subdat$exp_env_factor)[b]
      gdat <- filter(subdat, gen_factor == g)
      edat <- filter(gdat, exp_env_factor == e)
      
      sample_data. <- data.frame("index" = unique(edat$index),
                                 "phen_data" = (rnorm(1, edat$phen_data, edat$phen_mean_SE)),
                                 "gen_factor" = g,
                                 "exp_env_factor" = e)
      sample_data <- rbind(sample_data,sample_data.)
    }
  }
  
  # Generate slopes 
  slope_temp1 <- aov(phen_data ~ gen_factor * exp_env_factor, data = sample_data)
  simple_slopes(slope_temp1)
  G1_slope <- (slope_temp1[[1]][[3]])
  G2_slope <- (slope_temp1[[1]][[3]]+slope_temp1[[1]][[4]])
  slope_temp4 <- data.frame("gen_factor" = unique(sample_data$gen_factor),
                            "slope" = c(G1_slope,G2_slope))
  
  # Predict for Gmeans/Emeans
  GE_data <- data.frame("index" = sample_data$index,
                        "phen_data" = predict(slope_temp1),
                        "gen_factor" = sample_data$gen_factor,
                        "exp_env_factor" = sample_data$exp_env_factor)
  }
  return(list(slope_temp4,GE_data))                      
}
