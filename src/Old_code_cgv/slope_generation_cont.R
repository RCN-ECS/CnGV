
slope_generation_cont <- function(genphenenv_df){
  #Input: Takes simulated, continuous, MEAN data 
  #Output: Predicted slopes 
  
  # Output Dataframes
  sample_data <- data.frame()
  
  # Sampled dataset from means and SE
  for(a in 1:length(unique(genphenenv_df$gen_factor))){
    for(b in 1:length(unique(genphenenv_df$exp_env_cont))){
      g <- unique(genphenenv_df$gen_factor)[a]
      e <- unique(genphenenv_df$exp_env_cont)[b]
      gdat <- filter(genphenenv_df, gen_factor == g)
      edat <- filter(gdat, exp_env_cont == e)
      
      sample_data. <- data.frame("phen_data" = (rnorm(1, edat$phen_data, edat$phen_mean_SE)),
                                 "gen_factor" = g,
                                 "exp_env_cont" = e)
      sample_data <- rbind(sample_data,sample_data.)
    }
  }
  
  # Generate slopes 
  slope_temp1 <- lm(phen_data ~ gen_factor * exp_env_cont, data = sample_data)
  slope_temp2 <- simple_slopes(slope_temp1)
  slope_temp3 <- filter(slope_temp2, gen_factor != "sstest")
  slope_temp4 <- data.frame("gen_factor" = slope_temp3[,1],
                            "slope" = slope_temp3[,3],
                            "intercept" = rep(coef(slope_temp1)[[1]],length(unique(slope_temp3[,1]))))
  
  # Predict for Gmeans/Emeans
  GE_data <- data.frame("phen_data" = predict(slope_temp1),
                        "gen_factor" = sample_data$gen_factor,
                        "exp_env_cont" = sample_data$exp_env_cont)
  
  return(list(slope_temp4,GE_data))                      
}
