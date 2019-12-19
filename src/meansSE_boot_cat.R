
meansSE_boot_cat <- function(genphenenv_df,iterations){
  
  # Bootstrap replicates
  iter <- iterations
  
  # Output dfs
  sample_data <- data.frame()
  temp_slope_dat <- data.frame()
  model_data = data.frame()
  temp_GE_dat = data.frame()
  pred_dat = data.frame()
  
  # Bootstrap 
  for(d in 1:iter){
    temp_slope_dat. = slope_generation_cat(genphenenv_df)
    temp_slope_dat = rbind(temp_slope_dat.[[1]],temp_slope_dat)
    temp_GE_dat = rbind(temp_slope_dat.[[2]],temp_GE_dat)
  }
  
  # Standardize bootstrapped data
  phen_mean = mean(temp_GE_dat$phen_data) 
  phen_sd = sd(temp_GE_dat$phen_data) 
  temp_GE_dat$phen_corrected = ((temp_GE_dat$phen_data-phen_mean)/phen_sd)
  
  # Model coefficient confidence intervals
  for(e in 1:length(unique(temp_slope_dat$gen_factor))){
    gen <- unique(temp_slope_dat$gen_factor)[e]
    gen_dat <- filter(temp_slope_dat, gen_factor == gen)
    slope_CI <- quantile(gen_dat[,2], probs=c(0.025, 0.975), type=1) # Slope confidence intervals 
    slope_avg <- mean(gen_dat[,2]) # Mean should be similar to original slope estimate
    model_data. = data.frame("gen_factor" = gen,
                             "slope_average" = slope_avg,
                             "slope_lwrCI" = slope_CI[[1]],
                             "slope_uprCI" = slope_CI[[2]])
    model_data = rbind(model_data,model_data.)
  }
  
  # Predicted data
  for(f in 1:length(unique(temp_GE_dat$gen_factor))){
    for(g in 1:length(unique(temp_GE_dat$exp_env_cat))){
      gen <- unique(temp_GE_dat$gen_factor)[f]
      env <- unique(temp_GE_dat$exp_env_cat)[g]
      gen_dat <- filter(temp_GE_dat, gen_factor == gen)
      env_dat <- filter(gen_dat, exp_env_cat == env)
      CI <- quantile(env_dat$phen_corrected, probs=c(0.025, 0.975), type=1)
      avg_val <- mean(env_dat$phen_corrected)
      
      pred_dat. = data.frame("gen_factor" = unique(env_dat$gen_factor),
                             "exp_env_cat" = unique(env_dat$exp_env_cat),
                             "mean_phen" = avg_val,
                             "mean_lwrCI" = CI[[1]],
                             "mean_uprCI" = CI[[2]])
      pred_dat = rbind(pred_dat, pred_dat.) 
    }
  }
  
  # Match with native environment
  pred_dat$nat_env_cat <- genphenenv_df$nat_env_factor[match(pred_dat$gen,genphenenv_df$gen_factor)]
  ngen <- length(unique(pred_dat$gen_factor))
  nenv <- length(unique(pred_dat$exp_env_cat))
  
  # Gmeans and Emeans
  G1mean = sum(pred_dat[pred_dat$gen_factor == "G1",3])/ngen
  G2mean = sum(pred_dat[pred_dat$gen_factor == "G2",3])/ngen
  E1mean = sum(pred_dat[pred_dat$exp_env_cat == "E1",3])/nenv
  E2mean = sum(pred_dat[pred_dat$exp_env_cat == "E2",3])/nenv
  
  Cov_matrix = data.frame("gen" = unique(pred_dat$gen_factor),
                          "native_env" = unique(pred_dat$nat_env_cat),
                          "G_means" = c(G1mean,G2mean),
                          "E_means" = c(E1mean,E2mean))
  (cov_est = cov(Cov_matrix$G_means,Cov_matrix$E_means))
  
  return(list(model_data,Cov_matrix,cov_est))
}

