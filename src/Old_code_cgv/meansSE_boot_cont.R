
meansSE_boot_cont <- function(genphenenv_df,iterations){
  
  # Bootstrap replicates
  iter <- iterations

  # Output dfs
  temp_slope_dat <- data.frame()
  model_data <- data.frame()
  temp_GE_dat <- data.frame()
  pred_dat <- data.frame()
  
  # Indexes
  for(x in 1:length(unique(genphenenv_df))){
    ind_dat <- filter(genphenenv_df, index == unique(genphenenv_df$index)[x])
  
    # Bootstrap 
  for(d in 1:iter){
    temp_slope_dat. = slope_generation_cont(genphenenv_df)
    temp_slope_dat = rbind(temp_slope_dat.[[1]],temp_slope_dat)
    temp_GE_dat = rbind(temp_slope_dat.[[2]],temp_GE_dat)
    }
  
  # Standardize bootstrapped data
  phen_mean = mean(temp_GE_dat$phen_data) 
  phen_sd = sd(temp_GE_dat$phen_data) 
  env_avg = mean(temp_GE_dat$exp_env_cont) 
  env_sd = sd(temp_GE_dat$exp_env_cont)
  nat_avg = mean(genphenenv_df$nat_env_cont)
  nat_sd = sd(genphenenv_df$nat_env_cont)
  temp_GE_dat$phen_corrected = ((temp_GE_dat$phen_data-phen_mean)/phen_sd)
  temp_GE_dat$env_corrected = ((temp_GE_dat$exp_env_cont-env_avg)/env_sd) 
  temp_GE_dat$native_env_corrected = ((genphenenv_df$nat_env_cont-nat_avg)/nat_sd) 
  
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
    for(g in 1:length(unique(temp_GE_dat$exp_env_cont))){
      gen <- unique(temp_GE_dat$gen_factor)[f]
      env <- unique(temp_GE_dat$exp_env_cont)[g]
      gen_dat <- filter(temp_GE_dat, gen_factor == gen)
      env_dat <- filter(gen_dat, exp_env_cont == env)
      CI <- quantile(env_dat$phen_corrected, probs=c(0.025, 0.975), type=1)
      avg_val <- mean(env_dat$phen_corrected)
      
      pred_dat. = data.frame("gen_factor" = unique(env_dat$gen_factor),
                             "exp_env_cont" = unique(env_dat$exp_env_cont),
                             "mean_phen" = avg_val,
                             "mean_lwrCI" = CI[[1]],
                             "mean_uprCI" = CI[[2]])
      pred_dat = rbind(pred_dat, pred_dat.) 
      }
      }
  
    # Match 
    pred_dat$nat_env_cont <- genphenenv_df$nat_env_cont[match(pred_dat$gen,genphenenv_df$gen_factor)]
    pred_dat$native_env_corrected <- temp_GE_dat$native_env_corrected[match(pred_dat$gen,temp_GE_dat$gen_factor)]
    
    G1mean = (pred_dat[1,3]+pred_dat[2,3])/nrow(pred_dat[pred_dat$gen_factor == "G1",])
    G2mean = (pred_dat[3,3]+pred_dat[4,3])/nrow(pred_dat[pred_dat$gen_factor == "G2",])
    E1mean = (pred_dat[1,3]+pred_dat[3,3])/nrow(pred_dat[pred_dat$exp_env_cont == "-2",])
    E2mean = (pred_dat[2,3]+pred_dat[4,3])/nrow(pred_dat[pred_dat$exp_env_cont == "-2",])
    
    Cov_matrix. = data.frame("gen" = unique(pred_dat$gen_factor),
                            "native_env" = unique(pred_dat$nat_env_cont),
                            "G_means" = c(G1mean,G2mean),
                            "E_means" = c(E1mean,E2mean))
    Cov_matrix = rbind()
    (cov_est = cov(Cov_matrix$G_means,Cov_matrix$E_means))
  }
  
    return(list(model_data,Cov_matrix,cov_est))
}

 