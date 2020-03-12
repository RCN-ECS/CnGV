# To create new df, run "table_fun" in Power_simulation_datagen2.0
  
  # Load packages
  library("emmeans")
  library("lme4")
  library("tidyverse")
  library("rlist")
  
  # Load Parameters
  args = commandArgs(trailingOnly = TRUE)
  
  row = as.numeric(args[1])
  replicate = as.numeric(args[2])
  delta_env = as.numeric(args[3])
  delta_gen = as.numeric(args[4])
  sample_size = as.numeric(args[5])
  n_genotypes = as.numeric(args[6])
  std_dev = as.numeric(args[7])
  interaction = as.numeric(args[8])
  n_boot = 100
  
  # Output dataframe
  output <- data.frame() 
  
  # For reproducibility
  set.seed = 999
    
  n_environments = n_genotypes
    
  # Approximate Cov(G,E)
  cov_GE_approx = delta_env * delta_gen
    
  # Dataframe foundations
  gen <- rep(1:n_genotypes, each = sample_size*n_environments)
  env <- rep(1:n_environments, each = sample_size, times = n_environments) 
    
  # Random Noise
  noise <- rnorm(sample_size *n_genotypes * n_environments, 0, sd = std_dev) 
    
  # Interaction Terms
  int <- rep(rnorm(n_genotypes * n_environments, 0, sd = interaction),each = sample_size) # interaction term - one for each GE level
  
  # Create the model dataframe 
    model_df <- data.frame(gen, env, noise, int)
    model_df$gen_factor = factor(paste("G", model_df$gen, sep = "_"))
    model_df$exp_env_factor = factor(paste("E", model_df$env, sep = "_"))
    
    # Generate phenotype data using regression equation
    phen = delta_env * model_df$env + delta_gen * model_df$gen  + model_df$noise + model_df$int 
    model_df$phen = phen
    
    # Standardize data
    dat_avg <- mean(phen) 
    dat_std <- sd(phen)
    model_df$phen_corrected <- ((phen - dat_avg)/dat_std)
    
    # Anova
    test_temp <- aov(phen_corrected ~ exp_env_factor * gen_factor, data = model_df)
    
    # Estimated Marginal Means
    emm_options(msg.interaction = FALSE)
    emm_E = as.data.frame(emmeans(test_temp,"exp_env_factor"))
    emm_G = as.data.frame(emmeans(test_temp, "gen_factor"))
    emm_GxE = as.data.frame(emmeans(test_temp, ~ exp_env_factor*gen_factor))
    
    # Gmeans
    G_matrix = data.frame()
    for(h in 1:length(unique(emm_GxE$gen_factor))){
      gtemp <- filter(emm_GxE, gen_factor == unique(emm_GxE$gen_factor)[h])
      gmean <- sum(gtemp[,3])/n_genotypes
      tempdat = data.frame("G_means" = gmean,
                           "gen_factor" = unique(model_df$gen_factor)[h])
      G_matrix = rbind(G_matrix,tempdat)
    }
    
    # Emeans
    E_matrix = data.frame()
    for(j in 1:length(unique(emm_GxE$exp_env_factor))){
      etemp <- filter(emm_GxE, exp_env_factor == unique(emm_GxE$exp_env_factor)[j])
      emean <- sum(etemp[,3])/n_environments
      tempdat. = data.frame("E_means" = emean,
                            "exp_env_factor" = unique(model_df$exp_env_factor)[j])
      E_matrix = rbind(E_matrix,tempdat.)
    }
    
    # Covariance
    Cov_matrix = data.frame()
    Cov_matrix <- cbind(G_matrix,E_matrix)
    cov_est = cov(Cov_matrix$G_means,Cov_matrix$E_means)
    
    # Magnitude of GxE using EMMs
    GxE_emm <- abs(mean(model_df$phen_corrected) - # Overall mean
                     (emm_G$emmean[emm_G$gen_factor=="G_1"])- # G
                     (emm_E$emmean[emm_E$exp_env_factor=="E_1"])+ # E
                     (emm_GxE[1,3])) # GxE
    
    #############################
    ## True Covariance and GxE ##
    #############################
    
    # Generate phenotype data with no error using same regression equation
    no_err_phen = delta_env * model_df$env + delta_gen * model_df$gen + model_df$int 
    
    # Standardize no error data
    dat_avg_noerror <- mean(no_err_phen) 
    dat_std_noerror <- sd(no_err_phen)
    model_df$no_err_phen_corrected <- ((no_err_phen - dat_avg_noerror)/dat_std_noerror)
    
    # Anova with no error
    test_temp_noerror <- aov(no_err_phen_corrected ~ exp_env_factor * gen_factor, data = model_df)
    
    # Estimated Marginal Means with no error
    emm_options(msg.interaction = FALSE)
    emm_options(quietly = TRUE) # Will get perfect fit warning
    emm_E = as.data.frame(emmeans(test_temp_noerror,"exp_env_factor"))
    emm_G = as.data.frame(emmeans(test_temp_noerror, "gen_factor"))
    emm_GxE = as.data.frame(emmeans(test_temp_noerror, ~ exp_env_factor*gen_factor))
    
    # Gmeans with no error
    G_matrix = data.frame()
    for(y in 1:length(unique(emm_GxE$gen_factor))){
      gtemp <- filter(emm_GxE, gen_factor == unique(emm_GxE$gen_factor)[y])
      gmean <- sum(gtemp[,3])/n_genotypes
      tempdat = data.frame("G_means" = gmean,
                           "gen_factor" = unique(model_df$gen_factor)[y])
      G_matrix = rbind(G_matrix,tempdat)
    }
    
    # Emeans with no error
    E_matrix = data.frame()
    for(f in 1:length(unique(emm_GxE$exp_env_factor))){
      etemp <- filter(emm_GxE, exp_env_factor == unique(emm_GxE$exp_env_factor)[f])
      emean <- sum(etemp[,3])/n_environments
      tempdat. = data.frame("E_means" = emean,
                            "exp_env_factor" = unique(model_df$exp_env_factor)[f])
      E_matrix = rbind(E_matrix,tempdat.)
    }
    
    # Covariance with no error
    Cov_matrix = data.frame()
    Cov_matrix <- cbind(G_matrix,E_matrix)
    true_cov = cov(Cov_matrix$G_means,Cov_matrix$E_means)
    
    # Magnitude of GxE with no error 
    true_GxE <- abs(mean(model_df$phen_corrected) - # Overall mean
                      (emm_G$emmean[emm_G$gen_factor=="G_1"])- # G
                      (emm_E$emmean[emm_E$exp_env_factor=="E_1"])+ # E
                      (emm_GxE[1,3])) # GxE
    
    ###############
    ## Bootstrap ##
    ###############
    
    # Output Dataframe
    boot_df = data.frame()
    
    # Resampling Loop
    for(a in 1:n_boot){
      new_phen <- NULL
      shuffle_dat <- data.frame()
      
      # Each genotype and environment
      for (l in 1:nlevels(model_df$gen_factor)){
        for (j in 1:nlevels(model_df$exp_env_factor)){
          cond_G <- filter(model_df, gen_factor == unique(model_df$gen_factor)[l])
          cond_E <- filter(cond_G, exp_env_factor == unique(model_df$exp_env_factor)[j])
          
          # Shuffle data 
          new_phen <- sample(cond_E$phen_corrected, size=nrow(cond_E), replace=TRUE)
          
          # Output    
          shuffle_dat_temp <- data.frame("gen_factor" = cond_E$gen_factor,
                                         "exp_env_factor" = cond_E$exp_env_factor,
                                         "phen_corrected" = new_phen)
          shuffle_dat <- rbind(shuffle_dat, shuffle_dat_temp)
        }
      }
      
      # Anova
      test_boot <- aov(phen_corrected ~ exp_env_factor * gen_factor, data = shuffle_dat)
      
      # Estimated Marginal Means
      emm_options(msg.interaction = FALSE)
      emm_E_boot = as.data.frame(emmeans(test_boot,"exp_env_factor"))
      emm_G_boot = as.data.frame(emmeans(test_boot, "gen_factor"))
      emm_GxE_boot = as.data.frame(emmeans(test_boot, ~ exp_env_factor*gen_factor))
      
      # Gmeans - Bootstrap
      G_matrix_boot = data.frame()
      for(p in 1:length(unique(emm_GxE_boot$gen_factor))){
        g_boot <- filter(emm_GxE_boot, gen_factor == unique(emm_GxE_boot$gen_factor)[p])
        g_mean_boot <- sum(g_boot[,3])/length(unique(emm_GxE_boot$gen_factor))
        tempdat = data.frame("G_means" = g_mean_boot,
                             "gen_factor" = unique(g_boot$gen_factor))
        G_matrix_boot = rbind(G_matrix_boot,tempdat)
      }
      
      # Emeans - Bootstrap
      E_matrix_boot = data.frame()
      for(q in 1:length(unique(emm_GxE_boot$exp_env_factor))){
        e_boot <- filter(emm_GxE_boot, exp_env_factor == unique(emm_GxE_boot$exp_env_factor)[q])
        e_mean_boot <- sum(e_boot[,3])/length(unique(emm_GxE_boot$exp_env_factor))
        tempdat. = data.frame("E_means" = e_mean_boot,
                              "exp_env_factor" = unique(e_boot$exp_env_factor))
        E_matrix_boot = rbind(E_matrix_boot,tempdat.)
      }
      
      # Covariance - Bootstrap
      Cov_matrix_boot = data.frame()
      Cov_matrix_boot <- cbind(G_matrix_boot,E_matrix_boot)
      cov_est_boot = cov(Cov_matrix_boot$G_means,Cov_matrix_boot$E_means)
      
      # Magnitude of GxE - Bootstrap
      GxE_emm_boot <- abs(mean(shuffle_dat$phen_corrected) - # Overall mean
                            (emm_G_boot$emmean[emm_G_boot$gen_factor=="G_1"])- # G
                            (emm_E_boot$emmean[emm_E_boot$exp_env_factor=="E_1"])+ # E
                            (emm_GxE_boot[1,3])) # GxE
      
      boot_dat. <- data.frame("covariance" = cov_est_boot,
                              "GxE_mag" = GxE_emm_boot)
      
      boot_df <- rbind(boot_df,boot_dat.)
    }
    
    # Confidence Intervals
    GxE_CI = quantile(boot_df$GxE_mag, probs=c(0.025, 0.975), type=1) 
    GxE_avg = mean(boot_df$GxE_mag) 
    cov_CI = quantile(boot_df$covariance, probs=c(0.025, 0.975), type=1) 
    cov_avg = mean(boot_df$covariance) 
    
    #################
    ## Permutation ##
    #################
    
    # Output dataframe
    perm_df <- data.frame()
    
    # Shuffling loop
    for(b in 1:n_boot){
      
      # Shuffle data
      null_temp <- sample(model_df$phen_corrected, size=nrow(model_df), replace=FALSE)
      
      perm_dat <- data.frame("gen_factor" = model_df$gen_factor,
                             "exp_env_factor" = model_df$exp_env_factor,
                             "phen_corrected" = null_temp)
      
      # Anova
      test_perm <- aov(phen_corrected ~ exp_env_factor * gen_factor, data = perm_dat)
      
      # Estimated Marginal Means - Permutation
      emm_options(msg.interaction = FALSE)
      emm_E_perm = as.data.frame(emmeans(test_perm,"exp_env_factor"))
      emm_G_perm = as.data.frame(emmeans(test_perm, "gen_factor"))
      emm_GxE_perm = as.data.frame(emmeans(test_perm, ~ exp_env_factor*gen_factor))
      
      # Gmeans - Permutation
      G_matrix_perm = data.frame()
      for(r in 1:length(unique(emm_GxE_perm$gen_factor))){
        g_perm <- filter(emm_GxE_perm, gen_factor == unique(emm_GxE_perm$gen_factor)[r])
        g_mean_perm <- sum(g_perm[,3])/length(unique(emm_GxE_perm$gen_factor))
        tempdat = data.frame("G_means" = g_mean_perm,
                             "gen_factor" = unique(g_perm$gen_factor))
        G_matrix_perm = rbind(G_matrix_perm,tempdat)
      }
      
      # Emeans - Permutation
      E_matrix_perm = data.frame()
      for(s in 1:length(unique(emm_GxE_perm$exp_env_factor))){
        e_perm <- filter(emm_GxE_perm, exp_env_factor == unique(emm_GxE_perm$exp_env_factor)[s])
        e_mean_perm <- sum(e_perm[,3])/length(unique(emm_GxE_perm$exp_env_factor))
        tempdat. = data.frame("E_means" = e_mean_perm,
                              "exp_env_factor" = unique(e_perm$exp_env_factor))
        E_matrix_perm = rbind(E_matrix_perm,tempdat.)
      }
      
      # Covariance - Permutation
      Cov_matrix_perm = data.frame()
      Cov_matrix_perm <- cbind(G_matrix_perm,E_matrix_perm)
      cov_est_perm = cov(Cov_matrix_perm$G_means,Cov_matrix_perm$E_means)
      
      # Magnitude of GxE - Permutation
      GxE_emm_perm <- abs(mean(perm_dat$phen_corrected) - # Overall mean
                            (emm_G_perm$emmean[emm_G_perm$gen_factor=="G_1"])- # G
                            (emm_E_perm$emmean[emm_E_perm$exp_env_factor=="E_1"])+ # E
                            (emm_GxE_perm[1,3])) # GxE
      
      perm_dat. <- data.frame("covariance" = cov_est_perm,
                              "GxE_mag" = GxE_emm_perm)
      perm_df <- rbind(perm_df,perm_dat.)
    }
    
    # Covariance P-value
    ptemp = (rank(c(cov_est,perm_df[,1]))[1])/(n_boot+1) 
    cov_pvalue = NULL
    if(ptemp < 0.5){cov_pvalue = ptemp}else{cov_pvalue = (1-ptemp)} # 2-tailed
    
    # GxE P-value
    ptemp1 = (rank(c(GxE_emm,perm_df[,2]))[1])/(n_boot+1) 
    GxE_pvalue = NULL
    if(ptemp1 < 0.5){GxE_pvalue = ptemp1}else{GxE_pvalue = (1-ptemp1)} # 2-tailed
    
    # Generate Outputs
    output <- data.frame("row" = row,
                           "delta_env" = delta_env,
                           "delta_gen" = delta_gen,
                           "sample_size" = sample_size,
                           "n_genotypes" = n_genotypes,
                           "std_dev" = std_dev,
                           "interaction" = interaction,
                           "true_cov" = true_cov,
                           "cov_estimate" = cov_est,
                           "cov_lwrCI" = cov_CI[[1]],
                           "cov_uprCI" = cov_CI[[2]],
                           "cov_pvalue" = cov_pvalue,
                           "true_GxE" = true_GxE,
                           "GxE_estimate" = GxE_emm,
                           "GxE_lwrCI" = GxE_CI[[1]],
                           "GxE_uprCI" = GxE_CI[[2]],
                           "GxE_pvalue" = GxE_pvalue)


# Output
write.csv(output,paste0("Power_data_",row,"_output.csv"))

