# Starting list of parameters
param_list <- list(
  reps = 7,
  delta_env = c(0.01),#,1), # the amount the phenotype changes across 1 value of the environment (i.e., the slope). This is essentially the amount/degree of phenotypic plasticity that is the same across genotypes.
  delta_gen = c(-1),#,0,1), # the amount the phenotype changes from one genotype to the next. This is essitially the increase intercept from one genotype to the next.
  sample_size = c(5), 
  n_genotypes = c(2),
  n_environments = NULL,
  std_dev= c(0.5), # Random noise, with standard deviation of 1,
  interaction= c(0,20)) # this sd determines the amount of GxE)


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
                            param_temp)
  
  return(param_table)
}
df = table_fun(param_list)
dim(df)
write.csv(df, file = "~/Desktop/df.csv")


ring <- function(param_table, n_boot){
  
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
    
    # For reproducibility
    set.seed = 999
    
    n_environments = param_table$n_genotypes[i]
    
    # Approximate Cov(G,E)
    cov_GE_approx = param_table$delta_env[i] * param_table$delta_gen[i]
    
    gen <- rep(1:param_table$n_genotypes[i], each = param_table$sample_size[i]*n_environments)
    env <- rep(1:n_environments, each = param_table$sample_size[i],times = n_environments) 
    
    noise <- rnorm(param_table$sample_size[i] * param_table$n_genotypes[i] * n_environments, 0, sd = param_table$std_dev[i]) # Random noise
    
    # Create Interactions
    int <- rnorm(param_table$n_genotypes[i] * n_environments * param_table$sample_size[i], 0, sd = param_table$interaction[i]) # sd determines the amount of GxE
    #int_df <- data.frame(expand.grid(G = 1:param_table$n_genotypes[i], E = 1:n_environments), int)
    
    # Create the model dataframe 
    model_df <- data.frame(gen, env, noise, int)
    #model_df <- merge(model_df, int_df)
    model_df$gen_factor = factor(paste("G", model_df$gen, sep = "_"))
    model_df$exp_env_factor = factor(paste("E", model_df$env, sep = "_"))
    
    # Generate phenotype data using regression equation
    phen = param_table$delta_env[i] * model_df$env + param_table$delta_gen[i] * model_df$gen  + model_df$noise + model_df$int
    no_err_phen = param_table$delta_env[i] * model_df$env + param_table$delta_gen[i] * model_df$gen  + 0 + model_df$int
    model_df$phen = phen
    
    # Standardize data
    dat_avg <- mean(phen) 
    dat_std <- sd(phen)
    model_df$phen_corrected <- ((phen - dat_avg)/dat_std)
    
      # New dataset from means and SE
    mean_df = model_df %>%
        group_by(gen_factor, exp_env_factor) %>%
        summarise(new_mean = mean(phen_corrected),
                  new_sd = sd(phen_corrected))
    mean_df$se = mean_df$new_sd/sqrt(unique(model_df$sample_size))

    # Marginal Means
    overall_mean <- mean(mean_df$new_mean)
    
    mm_df = data.frame()
    for(a in 1:nrow(mean_df)){
      name = paste(mean_df$gen_factor[a],mean_df$exp_env_factor[a],sep = ":") 
      value = mean_df$new_mean[a] - overall_mean
      mm_df. = data.frame(name, value)
      mm_df = rbind(mm_df, mm_df.)
    }
    
    G11m <- mean_df$new_mean[1] - overall_mean
    G12m <- mean_df$new_mean[2] - overall_mean
    G21m <- mean_df$new_mean[3] - overall_mean
    G22m <- mean_df$new_mean[4] - overall_mean

    G1_mean <- mean(c(G11m,G12m)) ## THIS IS WHERE YOU ARE = NEED TO SOFT CODE THIS IN
    E1_mean <- mean(c(G11m,G21m))
    gxe = abs(overall_mean - G1_mean - E1_mean + mean_df$new_mean[1])
  
    # Gmeans and Emeans  
    G_matrix = data.frame()
    E_matrix = data.frame()
    Cov_matrix = data.frame()
    
    for(h in 1:length(unique(emm_GxE$gen_factor))){
      gtemp <- filter(emm_GxE, gen_factor == unique(emm_GxE$gen_factor)[h])
      gmean <- sum(gtemp[,3])/param_table$n_genotypes[h]
      tempdat = data.frame("G_means" = gmean,
                           "gen_factor" = unique(model_df$gen_factor)[h])
      G_matrix = rbind(G_matrix,tempdat)
    }
    
    for(j in 1:length(unique(emm_GxE$exp_env_factor))){
      etemp <- filter(emm_GxE, exp_env_factor == unique(emm_GxE$exp_env_factor)[j])
      emean <- sum(etemp[,3])/n_environments
      tempdat. = data.frame("E_means" = emean,
                            "exp_env_factor" = unique(model_df$exp_env_factor)[j])
      E_matrix = rbind(E_matrix,tempdat.)
    }
    
    # Covariance
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
    
    # Generate phenotype data using regression equation
    no_err_phen = param_table$delta_env[i] * model_df$env + param_table$delta_gen[i] * model_df$gen + model_df$int # no error
    
    # Standardize no error data
    dat_avg_noerror <- mean(no_err_phen) 
    dat_std_noerror <- sd(no_err_phen)
    model_df$no_err_phen_corrected <- ((no_err_phen - dat_avg_noerror)/dat_std_noerror)
    
    # Anova
    test_temp_noerror <- aov(no_err_phen_corrected ~ exp_env_factor * gen_factor, data = model_df)
    
    # Estimated Marginal Means
    emm_options(msg.interaction = FALSE)
    emm_E = as.data.frame(emmeans(test_temp_noerror,"exp_env_factor"))
    emm_G = as.data.frame(emmeans(test_temp_noerror, "gen_factor"))
    emm_GxE = as.data.frame(emmeans(test_temp_noerror, ~ exp_env_factor*gen_factor))
    
    # Gmeans and Emeans  
    G_matrix = data.frame()
    E_matrix = data.frame()
    Cov_matrix = data.frame()
    
    for(y in 1:length(unique(emm_GxE$gen_factor))){
      gtemp <- filter(emm_GxE, gen_factor == unique(emm_GxE$gen_factor)[y])
      gmean <- sum(gtemp[,3])/param_table$n_genotypes[y]
      tempdat = data.frame("G_means" = gmean,
                           "gen_factor" = unique(model_df$gen_factor)[y])
      G_matrix = rbind(G_matrix,tempdat)
    }
    
    for(f in 1:length(unique(emm_GxE$exp_env_factor))){
      etemp <- filter(emm_GxE, exp_env_factor == unique(emm_GxE$exp_env_factor)[f])
      emean <- sum(etemp[,3])/n_environments
      tempdat. = data.frame("E_means" = emean,
                            "exp_env_factor" = unique(model_df$exp_env_factor)[f])
      E_matrix = rbind(E_matrix,tempdat.)
    }
    
    # Covariance
    Cov_matrix <- cbind(G_matrix,E_matrix)
    true_cov = cov(Cov_matrix$G_means,Cov_matrix$E_means)
    
    # Magnitude of GxE using EMMs
    true_GxE <- abs(mean(model_df$phen_corrected) - # Overall mean
                      (emm_G$emmean[emm_G$gen_factor=="G_1"])- # G
                      (emm_E$emmean[emm_E$exp_env_factor=="E_1"])+ # E
                      (emm_GxE[1,3])) # GxE
    
    ###############
    ## Bootstrap ##
    ###############
    boot_df = data.frame()
    
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
      
      # Gmeans and Emeans  
      G_matrix_boot = data.frame()
      E_matrix_boot = data.frame()
      Cov_matrix_boot = data.frame()
      
      for(p in 1:length(unique(emm_GxE_boot$gen_factor))){
        g_boot <- filter(emm_GxE_boot, gen_factor == unique(emm_GxE_boot$gen_factor)[p])
        g_mean_boot <- sum(g_boot[,3])/length(unique(emm_GxE_boot$gen_factor))
        tempdat = data.frame("G_means" = g_mean_boot,
                             "gen_factor" = unique(g_boot$gen_factor))
        G_matrix_boot = rbind(G_matrix_boot,tempdat)
      }
      
      for(q in 1:length(unique(emm_GxE_boot$exp_env_factor))){
        e_boot <- filter(emm_GxE_boot, exp_env_factor == unique(emm_GxE_boot$exp_env_factor)[q])
        e_mean_boot <- sum(e_boot[,3])/length(unique(emm_GxE_boot$exp_env_factor))
        tempdat. = data.frame("E_means" = e_mean_boot,
                              "exp_env_factor" = unique(e_boot$exp_env_factor))
        E_matrix_boot = rbind(E_matrix_boot,tempdat.)
      }
      
      # Covariance
      Cov_matrix_boot <- cbind(G_matrix_boot,E_matrix_boot)
      cov_est_boot = cov(Cov_matrix_boot$G_means,Cov_matrix_boot$E_means)
      
      # Magnitude of GxE using EMMs
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
    
    perm_df <- data.frame()
    
    for(b in 1:n_boot){
      
      # Shuffle data
      null_temp <- sample(model_df$phen_corrected, size=nrow(model_df), replace=FALSE)
      
      perm_dat <- data.frame("gen_factor" = model_df$gen_factor,
                             "exp_env_factor" = model_df$exp_env_factor,
                             "phen_corrected" = null_temp)
      
      # Anova
      test_perm <- aov(phen_corrected ~ exp_env_factor * gen_factor, data = perm_dat)
      
      # Estimated Marginal Means
      emm_options(msg.interaction = FALSE)
      emm_E_perm = as.data.frame(emmeans(test_perm,"exp_env_factor"))
      emm_G_perm = as.data.frame(emmeans(test_perm, "gen_factor"))
      emm_GxE_perm = as.data.frame(emmeans(test_perm, ~ exp_env_factor*gen_factor))
      
      # Gmeans and Emeans  
      G_matrix_perm = data.frame()
      E_matrix_perm = data.frame()
      Cov_matrix_perm = data.frame()
      
      for(r in 1:length(unique(emm_GxE_perm$gen_factor))){
        g_perm <- filter(emm_GxE_perm, gen_factor == unique(emm_GxE_perm$gen_factor)[r])
        g_mean_perm <- sum(g_perm[,3])/length(unique(emm_GxE_perm$gen_factor))
        tempdat = data.frame("G_means" = g_mean_perm,
                             "gen_factor" = unique(g_perm$gen_factor))
        G_matrix_perm = rbind(G_matrix_perm,tempdat)
      }
      
      for(s in 1:length(unique(emm_GxE_perm$exp_env_factor))){
        e_perm <- filter(emm_GxE_perm, exp_env_factor == unique(emm_GxE_perm$exp_env_factor)[s])
        e_mean_perm <- sum(e_perm[,3])/length(unique(emm_GxE_perm$exp_env_factor))
        tempdat. = data.frame("E_means" = e_mean_perm,
                              "exp_env_factor" = unique(e_perm$exp_env_factor))
        E_matrix_perm = rbind(E_matrix_perm,tempdat.)
      }
      
      # Covariance
      Cov_matrix_perm <- cbind(G_matrix_perm,E_matrix_perm)
      cov_est_perm = cov(Cov_matrix_perm$G_means,Cov_matrix_perm$E_means)
      
      # Magnitude of GxE using EMMs
      GxE_emm_perm <- abs(mean(perm_dat$phen_corrected) - # Overall mean
                            (emm_G_perm$emmean[emm_G_perm$gen_factor=="G_1"])- # G
                            (emm_E_perm$emmean[emm_E_perm$exp_env_factor=="E_1"])+ # E
                            (emm_GxE_perm[1,3])) # GxE
      
      perm_dat. <- data.frame("covariance" = cov_est_perm,
                              "GxE_mag" = GxE_emm_perm)
      perm_df <- rbind(perm_df,perm_dat.)
    }
    
    # Covariance p-value
    ptemp = (rank(c(cov_est,perm_df[,1]))[1])/(n_boot+1) 
    cov_pvalue = NULL
    if(ptemp < 0.5){cov_pvalue = ptemp}else{cov_pvalue = (1-ptemp)}
    
    # GxE p-value
    ptemp1 = (rank(c(GxE_emm,perm_df[,2]))[1])/(n_boot+1) 
    GxE_pvalue = NULL
    if(ptemp1 < 0.5){GxE_pvalue = ptemp1}else{GxE_pvalue = (1-ptemp1)}
    
    # Generate Outputs
    temp_out <- data.frame("row" = param_table$row[i],
                           "delta_env" = param_table$delta_env[i],
                           "delta_gen" = param_table$delta_gen[i],
                           "sample_size" = param_table$sample_size[i],
                           "n_genotypes" = param_table$n_genotypes[i],
                           "std_dev" = param_table$std_dev[i],
                           "interaction" = param_table$interaction[i],
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
    output = rbind(output, temp_out)
  }
  return(output)
}

test = ring(df,25) # Parameter table, then number of bootstraps/perms    
