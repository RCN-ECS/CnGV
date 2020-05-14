
# Function to analyze data from empirical work/ studies

amarillo_armadillo <- function(input_df, n_boot, data_type){ # Data, Number of bootstraps, data_type = c("means", "raw")
  
  # Load packages
  library("emmeans")
  library("lme4")
  library("tidyverse")
  
  #if(data_type == "raw"){
    
    # Output 
    output = data.frame()
    
    # Standardize data
    input_df$phen_corrected = (input_df$phen_data - mean(input_df$phen_data))/sd(input_df$phen_data)
    
    # Sanity Check 
    #ggplot(input_df,aes(x=exp_env_factor,y=phen_corrected, group = gen_factor, colour=gen_factor))+geom_point()+stat_smooth(method="glm")+theme_classic()
    
    # Anova
    test_temp <- lm(phen_corrected ~ exp_env_factor * gen_factor, data = input_df)
    
    # Estimated Marginal Means
    emm_options(msg.interaction = FALSE)
    emm_E = as.data.frame(emmeans(test_temp,"exp_env_factor"))
    emm_G = as.data.frame(emmeans(test_temp, "gen_factor"))
    emm_GxE = as.data.frame(emmeans(test_temp, ~ exp_env_factor*gen_factor))
    
    # Gmeans
    G_matrix = data.frame()
    for(h in 1:length(unique(emm_GxE$gen_factor))){
      gtemp <- filter(emm_GxE, gen_factor == unique(emm_GxE$gen_factor)[h])
      gmean <- mean(gtemp[,3])
      tempdat = data.frame("G_means" = gmean,
                           "gen_factor" = unique(emm_GxE$gen_factor)[h])
      G_matrix = rbind(G_matrix,tempdat)
    }
    
    # Emeans
    E_matrix = data.frame()
    for(j in 1:length(unique(emm_GxE$exp_env_factor))){
      etemp <- filter(emm_GxE, exp_env_factor == unique(emm_GxE$exp_env_factor)[j])
      emean <- mean(etemp[,3])
      tempdat. = data.frame("E_means" = emean,
                            "exp_env_factor" = unique(emm_GxE$exp_env_factor)[j])
      E_matrix = rbind(E_matrix,tempdat.)
    }
    
    # Match Genotypes to Native Environment
    G_matrix$exp_env_factor <- unique(input_df$nat_env_factor)[match(G_matrix$gen_factor,factor(unique(input_df$gen_factor)))]
    G_matrix$E_means <- E_matrix$E_means[match(E_matrix$exp_env_factor,G_matrix$exp_env_factor)]
    
    # Covariances
    correction_ratio = max(sd(G_matrix$E_means),sd(G_matrix$G_means))
    cov_corrected = round(cov(G_matrix$E_means, G_matrix$G_means)/(correction_ratio^2),2)
    
    # Magnitude of GxE -- EMMs
    GxE_emm <- abs(mean(input_df$phen_corrected) - # Overall mean
                     (emm_G$emmean[emm_G$gen_factor == "G_1"])- # G
                     (emm_E$emmean[emm_E$exp_env_factor == "E_1"])+ # E
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
      
      # Resample data within each genotype and environment
      for (l in 1:nlevels(input_df$gen_factor)){
        for (j in 1:nlevels(input_df$exp_env_factor)){
          cond_G <- filter(input_df, gen_factor == unique(input_df$gen_factor)[l])
          cond_E <- filter(cond_G, exp_env_factor == unique(input_df$exp_env_factor)[j])
          
          # Shuffle data 
          new_phen <- sample(cond_E$phen_corrected, size=nrow(cond_E), replace=TRUE)
          
          # Output    
          shuffle_dat_temp <- data.frame("gen_factor" = cond_E$gen_factor,
                                         "exp_env_factor" = cond_E$exp_env_factor,
                                         "phen_corrected" = new_phen)
          shuffle_dat <- rbind(shuffle_dat, shuffle_dat_temp)
        }
      }
      # Bootstrap Anova
      test_boot <- lm(phen_corrected ~ exp_env_factor * gen_factor, data = shuffle_dat)
      
      # Estimated Marginal Means
      emm_options(msg.interaction = FALSE)
      emm_E_boot = as.data.frame(emmeans(test_boot,"exp_env_factor"))
      emm_G_boot = as.data.frame(emmeans(test_boot, "gen_factor"))
      emm_GxE_boot = as.data.frame(emmeans(test_boot, ~ exp_env_factor*gen_factor))
      
      # Gmeans
      G_matrix_boot = data.frame()
      for(h in 1:length(unique(emm_GxE_boot$gen_factor))){
        gtemp <- filter(emm_GxE_boot, gen_factor == unique(emm_GxE_boot$gen_factor)[h])
        gmean <- mean(gtemp[,3])
        tempdat = data.frame("G_means" = gmean,
                             "gen_factor" = unique(emm_GxE$gen_factor)[h])
        G_matrix_boot = rbind(G_matrix_boot,tempdat)
      }
      
      # Emeans
      E_matrix_boot = data.frame()
      for(j in 1:length(unique(emm_GxE_boot$exp_env_factor))){
        etemp <- filter(emm_GxE_boot, exp_env_factor == unique(emm_GxE_boot$exp_env_factor)[j])
        emean <- mean(etemp[,3])
        tempdat. = data.frame("E_means" = emean,
                              "exp_env_factor" = unique(emm_GxE_boot$exp_env_factor)[j])
        E_matrix_boot = rbind(E_matrix_boot,tempdat.)
      }
      
      # Match Genotypes to Native Environment
      G_matrix_boot$exp_env_factor <- unique(input_df$nat_env_factor)[match(G_matrix_boot$gen_factor,factor(unique(input_df$gen_factor)))]
      G_matrix_boot$E_means <- E_matrix_boot$E_means[match(E_matrix_boot$exp_env_factor,G_matrix_boot$exp_env_factor)]
      
      # Covariances
      correction_ratio_boot = max(sd(G_matrix_boot$E_means),sd(G_matrix_boot$G_means))
      cov_corrected_boot = round(cov(G_matrix_boot$E_means, G_matrix_boot$G_means)/(correction_ratio_boot^2),2)
      
      # Magnitude of GxE -- EMMs
      GxE_emm_boot <- abs(mean(input_df$phen_corrected) - # Overall mean
                       (emm_G_boot$emmean[emm_G_boot$gen_factor == "G_1"])- # G
                       (emm_E_boot$emmean[emm_E_boot$exp_env_factor == "E_1"])+ # E
                       (emm_GxE_boot[1,3])) # GxE
      temp_boot <- data.frame("Covariance" = cov_corrected_boot,
                              "GxE" = round(GxE_emm_boot,2))
      boot_df <-  rbind(boot_df, temp_boot)
    } 
    
    # 95% Confidence Intervals
    cov_corrected_CI = quantile(boot_df$Covariance, probs=c(0.025, 0.975), type=1)
    GxE_CI = quantile(boot_df$GxE, probs=c(0.025, 0.975), type=1)
    
    #################
    ## Permutation ##
    #################
    
    # Output dataframe
    perm_df <- data.frame()
    
    # Shuffling loop
    for(b in 1:n_boot){
      
      # Shuffle data
      null_temp <- sample(input_df$phen_data, size=nrow(input_df), replace=FALSE)

      perm_dat <- data.frame("gen_factor" = input_df$gen_factor,
                             "exp_env_factor" = input_df$exp_env_factor,
                             "phen" = null_temp)
      # Re-Standardize 
      perm_dat$phen_corrected = (perm_dat$phen - mean(perm_dat$phen))/sd(perm_dat$phen)
      
      # Permutation Anova
      test_perm <- lm(phen_corrected ~ exp_env_factor * gen_factor, data = perm_dat)
    
      # Estimated Marginal Means -- Permutation
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
      
      # Match Genotypes to Native Environment
      G_matrix_perm$exp_env_factor <- unique(input_df$nat_env_factor)[match(G_matrix_perm$gen_factor,factor(unique(input_df$gen_factor)))]
      G_matrix_perm$E_means <- E_matrix_perm$E_means[match(E_matrix_perm$exp_env_factor,G_matrix$exp_env_factor)]
      
      # Covariance- permutation
      correction_ratio_perm = max(sd(G_matrix_perm$E_means),sd(G_matrix_perm$G_means))
      cov_corrected_perm = round(cov(G_matrix_perm$E_means, G_matrix_perm$G_means)/(correction_ratio_perm^2),2)
      
      # Magnitude of GxE -- EMMs -- Permutation
      GxE_emm_perm <- abs(mean(perm_dat$phen_corrected) - # Overall mean
                            (emm_G_perm$emmean[emm_G_perm$gen_factor=="G_1"])- # G
                            (emm_E_perm$emmean[emm_E_perm$exp_env_factor=="E_1"])+ # E
                            (emm_GxE_perm[1,3])) # GxE
      
      perm_dat <- data.frame("Covariance" = cov_corrected_perm,
                             "GxE" = GxE_emm_perm)
      perm_df <- rbind(perm_df,perm_dat)
    }
    
    # Covariance p-value
    ptemp1 = (rank(c(cov_corrected,perm_df$Covariance))[1])/(n_boot+1) 
    cov_pvalue = NULL
    if(ptemp1 < 0.5){cov_pvalue = ptemp1}else{cov_pvalue = (1-ptemp1)} # 2-tailed
    
    # GxE p-value
    ptemp2 = (rank(c(GxE_emm,perm_df$GxE))[1])/(n_boot+1) 
    GxE_pvalue = 1-ptemp2 # Right-tailed
    
    # Output
    output = data.frame("Covariance Estimate" = cov_corrected,
                        "Covariance Lower CI" = cov_corrected_CI[[1]],
                        "Covariance Upper CI" = cov_corrected_CI[[2]],
                        "Covariance p-value" = cov_pvalue,
                        "GxE Estimate" = GxE_emm,
                        "GxE Lower CI" = GxE_CI[[1]],
                        "GxE Upper CI" = GxE_CI[[2]],
                        "GxE p-value" = GxE_pvalue)
    return(output)
    
  #}else{ # ENTER MEANS CODE HERE
   #NULL} 
    }
  
testset = amarillo_armadillo(input_df, 100, raw)
