###########################################################################################
##              Functions for co/counter gradient meta-analysis                          ##
##        Authors: Molly Albecker, Thais Bittar, Geoff Trussell, Katie Lotterhos         ##
###########################################################################################


amarillo_armadillo <- function(input_df, n_boot, data_type){ # Data, Number of bootstraps, data_type = "raw" or "means"
  
  library(tibble)
  library(ggplot2)
  library(gridExtra)
  
  if(data_type == "raw"){
    
    # Output 
    output = data.frame()
    
    # Anova model fit & GxE estimates
    m1 <- mod.GxE(input_df, is.perm = FALSE) # Raw phenotype dataframe
    
    # GxE Estimates
    cov_matrix <- m1[[1]]
    GxE_emm <- m1[[3]]
    GxE_loop_output <- m1[[4]] # All GxE estimates from loop 
    omega2 <- m1[[5]]
    emm_df <- m1[[12]]
    
    # Covariance Estimates
    cov_corrected = round(cov.function(cov_matrix,emm_df),3)
    
    ###############
    ## Bootstrap ##
    ###############
    
    boot_dat_raw = boot_df_raw = data.frame()
    
    for(i in 1:n_boot){
      
      # Shuffle Data
      shuffle_dat <- bootstrap_raw(input_df) 
      
      # Anova model fit & GxE estimates
      m2 <- mod.GxE(shuffle_dat) # Insert shuffled raw phenotype dataframe
      
      # GxE Estimates
      cov_matrix_boot <- m2[[1]]
      GxE_emm_boot <- m2[[3]]
      GxE_loop_output_boot <- m2[[4]] # GxE output 
      omega2_boot <- m2[[5]]
      emm_df_boot <- m2[[12]]
      
      # Covariance Estimates
      cov_corrected_boot = round(cov.function(cov_matrix_boot,emm_df_boot),3)
      
      # Bootstrap dataframe
      boot_dat_raw <- data.frame("covariance_boot" = cov_corrected_boot,
                                 "GxE_emm_boot" = GxE_emm_boot,
                                 "GxE_omega_boot" = omega2_boot)
      boot_df_raw <- rbind(boot_df_raw,boot_dat_raw)
    }
    
    # Check: Histograms of Bootstrap
    # ggplot(boot_df_raw, aes(x = covariance_boot), alpha = 0.5)+
    #  geom_histogram()+ geom_vline(aes(xintercept = cov_corrected))+theme_classic()+ ggtitle("Bootstrap: Raw Data")
    
    
    # Covariance Confidence Intervals 
    cov_CI = quantile(boot_df_raw$covariance_boot, probs=c(0.025, 0.975), type=1) 
    
    # GxE Confidence Intervals
    GxE_emm_CI = quantile(boot_df_raw$GxE_emm_boot, probs = c(0.025, 0.975), type=1)
    GxE_omega_CI = quantile(boot_df_raw$GxE_omega_boot, probs=c(0.025, 0.975), type=1)
    
    
    #######################################
    #####   Permutation -- Raw Data   #####
    #######################################
    
    # Output dataframe
    perm_df_raw = perm_dat_raw = data.frame()
    
    for(i in 1:n_boot){
      
      # Resample Data
      perm_dat <- permutation_raw(input_df)
      
      # Anova model fit & GxE estimates
      m3 <- mod.GxE(perm_dat, is.perm = TRUE) # Insert permuted data
      
      # GxE Estimates
      cov_matrix_perm <- m3[[1]]
      GxE_emm_perm <- m3[[3]]
      GxE_loop_output_perm <- m3[[4]] # GxE output 
      omega2_perm <- m3[[5]]
      emm_df_perm <- m3[[12]]
      
      # Covariance Estimates
      cov_corrected_perm = round(cov.function(cov_matrix_perm,emm_df_perm),3)
      
      # Permutation dataframe
      perm_dat_raw <- data.frame("covariance_perm" = cov_corrected_perm,
                                 "GxE_emm_perm" = GxE_emm_perm,
                                 "GxE_omega_perm" = omega2_perm)
      perm_df_raw <- rbind(perm_df_raw,perm_dat_raw)
    }
    
    # Check: Permutation histogram
    # ggplot(perm_df_raw, aes(x = GxE_emm_perm), alpha = 0.5)+
    #  geom_histogram()+ geom_vline(aes(xintercept = GxE_emm))+theme_classic()+ ggtitle("Permutation: Raw Data")
    
    # Covariance P-values
    cov_pvalue <- pvalue_fun(cov_corrected,perm_df_raw$covariance_perm,"twotail", n_boot)
    
    # GxE P-values
    GxE_emm_pvalue <- pvalue_fun(GxE_emm,perm_df_raw$GxE_emm_perm,"righttail",n_boot)
    GxE_omega_pvalue <- pvalue_fun(omega2,perm_df_raw$GxE_omega_perm,"righttail",n_boot)
    
    # Output
    output = data.frame("Covariance Estimate" = cov_corrected,
                        "Covariance Lower CI" = cov_CI[[1]],
                        "Covariance Upper CI" = cov_CI[[2]],
                        "Covariance p-value" = cov_pvalue,
                        "GxE Estimate" = GxE_emm,
                        "GxE Lower CI" = GxE_emm_CI[[1]],
                        "GxE Upper CI" = GxE_emm_CI[[2]],
                        "GxE p-value" = GxE_emm_pvalue,
                        "Omega2" = omega2,
                        "Omega2 Lower CI" = GxE_omega_CI[[1]],
                        "Omega2 Upper CI" = GxE_omega_CI[[2]],
                        "Omega2 p-value" = GxE_omega_pvalue)
    return(output)
    
  }else{ 
    
    # Output 
    output = data.frame()
    
    # Plot (Should look very similar but shifted on Y-axis)
    # ggplot(input_df, aes(x = exp_env_factor, y = phen_data, group = gen_factor,colour=nat_env_factor)) + geom_point() + geom_line()
    # ggplot(input_df, aes(x = exp_env_factor, y = avg_phen_corrected, group = gen_factor,colour=nat_env_factor))+geom_point() + geom_line()
    
    # GxE estimates
    m4 <- mean.GxE(input_df) # Insert means data frame 
    
    # GxE 
    Cov_mean_matrix <- m4[[1]]
    GxE_means <- m4[[2]]
    GxE_means_loop_output <- m4[[3]]
    
    # Covariance
    cov_means = round(cov.function_means(Cov_mean_matrix,phen_df = input_df),3)
    
    ###################################
    ##### BOOTSTRAP -- MEANS DATA #####
    ###################################
    
    # Output Dataframes
    boot_df_means = boot_dat_means = data.frame()
    
    for(i in 1:n_boot){
      
      # Shuffle Data
      shuffle_means <- bootstrap_means(input_df) # Insert means data, Need n_boot seeds
      
      # GxE :: Covariance Matrix
      m5 <- mean.GxE(shuffle_means) # Insert shuffled up means data frame
      
      # GxE Estimates
      Cov_mean_matrix_boot <- m5[[1]]
      GxE_means_boot <- m5[[2]]
      
      # Covariance Estimates
      cov_corrected_mean_boot = round(cov.function_means(Cov_mean_matrix_boot, shuffle_means),3)
      
      # Bootstrap dataframe
      boot_dat_means <- data.frame("cov_means_boot" = cov_corrected_mean_boot,
                                   "GxE_means_boot" = GxE_means_boot)
      boot_df_means <- rbind(boot_df_means,boot_dat_means)
    }
    
    # Covariance Confidence Intervals -- Means
    cov_means_CI = quantile(boot_df_means$cov_means_boot, probs=c(0.025, 0.975), na.rm = TRUE, type=1) 
    
    # GxE Confidence Intervals -- Means
    GxE_means_CI = quantile(boot_df_means$GxE_means_boot, probs=c(0.025, 0.975), na.rm = TRUE, type=1) 
    
    #######################################
    #####  Permutation -- Means Data  #####
    #######################################
    
    # Output
    perm_df_means = perm_dat_means = data.frame()
    
    for(i in 1:n_boot){
      
      # Set seeds for perm_means and mean.GxE
      perm.seeds1 = i
      
      # Resample Data
      perm_means <- permutation_means(input_df,perm.seeds1) # Need a seed to ensure correct matching of genotype to native
      
      # GxE :: Covariance Matrix
      m6 <- mean.GxE(perm_means, is.perm = TRUE) # Insert resampled mean phenotype dataframe
      
      # Pull out Estimates
      Cov_mean_matrix_perm <- m6[[1]]
      GxE_means_perm <- m6[[2]]
      GxE_means_output_perm <- m6[[3]]
      
      # Covariance Estimates 
      cov_corrected_mean_perm = round(cov.function_means(Cov_mean_matrix_perm, perm_means),3)
      
      # Permutation dataframe -- Means
      perm_dat_means <- data.frame("cov_means_perm" = cov_corrected_mean_perm,
                                   "GxE_means_perm" = GxE_means_perm)
      perm_df_means <- rbind(perm_df_means,perm_dat_means)
    }
    
    # Check: Histogram
    # perm_df_means$line = cov_means
    # ggplot(perm_df_means, aes(x = cov_means_perm)) + geom_histogram() + geom_vline(aes(xintercept = line),colour = "red")
    
    # Covariance P-values
    cov_mean_pvalue <- pvalue_fun(cov_means,perm_df_means$cov_means_perm,"twotail", n_boot)
    
    # GxE P-values
    GxE_mean_pvalue <- pvalue_fun(GxE_means,perm_df_means$GxE_means_perm,"righttail",n_boot)
    
    # Output
    output = data.frame("Covariance Estimate" = cov_means,
                        "Covariance Lower CI" = cov_means_CI[[1]],
                        "Covariance Upper CI" = cov_means_CI[[2]],
                        "Covariance p-value" = cov_mean_pvalue,
                        "GxE Estimate" = GxE_means,
                        "GxE Lower CI" = GxE_means_CI[[1]],
                        "GxE Upper CI" = GxE_means_CI[[2]],
                        "GxE p-value" = GxE_mean_pvalue,
                        "Omega2" = NA,
                        "Omega2 Lower CI" = NA,
                        "Omega2 Upper CI" = NA,
                        "Omega2 p-value" = NA)
    return(output)
  } 
}

mod.GxE <- function(input_df,is.perm = FALSE){ # input is model_df
  
  # Outputs
  allGE <- c()
  loopGxE <- c()
  
  # Establish Native Environments
  native_df = data.frame("gen_factor" = unique(input_df$gen_factor))
  native_df$nat_env_factor = input_df$nat_env_factor[match(native_df$gen_factor,input_df$gen_factor)]
  
  # Anova
  aov.test <- lm(phen_corrected ~ exp_env_factor * gen_factor, data = input_df) 
  
  # Estimated Marginal Means
  emm_GxE = as.data.frame(emmeans(aov.test, ~ exp_env_factor*gen_factor))
  
  # Gmeans
  G_matrix <- data.frame("G_means" = tapply(emm_GxE$emmean, emm_GxE$gen_factor, mean, na.rm=TRUE),
                         "gen_factor" = unique(emm_GxE$gen_factor))
  E_matrix <- data.frame("E_means" = tapply(emm_GxE$emmean, emm_GxE$exp_env_factor, mean, na.rm=TRUE),
                         "exp_env_factor" = unique(emm_GxE$exp_env_factor))
  
  # Covariance 
  Cov_matrix = G_matrix
  Cov_matrix$exp_env_factor = native_df$nat_env_factor[match(G_matrix$gen_factor,native_df$gen_factor)] # Indicator Variable
  Cov_matrix$E_means = E_matrix$E_means[match(Cov_matrix$exp_env_factor,E_matrix$exp_env_factor)]
  
  # Output based on Stamps/Hadfield approach
  delta_E = ((emm_GxE$emmean[3]-emm_GxE$emmean[4])+(emm_GxE$emmean[1]-emm_GxE$emmean[2]))/2 
  delta_H = (emm_GxE$emmean[1]-emm_GxE$emmean[4])
  aov.coefs = coef(aov.test)
  
  # Omega^2
  w2_GxE = (summary(aov(aov.test))[[1]][3,2] - # (SS_effect -
              (summary(aov(aov.test))[[1]][3,1]*summary(aov(aov.test))[[1]][4,3])) / # (Df_effect * MS_error))/
              (sum(summary(aov(aov.test))[[1]][,2]) + # (SS_total+
              (summary(aov(aov.test))[[1]][4,3])) # MS_error)
  
  # Eta^2 
  eta_GxE = summary(aov(aov.test))[[1]][3,2]/sum(summary(aov(aov.test))[[1]][,2])
  
  # Proportion of GxE sums of squares without error
  GxE_SumsSquares = summary(aov(aov.test))[[1]][3,2]/sum(summary(aov(aov.test))[[1]][c(1:3),2])
  
  # Output model data
  mod_df <- as.data.frame(summary(aov(aov.test))[[1]])
  mod_df <- rownames_to_column(mod_df) 
  colnames(mod_df)[1] <- "Fixed_effect"
  
  if(is.perm == FALSE){ 
    
    # Magnitude of GxE -- EMMs
    GxE_emm_original<- abs(emm_GxE$emmean[emm_GxE$gen_factor == "G_1" & emm_GxE$exp_env_factor == "E_1"] - # GxE (Phenotype of ith genotype in jth environment)
                             mean(emm_GxE$emmean[emm_GxE$gen_factor == "G_1"]) - # phenotype of ith Genotype
                             mean(emm_GxE$emmean[emm_GxE$exp_env_factor == "E_1"]) + # phenotype of jth Environment
                             mean(emm_GxE$emmean)) # Overall mean phenotype
    
    # Magnitude of GxE -- Loop
    allGE = c()
    loopGxE = NULL
    for (i in 1:nlevels(emm_GxE$gen_factor)){
      for (j in 1:nlevels(emm_GxE$exp_env_factor)){
        G_levels <- levels(emm_GxE$gen_factor)
        E_levels <- levels(emm_GxE$exp_env_factor)
        loopGxE <- abs(emm_GxE$emmean[emm_GxE$gen_factor == G_levels[i] & emm_GxE$exp_env_factor == E_levels[j]] - # GxE (Phenotype of ith genotype in jth environment)
                         mean(emm_GxE$emmean[emm_GxE$gen_factor == G_levels[i]],na.rm=TRUE) - # mean phenotype of ith Genotype
                         mean(emm_GxE$emmean[emm_GxE$exp_env_factor == E_levels[j]],na.rm=TRUE) + # mean phenotype of jth Environment
                         mean(emm_GxE$emmean,na.rm=TRUE)) # Overall mean
        allGE <- c(allGE, loopGxE)
      }
    }
    
    GxE_emm_loop = mean(allGE,na.rm=TRUE)
    
  }else{ # Below generates null distribution for null of G+E means
    
    GxE_emm_original<- NULL # Placeholder, may delete
    
    allGE <- NULL
    for (i in 1:nlevels(input_df$gen_factor)){
      for (j in 1:nlevels(input_df$exp_env_factor)){
        
        G_levels <- levels(input_df$gen_factor)
        E_levels <- levels(input_df$exp_env_factor)
        
        Gi_mean <- mean(sample(input_df$phen_corrected[input_df$gen_factor == G_levels[i]], 
                               size = length(input_df$phen_corrected[input_df$gen_factor == G_levels[i]]),
                               replace = TRUE), na.rm=TRUE)
        Ej_mean <- mean(sample(input_df$phen_corrected[input_df$exp_env_factor == E_levels[j]],
                               size = length(input_df$phen_corrected[input_df$exp_env_factor == E_levels[j]]),
                               replace = TRUE), na.rm=TRUE)
        GE_sd <- input_df %>%
          filter(gen_factor == G_levels[i]) %>%
          filter(exp_env_factor == E_levels[j]) %>%
          summarize("GEsd" = sd(phen_corrected, na.rm = TRUE))
        
        # Create a sample of the null expectation for the Gi+Ej
        GiEj_null_samp <- rnorm(1, mean = (Gi_mean + Ej_mean), sd = abs(GE_sd[[1]]))
        
        # Estimate 
        GxE_mean.temp <- abs(GiEj_null_samp - # G+E (Phenotype of ith genotype in jth environment)
                               Gi_mean - # mean phenotype of ith Genotype
                               Ej_mean + # mean phenotype of jth Environment
                               mean(input_df$phen_corrected, na.rm=TRUE)) # Overall mean
        allGE <- c(allGE, GxE_mean.temp)
      }
    }
    
    GxE_emm_loop = mean(allGE,na.rm = TRUE)
    
  }
  return(list(Cov_matrix, GxE_emm_original, GxE_emm_loop, allGE, w2_GxE, eta_GxE, GxE_SumsSquares, mod_df, delta_E, delta_H, aov.coefs,emm_GxE))
}

mean.GxE <- function(input_df,is.perm = FALSE){ # input is mean_df
  
  # Clear outputs
  allGEmeans <- c()
  GxE_mean.temp <- c()
  GiEj_mean = Gi_mean = Ej_mean = GiEj_null_samp = NULL
  
  if(is.perm == FALSE){ 
    
    # Means of Means
    E_means <- tapply(input_df$avg_phen_corrected, input_df$exp_env_factor, mean, na.rm=TRUE)
    G_means <- tapply(input_df$avg_phen_corrected, input_df$gen_factor, mean, na.rm=TRUE)
    Gmean_mat <- data.frame("G_means" = G_means, "gen_factor" = unique(input_df$gen_factor))
    Emean_mat <- data.frame("E_means" = E_means, "exp_env_factor" = unique(input_df$exp_env_factor))
    
    # Match means to native
    Cov_mean_matrix = Gmean_mat
    Cov_mean_matrix$exp_env_factor <- input_df$nat_env_factor[match(Cov_mean_matrix$gen_factor,input_df$gen_factor)]
    Cov_mean_matrix$E_means <- Emean_mat$E_means[match(Cov_mean_matrix$exp_env_factor,Emean_mat$exp_env_factor)]
    
    # Magnitude of GxE -- Loop -- Means
    for (p in 1:nlevels(input_df$gen_factor)){
      for (q in 1:nlevels(input_df$exp_env_factor)){
        G_levels <- levels(input_df$gen_factor)
        E_levels <- levels(input_df$exp_env_factor)
        GxE_mean.temp <- abs(input_df$avg_phen_corrected[input_df$gen_factor == G_levels[p] & input_df$exp_env_factor == E_levels[q]] - # GxE (Phenotype of ith genotype in jth environment)
                               mean(input_df$avg_phen_corrected[input_df$gen_factor == G_levels[p]], na.rm=TRUE)- # mean phenotype of ith Genotype
                               mean(input_df$avg_phen_corrected[input_df$exp_env_factor == E_levels[q]], na.rm=TRUE)+ # mean phenotype of jth Environment
                               mean(input_df$avg_phen_corrected, na.rm=TRUE)) # Overall mean
        allGEmeans <- c(allGEmeans, GxE_mean.temp)
      }
    }
    
    
    GxE_means = mean(allGEmeans, na.rm=TRUE)
    
  }else{ # Below generates null distribution for null of G+E means
    
    # Means of Means
    E_means <- tapply(input_df$avg_phen_corrected, input_df$exp_env_factor, mean, na.rm=TRUE)
    G_means <- tapply(input_df$avg_phen_corrected, input_df$gen_factor, mean, na.rm=TRUE)
    Gmean_mat <- data.frame("G_means" = G_means, "gen_factor" = unique(input_df$gen_factor))
    Emean_mat <- data.frame("E_means" = E_means, "exp_env_factor" = unique(input_df$exp_env_factor))
    
    # Match means to native
    Cov_mean_matrix = Gmean_mat
    Cov_mean_matrix$exp_env_factor <- input_df$nat_env_factor[match(Cov_mean_matrix$gen_factor,input_df$gen_factor)]
    Cov_mean_matrix$E_means <- Emean_mat$E_means[match(Cov_mean_matrix$exp_env_factor,Emean_mat$exp_env_factor)]
    
    for (r in 1:nlevels(input_df$gen_factor)){
      for (s in 1:nlevels(input_df$exp_env_factor)){
        
        G_levels <- levels(input_df$gen_factor)
        E_levels <- levels(input_df$exp_env_factor)
        
        GiEj_mean <- input_df$avg_phen_corrected[input_df$gen_factor == G_levels[r] & input_df$exp_env_factor == E_levels[s]]
        Gi_mean <- mean(input_df$avg_phen_corrected[input_df$gen_factor == G_levels[r]], na.rm=TRUE)
        Ej_mean <- mean(input_df$avg_phen_corrected[input_df$exp_env_factor == E_levels[s]], na.rm=TRUE)
        
        # Create a sample of the null expectation for the GiEj
        GiEj_null_samp <- rnorm(1, mean = (Gi_mean + Ej_mean), sd = mean(input_df$error[input_df$gen_factor == G_levels[r] & input_df$exp_env_factor == E_levels[s]], na.rm=TRUE))
        
        # Estimate 
        GxE_mean.temp <- abs(GiEj_null_samp - # GxE (Phenotype of ith genotype in jth environment)
                               Gi_mean - # mean phenotype of ith Genotype
                               Ej_mean + # mean phenotype of jth Environment
                               mean(input_df$avg_phen_corrected, na.rm=TRUE)) # Overall mean
        allGEmeans <- c(allGEmeans, GxE_mean.temp)
      }
    }
    
    GxE_means = mean(allGEmeans, na.rm=TRUE)
    
  }
  
  return(list(Cov_mean_matrix, GxE_means, allGEmeans))
}

bootstrap_raw <- function(input_df){ # input is model_df
  
  # Clear outputs
  new_phen <- NULL
  shuffle_dat <- data.frame()
  shuffle_dat_temp <- data.frame()
  
  # Resample data within each genotype and environment
  for (l in 1:nlevels(input_df$gen_factor)){
    for (j in 1:nlevels(input_df$exp_env_factor)){
      
      cond <- input_df %>%
        filter(gen_factor == unique(input_df$gen_factor)[l]) %>%
        filter(exp_env_factor == unique(input_df$exp_env_factor)[j])
      
      # Shuffle data 
      new_phen <- sample(cond$phen_corrected, size=nrow(cond), replace=TRUE)
      
      # Output    
      shuffle_dat_temp <- data.frame("gen_factor" = cond$gen_factor,
                                     "exp_env_factor" = cond$exp_env_factor,
                                     "nat_env_factor" = cond$nat_env_factor,
                                     "phen_corrected" = new_phen)
      shuffle_dat <- rbind(shuffle_dat, shuffle_dat_temp)
    }
  }
  return(shuffle_dat)
}

bootstrap_means <- function(input_df){ # input is means_df
  
  # Clear outputs
  new_phen.<- new_phen <- NULL
  new_mean_temp <- data.frame()
  new_means <- data.frame()
  
  # Resample means dataframe
  for (u in 1:nlevels(input_df$gen_factor)){
    for (r in 1:nlevels(input_df$exp_env_factor)){
      
      # Retain levels
      cond <- input_df %>%
        filter(gen_factor == unique(input_df$gen_factor)[u]) %>%
        filter(exp_env_factor == unique(input_df$exp_env_factor)[r])
      
      # Create new means data
      new_phen <- rnorm(nrow(cond), mean = cond$avg_phen_corrected, sd = cond$error)
      
      # Output
      new_mean_temp <- data.frame("gen_factor" = cond$gen_factor,
                                  "exp_env_factor" = cond$exp_env_factor,
                                  "nat_env_factor" = cond$nat_env_factor,
                                  "avg_phen_corrected" = new_phen)
      new_means <- rbind(new_means, new_mean_temp)
    }
  }
  
  return(new_means)
}

permutation_raw <- function(input_df){ # input is model_df
  
  # Clear outputs
  perm_dat = data.frame()
  null_temp <- NULL
  
  # Shuffle raw data
  null_temp <- sample(input_df$phen_corrected, size=nrow(input_df), replace=FALSE)
  
  perm_dat <- data.frame("gen_factor" = input_df$gen_factor,
                         "exp_env_factor" = input_df$exp_env_factor,
                         "nat_env_factor" = input_df$nat_env_factor,
                         "phen_corrected" = null_temp)
  return(perm_dat)
}

permutation_means <- function(input_df, perm.seed){ # means dataframe (mean_df)
  
  # Clear outputs
  perm_means <- data.frame()
  null_gen = null_env = null_means = NULL
  
  # Shuffle means data (same set.seed keeps phen and corresponding se matched)
  set.seed(perm.seed)
  null_means <- sample(input_df$avg_phen_corrected, size = length(input_df$avg_phen_corrected), replace = FALSE)
  
  # Create new means data
  set.seed(perm.seed)
  null_se <- sample(input_df$error, size = length(input_df$error), replace = FALSE)

  perm_means <- data.frame("gen_factor" = input_df$gen_factor,
                           "exp_env_factor" = input_df$exp_env_factor,
                           "nat_env_factor" = input_df$nat_env_factor,
                           "avg_phen_corrected" = null_means,
                           "error" = null_se)
 
  return(perm_means)
}

pvalue_fun <- function(estimate, rankdat, test, n_boot){ #Test = "twotail" or "righttail"
  
  p.value = NULL
  
  if(test == "twotail"){
    p.value = sum(abs(rankdat) >= abs(estimate))/(n_boot+1) # Two-tailed
  }else{
    p.value = sum(rankdat >= estimate)/(n_boot+1) # Right-tailed
  }
  return(p.value)
}

cov.function <- function(cov_mat, emm_df, is.sample = TRUE){ # input_df = cov matrix
  
  pvar <- function(x) { #population variance - for sample variance can use var() in R
    sum((x - mean(x))**2) / length(x)
  }
  
  N = length(cov_mat$gen_factor)
  overallmean = mean(emm_df$emmean,na.rm=TRUE) 
  numerator = sum((cov_mat$G_means - overallmean)*(cov_mat$E_means - overallmean))
  
  if(is.sample == TRUE){
    sample_correcter = max(var(cov_mat$E_means),var(cov_mat$G_means))
    cv = (1/(N-1))*(numerator/sample_correcter)
  }else{
    population_correcter = max(pvar(cov_mat$E_means),pvar(cov_mat$G_means))
    cv = (1/(N))*(numerator/population_correcter)
  }
  return(cv)
} 

cov.function_means <- function(cov_mat, phen_df, is.sample = TRUE){ # input_df = cov_matrix of G_means and E_means
  
  pvar <- function(x) { #population variance - for sample variance can use var() in R
    sum((x - mean(x))**2) / length(x)
  }
  
  N = length(cov_mat$gen_factor)
  
  overallmean = mean(phen_df$avg_phen_corrected,na.rm = TRUE)
  numerator = sum((cov_mat$G_means - overallmean)*(cov_mat$E_means - overallmean))
  
  if(is.sample == TRUE){
    sample_correcter = max(var(cov_mat$E_means),var(cov_mat$G_means))
    cv = (1/(N-1))*(numerator/sample_correcter)
  }else{
    population_correcter = max(pvar(cov_mat$E_means),pvar(cov_mat$G_means))
    cv = (1/(N))*(numerator/population_correcter)
  }
  return(cv)
}

empty_as_na <- function(x){
  if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
  ifelse(as.character(x)!="", x, NA)
}

is.empty <- function(x, mode=NULL){
  if (is.null(mode)) mode <- class(x)
  identical(vector(mode,1),c(x,vector(class(x),1)))
}
