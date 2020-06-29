###################################################################################
##              Functions for co/counter gradient simulations                    ##
##        Authors: Molly Albecker, Geoff Trussell, Katie Lotterhos               ##
###################################################################################

df.foundations <- function(delta_env, delta_gen, sample_size, n_env, std_dev, n_pop, interaction){
  
  # Dataframe foundations
  n_environments <- n_env 
  gen <- rep(1:n_pop, each = sample_size*n_environments)
  env <- rep(1:n_environments, times = n_pop, each = sample_size) 
  
  # Random Noise
  noise <- rnorm(sample_size * n_pop * n_environments, 0, sd = std_dev) 
  
  # Interaction Term
  int <- rep(rnorm(n_pop * n_environments, 0, sd = interaction), each = sample_size) # interaction term - one for each GE level
  
  ############################################
  ## Phenotype data with standard deviation ##
  ############################################
  
  # Create model dataframe 
  model_df <- data.frame(gen, env, noise, int)
  model_df$gen_factor <- factor(paste("G", model_df$gen, sep = "_"))
  model_df$exp_env_factor <- factor(paste("E", model_df$env, sep = "_"))
  
  # Native environments
  model_df$nat_env_factor = rep(unique(model_df$exp_env_factor), each = (sample_size*n_pop))
  
  # Phenotype data using regression equation
  phen = delta_env * model_df$env + delta_gen * model_df$gen  + model_df$noise + model_df$int 
  model_df$phen = phen
  
  # Generate means dataframe
  mean_df = model_df %>%
    group_by(gen_factor, exp_env_factor) %>%
    summarize("avg_phen" = mean(phen),
              "deviation" = sd(phen))
  mean_df$se <- mean_df$deviation/(sqrt(sample_size))
  mean_df$nat_env_factor <- model_df$nat_env_factor[match(mean_df$gen_factor,model_df$gen_factor)]
  
  # Standardize data
  model_df$phen_corrected <- (phen - mean(phen))/sd(phen)
  mean_df$avg_phen_corrected <- (mean_df$avg_phen - mean(mean_df$avg_phen))/sd(mean_df$avg_phen) 
  
  ###############################################
  ## Phenotype data with no standard deviation ##
  ###############################################
  
  # Create model dataframe 
  df.ne <- data.frame(gen, env, int)
  df.ne$gen_factor <- factor(paste("G", df.ne$gen, sep = "_"))
  df.ne$exp_env_factor <- factor(paste("E", df.ne$env, sep = "_"))
  
  # Native environments
  df.ne$nat_env_factor = rep(unique(df.ne$exp_env_factor), each = (sample_size*n_pop))
  
  # Simulate phenotype
  no_err_phen = (delta_env * df.ne$env) + (delta_gen * df.ne$gen) + df.ne$int
  df.ne$no_err_phen = no_err_phen
  
  # Generate means dataframe
  mean_df_ne = df.ne %>%
    group_by(gen_factor, exp_env_factor) %>%
    summarize(avg_phen_ne = mean(no_err_phen),
              deviation_ne = sd(no_err_phen))
  mean_df_ne$se_ne = mean_df_ne$deviation_ne/(sqrt(sample_size)) # Should be zero 
  
  # Add native environments
  mean_df_ne$nat_env_factor = df.ne$nat_env_factor[match(mean_df_ne$gen_factor,df.ne$gen_factor)]
  
  # Standardize data
  df.ne$phen_corrected = (no_err_phen - mean(no_err_phen))/sd(no_err_phen)
  mean_df_ne$avg_phen_corrected = (mean_df_ne$avg_phen_ne - mean(mean_df_ne$avg_phen_ne))/sd(mean_df_ne$avg_phen_ne) 
  
  return(list(model_df,mean_df,df.ne,mean_df_ne))
}

mod.GxE <- function(input_df){ # input is model_df
  
  # Clear outputs
  allGE <- c()
  loopGxE <- c()
  
  # Anova
  aov.test <- lm(phen_corrected ~ exp_env_factor * gen_factor, data = input_df)
  
  # Estimated Marginal Means
  emm_E = as.data.frame(emmeans(aov.test,"exp_env_factor"))
  emm_G = as.data.frame(emmeans(aov.test, "gen_factor"))
  emm_GxE = as.data.frame(emmeans(aov.test, ~ exp_env_factor*gen_factor))
  
  # Gmeans
  E_means <- tapply(emm_GxE$emmean, emm_GxE$exp_env_factor, mean)
  G_means <- tapply(emm_GxE$emmean, emm_GxE$gen_factor, mean)
  G_matrix <- data.frame("G_means" = G_means, "gen_factor" = unique(emm_GxE$gen_factor))
  E_matrix <- data.frame("E_means" = E_means, "exp_env_factor" = unique(emm_GxE$exp_env_factor))
  
  # Match Genotypes to Native Environment
  Cov_matrix = G_matrix
  Cov_matrix$exp_env_factor <- input_df$nat_env_factor[match(G_matrix$gen_factor,input_df$gen_factor)]
  Cov_matrix$E_means <- E_matrix$E_means[match(Cov_matrix$exp_env_factor,E_matrix$exp_env_factor)]
  
  # Magnitude of GxE -- EMMs
  GxE_emm_original<- abs(emm_GxE$emmean[emm_GxE$gen_factor == "G_1" & emm_GxE$exp_env_factor == "E_1"] - # GxE (Phenotype of ith genotype in jth environment)
                         emm_G$emmean[emm_G$gen_factor == "G_1"] - # phenotype of ith Genotype
                         emm_E$emmean[emm_E$exp_env_factor == "E_1"] + # phenotype of jth Environment
                         mean(emm_GxE$emmean)) # Overall mean phenotype
  
  # Magnitude of GxE -- Loop
  allGE = c()
  loopGxE = NULL
  for (i in 1:nlevels(emm_GxE$gen_factor)){
    for (j in 1:nlevels(emm_GxE$exp_env_factor)){
      G_levels <- levels(emm_GxE$gen_factor)
      E_levels <- levels(emm_GxE$exp_env_factor)
      loopGxE <- abs(emm_GxE$emmean[emm_GxE$gen_factor == G_levels[i] & emm_GxE$exp_env_factor == E_levels[j]] - # GxE (Phenotype of ith genotype in jth environment)
                    emm_G$emmean[emm_G$gen_factor == G_levels[i]] - # phenotype of ith Genotype
                    emm_E$emmean[emm_E$exp_env_factor == E_levels[j]] + # phenotype of jth Environment
                    mean(emm_GxE$emmean)) # Overall mean
      allGE <- c(allGE, loopGxE)
    }
  }
  #hist(allGE)
  GxE_emm_loop = mean(allGE)
  
  # Omega^2
  w2_GxE = (summary(aov(aov.test))[[1]][3,2] - #(SS_effect -
           (summary(aov(aov.test))[[1]][3,1]*summary(aov(aov.test))[[1]][4,3])) / #(Df_effect * MS_error))/
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
  
  return(list(Cov_matrix, GxE_emm_original, GxE_emm_loop, allGE, w2_GxE, eta_GxE, GxE_SumsSquares, mod_df))
}

mean.GxE <- function(input_df,is.perm){ # input is mean_df
  
  # Clear outputs
  allGEmeans <- c()
  GxE_mean.temp <- c()
  GiEj_mean = Gi_mean = Ej_mean = GiEj_null_samp = NULL
  
  if(is.perm == FALSE){

  # Means of Means
  E_means <- tapply(input_df$avg_phen_corrected, input_df$exp_env_factor, mean)
  G_means <- tapply(input_df$avg_phen_corrected, input_df$gen_factor, mean)
  Gmean_mat <- data.frame("G_means" = G_means, "gen_factor" = unique(input_df$gen_factor))
  Emean_mat <- data.frame("E_means" = E_means, "exp_env_factor" = unique(input_df$exp_env_factor))
  
  # Match means to native
  Cov_mean_matrix = Gmean_mat
  Cov_mean_matrix$exp_env_factor <- input_df$nat_env_factor[match(Cov_mean_matrix$gen_factor,input_df$gen_factor)]
  Cov_mean_matrix$E_means <- Emean_mat$E_means[match(Cov_mean_matrix$exp_env_factor,Emean_mat$exp_env_factor)]
  
  # Magnitude of GxE -- Loop -- Means
  for (i in 1:nlevels(input_df$gen_factor)){
    for (j in 1:nlevels(input_df$exp_env_factor)){
      G_levels <- levels(input_df$gen_factor)
      E_levels <- levels(input_df$exp_env_factor)
      GxE_mean.temp <- abs(input_df$avg_phen_corrected[input_df$gen_factor == G_levels[i] & input_df$exp_env_factor == E_levels[j]] - # GxE (Phenotype of ith genotype in jth environment)
                            mean(input_df$avg_phen_corrected[input_df$gen_factor == G_levels[i]])- # mean phenotype of ith Genotype
                            mean(input_df$avg_phen_corrected[input_df$exp_env_factor == E_levels[j]])+ # mean phenotype of jth Environment
                            mean(input_df$avg_phen_corrected)) # Overall mean
      allGEmeans <- c(allGEmeans, GxE_mean.temp)
    }
  }
  #hist(allGEmeans)
  GxE_means = mean(allGEmeans)
  
  }else{ # Below is modified code to generate null distribution for GxE
   
    # Means of Means
    E_means <- tapply(input_df$avg_phen_corrected, input_df$exp_env_factor, mean)
    G_means <- tapply(input_df$avg_phen_corrected, input_df$gen_factor, mean)
    Gmean_mat <- data.frame("G_means" = G_means, "gen_factor" = unique(input_df$gen_factor))
    Emean_mat <- data.frame("E_means" = E_means, "exp_env_factor" = unique(input_df$exp_env_factor))
    
    # Match means to native
    Cov_mean_matrix = Gmean_mat
    Cov_mean_matrix$exp_env_factor <- input_df$nat_env_factor[match(Cov_mean_matrix$gen_factor,input_df$gen_factor)]
    Cov_mean_matrix$E_means <- Emean_mat$E_means[match(Cov_mean_matrix$exp_env_factor,Emean_mat$exp_env_factor)]
    
    for (i in 1:nlevels(input_df$gen_factor)){
      for (j in 1:nlevels(input_df$exp_env_factor)){
        
        G_levels <- levels(input_df$gen_factor)
        E_levels <- levels(input_df$exp_env_factor)
        
        GiEj_mean <- input_df$avg_phen_corrected[input_df$gen_factor == G_levels[i] & input_df$exp_env_factor == E_levels[j]]
        Gi_mean <- mean(input_df$avg_phen_corrected[input_df$gen_factor == G_levels[i]])
        Ej_mean <- mean(input_df$avg_phen_corrected[input_df$exp_env_factor == E_levels[j]])
        
        # Create a sample of the null expectation for the GiEj
        GiEj_null_samp <- rnorm(1, mean = (Gi_mean + Ej_mean), sd = mean(input_df$se[input_df$gen_factor == G_levels[i] & input_df$exp_env_factor == E_levels[j]]))
        
        # Estimate 
        GxE_mean.temp <- abs(GiEj_null_samp - # GxE (Phenotype of ith genotype in jth environment)
                             Gi_mean - # mean phenotype of ith Genotype
                             Ej_mean + # mean phenotype of jth Environment
                             mean(input_df$avg_phen_corrected)) # Overall mean
        allGEmeans <- c(allGEmeans, GxE_mean.temp)
      }
    }
    
    GxE_means = mean(allGEmeans)
    
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
  new_phen.<- new_phen <-NULL
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
      new_phen. <- rnorm(nrow(cond), mean = cond$avg_phen, sd = cond$se) # generate replicate mean
      new_phen <- sample(new_phen., size = length(new_phen.), replace = TRUE) # shuffle

      # Output
      new_mean_temp <- data.frame("gen_factor" = cond$gen_factor,
                                  "exp_env_factor" = cond$exp_env_factor,
                                  "nat_env_factor" = cond$nat_env_factor,
                                  "mean_phen" = new_phen.)
      new_means <- rbind(new_means, new_mean_temp)
    }
  }
  
  # Standardize resampled means
  new_means$avg_phen_corrected = (new_means$mean_phen - mean(new_means$mean_phen))/sd(new_means$mean_phen) 
  
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

permutation_means <- function(input_df,seed){ # means dataframe (mean_df)
  
  # Clear outputs
  perm_means <- data.frame()
  null_gen = null_env = null_means = NULL
  
  # Shuffle means data (same set.seed keeps phen and corresponding se matched)
  set.seed(seed)
  null_means. <- sample(input_df$avg_phen, size = length(input_df$avg_phen), replace = FALSE)
  
  set.seed(seed)
  null_se <- sample(input_df$se, size = length(input_df$se), replace = FALSE)
  
  null_means <- rnorm(length(null_means.), mean = null_means., sd = null_se) # create replicate mean
  
  perm_means <- data.frame("gen_factor" = input_df$gen_factor,
                           "exp_env_factor" = input_df$exp_env_factor,
                           "nat_env_factor" = input_df$nat_env_factor,
                           "avg_phen" = null_means,
                           "se" = null_se)
  # Restandardize
  perm_means$avg_phen_corrected = (perm_means$avg_phen - mean(perm_means$avg_phen))/sd(perm_means$avg_phen)
  
  return(perm_means)
}

pvalue_fun <- function(estimate, rankdat, test){ #Test = "twotail" or "righttail"
  
  p.value = NULL
  
  if(test == "twotail"){
    p.value = sum(abs(rankdat) >= abs(estimate))/(nperms+1) # Two-tailed

  }else if(test == "righttail"){
    p.value = sum(rankdat >= estimate)/(nperms+1) # Right-tailed
  
    }else{p.value = "Invalid test entry- do you mean twotail or righttail?"}
   
  return(p.value)
}


