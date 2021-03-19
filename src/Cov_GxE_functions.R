###################################################################################
##              Functions for co/counter gradient simulations                    ##
##        Authors: Molly Albecker, Geoff Trussell, Katie Lotterhos               ##
###################################################################################

df.foundations <- function(delta_env, delta_gen, sample_size, n_env, std_dev, n_pop, interaction, seed1, seed2){
  
  # Dataframe foundations
  n_environments <- n_env 
  gen <- rep(1:n_pop, each = sample_size*n_environments)
  env <- rep(1:n_environments, times = n_pop, each = sample_size) 
  
  # Interaction Term
  set.seed = seed2
  int <- rep(rnorm(n_pop * n_environments, 0, sd = interaction), each = sample_size) # interaction term - one for each GE level
  
  ############################################
  ## Phenotype data with standard deviation ##
  ############################################
  
  # Create model dataframe 
  model_df <- data.frame(gen, env, int)
  model_df$gen_factor <- factor(paste("G", model_df$gen, sep = "_"))
  model_df$exp_env_factor <- factor(paste("E", model_df$env, sep = "_"))
  model_df$GE_factor <- paste0("G",model_df$gen, "E", model_df$env)
  
  # Native environments
  model_df$nat_env_factor = rep(unique(model_df$exp_env_factor), each = (sample_size*n_pop))
  
  # True mean phenotype data using regression equation
  model_df$GE_true = delta_env * model_df$env + delta_gen * model_df$gen + model_df$int 
  
  # True means 
  G_true <- data.frame(G_true = tapply(model_df$GE_true, model_df$gen_factor, mean))
  G_true <- rownames_to_column(G_true, var = "gen_factor")
  E_true <- data.frame(E_true = tapply(model_df$GE_true, model_df$exp_env_factor, mean))
  E_true <- rownames_to_column(E_true, var = "exp_env_factor")
  
  model_df <- full_join(model_df, G_true, by = "gen_factor")
  model_df <- full_join(model_df, E_true, by = "exp_env_factor")
  model_df$mean_true <- mean(model_df$GE_true)
  
  # Calculate true interaction term
  model_df$int_true <- model_df$GE_true - model_df$G_true - model_df$E_true + model_df$mean_true # GEij - Gmean - Emean + Overall
  
  # Standardize true means 
  GE_true_means <- tapply(model_df$GE_true, model_df$GE_factor, mean)
  model_df$GE_stn_true <- (model_df$GE_true - mean(GE_true_means))/sd(GE_true_means)
  
  # Add random noise
  set.seed = seed1
  model_df$e <- rnorm(sample_size * n_pop * n_environments, 0, sd = std_dev) 
  
  # Add error to standardized phenotype
  model_df$phen_corrected <- model_df$GE_stn_true + model_df$e
  
  # Chr -> Factor
  model_df$gen_factor <- as.factor(model_df$gen_factor)
  model_df$exp_env_factor <- as.factor(model_df$exp_env_factor)
  
  # Sanity plots
  # ggplot(model_df, aes(x = exp_env_factor, y = GE_stn_true, group = gen_factor, colour = nat_env_factor))+ geom_point()+geom_line()
  # ggplot(model_df, aes(x = exp_env_factor, y = GE_true, group = gen_factor, colour = gen_factor))+ geom_point()+geom_line()
  # ggplot(model_df, aes(x = exp_env_factor, y = phen_corrected, group = gen_factor, colour = gen_factor))+ geom_point()+geom_smooth(se=FALSE)
  
  # Generate means dataframe
  GE_means <- data.frame(avg_phen_corrected = tapply(model_df$phen, model_df$GE_factor, mean))
  GE_means <- rownames_to_column(GE_means, "GE_factor")
  GE = data.frame(GE_factor = unique(model_df$GE_factor),
                  gen_factor = rep(unique(model_df$gen_factor), each = n_pop),
                  exp_env_factor = rep(unique(model_df$exp_env_factor), times = n_pop))
  mean_df <- full_join(GE, GE_means, by="GE_factor")
  sd_avg = model_df %>%
    group_by(gen_factor, exp_env_factor) %>%
    summarize("deviation" = mean(abs(e)))
  mean_df <- full_join(mean_df, sd_avg)
  mean_df$se <- mean_df$deviation/(sqrt(sample_size))
  mean_df$nat_env_factor <- model_df$nat_env_factor[match(mean_df$gen_factor,model_df$gen_factor)]
  
  # Means Plots
  #ggplot(mean_df, aes(x = exp_env_factor, y = avg_phen_corrected, group = gen_factor, colour = gen_factor))+geom_line()
  
  # Chr -> Factor
  mean_df$gen_factor <- as.factor(mean_df$gen_factor)
  mean_df$exp_env_factor <- as.factor(mean_df$exp_env_factor)
  
  ###############################################
  ## Phenotype data with no standard deviation ##
  ###############################################
  
  # Raw data
  df.ne = data.frame(gen_factor = model_df$gen_factor, 
                     exp_env_factor = model_df$exp_env_factor, 
                     nat_env_factor = model_df$nat_env_factor,
                     phen_corrected = model_df$GE_stn_true)
  
  # Generate means dataframe
  GE_means_ne <- data.frame(avg_phen_corrected = tapply(model_df$GE_stn_true, model_df$GE_factor, mean))
  GE_means_ne <- rownames_to_column(GE_means_ne, "GE_factor")
  mean_df_ne <- full_join(GE, GE_means_ne, by="GE_factor")
  mean_df_ne$se <- 0
  
  # Add native environments
  mean_df_ne$nat_env_factor = df.ne$nat_env_factor[match(mean_df_ne$gen_factor,df.ne$gen_factor)]
  
  # Chr -> Factor
  mean_df_ne$gen_factor <- as.factor(mean_df_ne$gen_factor)
  mean_df_ne$exp_env_factor <- as.factor(mean_df_ne$exp_env_factor)
  
  # Create Variance Partitioning dataframe
  varpar.df = data.frame(model_df)
  
  G_stn_est <- data.frame("G_stn_est" = tapply(model_df$phen_corrected, model_df$gen_factor, mean)) # Gmeans equivalent - yi_bar
  G_stn_est <- rownames_to_column(G_stn_est, "gen_factor")
  varpar.df <- full_join(varpar.df, G_stn_est)
  
  E_stn_est <- data.frame("E_stn_est" = tapply(model_df$phen_corrected, model_df$exp_env_factor, mean)) # Emeans equivalent - yj_bar
  E_stn_est <- rownames_to_column(E_stn_est, "exp_env_factor")
  varpar.df <- full_join(varpar.df, E_stn_est)
  
  GE_stn_est <- data.frame("GE_stn_est" = tapply(model_df$phen_corrected, model_df$GE_factor, mean))
  GE_stn_est <- rownames_to_column(GE_stn_est, "GE_factor")
  varpar.df <- full_join(varpar.df, GE_stn_est)
  
  varpar.df$mean_stn_est = mean(model_df$phen_corrected)
  varpar.df$int_stn_est = varpar.df$mean_stn_est + varpar.df$GE_stn_est - varpar.df$E_stn_est - varpar.df$G_stn_est
  
  varpar.df$ind <- varpar.df$env == varpar.df$gen
  for(i in 1:nrow(varpar.df)){
    if(varpar.df$ind[i]){varpar.df$I[i] = 1}else{varpar.df$I[i] = 0}
  }
  
  return(list(model_df,mean_df,df.ne,mean_df_ne,varpar.df))
}

df.foundations2 <- function(delta_env, delta_gen, sample_size, n_env, std_dev, n_pop, interaction, seed1, seed2, seed3){
  
  # Dataframe foundations
  n_environments <- n_env 
  gen <- rep(1:n_environments, each = sample_size*n_pop)
  env <- rep(1:n_environments, times = n_pop, each = sample_size) 
  
  # Interaction Term
  set.seed = seed2
  int <- rep(rnorm(n_pop * n_environments, 0, sd = interaction), each = sample_size) # interaction term - one for each GE level
  
  ############################################
  ## Phenotype data with standard deviation ##
  ############################################
  
  # Create model dataframe 
  model_df <- data.frame(gen, env, int)
  
  # True mean phenotype data using regression equation
  model_df$GE_true = delta_env * model_df$env + delta_gen * model_df$gen + model_df$int 
  
  # Now add in actual genotypes
  model_df$gen <- rep(1:n_pop, each = sample_size*n_env)
  model_df$env <- rep(1:n_environments, times = n_pop, each = sample_size) 
  model_df$gen_factor <- factor(paste("G", model_df$gen, sep = "_"))
  model_df$exp_env_factor <- factor(paste("E", model_df$env, sep = "_"))
  model_df$GE_factor <- paste0("G",model_df$gen, "E", model_df$env)
  
  # True means 
  G_true <- data.frame(G_true = tapply(model_df$GE_true, model_df$gen_factor, mean))
  G_true <- rownames_to_column(G_true, var = "gen_factor")
  E_true <- data.frame(E_true = tapply(model_df$GE_true, model_df$exp_env_factor, mean))
  E_true <- rownames_to_column(E_true, var = "exp_env_factor")
  
  model_df <- full_join(model_df, G_true, by = "gen_factor")
  model_df <- full_join(model_df, E_true, by = "exp_env_factor")
  model_df$mean_true <- mean(model_df$GE_true)
  
  # Calculate true interaction term
  model_df$int_true <- model_df$GE_true - model_df$G_true - model_df$E_true + model_df$mean_true # GEij - Gmean - Emean + Overall
  
  # Native environments
  model_df$nat_env_factor = rep(unique(model_df$exp_env_factor), each = (sample_size*n_pop))
  
  # Standardize true means 
  GE_true_means <- tapply(model_df$GE_true, model_df$GE_factor, mean)
  model_df$GE_stn_true <- (model_df$GE_true - mean(GE_true_means))/sd(GE_true_means)
  
  # Add random noise
  set.seed = seed1
  model_df$e <- rnorm(sample_size * n_pop * n_environments, 0, sd = std_dev) 
  
  # Phenotype with error
  model_df$phen_corrected <- model_df$GE_stn_true + model_df$e
  
  # Chr -> Factor
  model_df$gen_factor <- as.factor(model_df$gen_factor)
  model_df$exp_env_factor <- as.factor(model_df$exp_env_factor)
  
  # Sanity plots
  #ggplot(model_df, aes(x = exp_env_factor, y = GE_stn_true, group = gen_factor, colour = nat_env_factor))+ geom_point()+geom_line(size = 1) + theme_classic()
  #ggplot(model_df, aes(x = exp_env_factor, y = GE_true, group = gen_factor, colour = gen_factor))+ geom_point()+geom_line()
  #ggplot(model_df, aes(x = exp_env_factor, y = phen_corrected, group = gen_factor, colour = nat_env_factor))+ geom_point()+geom_smooth(se=FALSE,size = 1) + theme_classic()
  
  # Generate means dataframe
  GE_means <- data.frame(avg_phen_corrected = tapply(model_df$phen, model_df$GE_factor, mean))
  GE_means <- rownames_to_column(GE_means, "GE_factor")
  GE = data.frame(GE_factor = unique(model_df$GE_factor),
                  gen_factor = rep(unique(model_df$gen_factor), each = n_env),
                  exp_env_factor = rep(unique(model_df$exp_env_factor), times = n_env))
  mean_df <- full_join(GE, GE_means, by="GE_factor")
  sd_avg = model_df %>%
    group_by(gen_factor, exp_env_factor) %>%
    summarize("deviation" = mean(abs(e)))
  mean_df <- full_join(mean_df, sd_avg)
  mean_df$se <- mean_df$deviation/(sqrt(sample_size))
  mean_df$nat_env_factor <- model_df$nat_env_factor[match(mean_df$gen_factor,model_df$gen_factor)]
  
  # Chr -> Factor
  mean_df$gen_factor <- as.factor(mean_df$gen_factor)
  mean_df$exp_env_factor <- as.factor(mean_df$exp_env_factor)
  
  # Means Plot
  # ggplot(mean_df, aes(x = exp_env_factor, y = avg_phen_corrected, group = gen_factor, colour = nat_env_factor))+geom_line()
  
  ###############################################
  ## Phenotype data with no standard deviation ##
  ###############################################
  
  # Raw data
  df.ne = data.frame(gen_factor = model_df$gen_factor, 
                     exp_env_factor = model_df$exp_env_factor, 
                     nat_env_factor = model_df$nat_env_factor,
                     phen_corrected = model_df$GE_stn_true)
  
  # Generate means dataframe
  GE_means_ne <- data.frame(avg_phen_corrected = tapply(model_df$GE_stn_true, model_df$GE_factor, mean))
  GE_means_ne <- rownames_to_column(GE_means_ne, "GE_factor")
  mean_df_ne <- full_join(GE, GE_means_ne, by="GE_factor")
  mean_df_ne$se <- 0
  
  # Add native environments
  mean_df_ne$nat_env_factor = df.ne$nat_env_factor[match(mean_df_ne$gen_factor,df.ne$gen_factor)]
  
  # Chr -> Factor
  mean_df_ne$gen_factor <- as.factor(mean_df_ne$gen_factor)
  mean_df_ne$exp_env_factor <- as.factor(mean_df_ne$exp_env_factor)
  
  # Create Variance Partitioning dataframe
  varpar.df = data.frame(model_df)
  
  G_stn_est <- data.frame("G_stn_est" = tapply(model_df$phen_corrected, model_df$gen_factor, mean)) # Gmeans equivalent - yibar
  G_stn_est <- rownames_to_column(G_stn_est, "gen_factor")
  varpar.df <- full_join(varpar.df, G_stn_est)
  
  E_stn_est <- data.frame("E_stn_est" = tapply(model_df$phen_corrected, model_df$exp_env_factor, mean)) # Gmeans equivalent - yibar
  E_stn_est <- rownames_to_column(E_stn_est, "exp_env_factor")
  varpar.df <- full_join(varpar.df, E_stn_est)
  
  GE_stn_est <- data.frame("GE_stn_est" = tapply(model_df$phen_corrected, model_df$GE_factor, mean))
  GE_stn_est <- rownames_to_column(GE_stn_est, "GE_factor")
  varpar.df <- full_join(varpar.df, GE_stn_est)
  
  varpar.df$mean_stn_est = mean(model_df$phen_corrected)
  varpar.df$int_stn_est = varpar.df$mean_stn_est + varpar.df$GE_stn_est - varpar.df$E_stn_est - varpar.df$G_stn_est
  
  varpar.df$ind <- varpar.df$env == varpar.df$gen
  for(i in 1:nrow(varpar.df)){
    if(varpar.df$ind[i]){varpar.df$I[i] = 1}else{varpar.df$I[i] = 0}
  }
  
  return(list(model_df,mean_df,df.ne,mean_df_ne,varpar.df))
}

var.partition <- function(input_df){ 
  
  # Variance due to genotype
  V_g = sum((input_df$G_stn_est-input_df$mean_stn_est)^2)
  
  # Variance due to environment
  V_e = sum((input_df$E_stn_est-input_df$mean_stn_est)^2)
  
  # Variance due to GxE
  V_gxe = sum(input_df$int_stn_est^2)
  
  # Variance due to CovGE
  V_cov <-  nrow(input_df)/sum(input_df$I)*
    (sum((input_df$G_stn_est-input_df$mean_stn_est)*(input_df$E_stn_est-input_df$mean_stn_est)*input_df$I))
  
  # Variance due to error
  V_error = sum((input_df$phen_corrected - input_df$GE_stn_est)^2)
  
  SS <- round(rbind(V_g, V_e, V_gxe, V_cov, V_error),4)
  omega2 <- round(abs(SS)/sum(abs(SS)),4)
  var_part = data.frame(SS, abs(SS), omega2)
  var_part = rownames_to_column(var_part, "Variance_Component")
  
  return(var_part)
}

mod.GxE <- function(input_df,is.perm,seed){ # input is model_df
  
  # Outputs
  allGE <- c()
  loopGxE <- c()
  
  # Establish Native Environments
  native_df = data.frame("gen_factor" = unique(input_df$gen_factor))
  native_df$nat_env_factor = input_df$nat_env_factor[match(native_df$gen_factor,input_df$gen_factor)]
  
  # Anova
  aov.test <- lm(phen_corrected ~ exp_env_factor * gen_factor, data = input_df) 
  
  # Estimated Marginal Means
  emm_E = as.data.frame(emmeans(aov.test,"exp_env_factor"))
  emm_G = as.data.frame(emmeans(aov.test, "gen_factor"))
  emm_GxE = as.data.frame(emmeans(aov.test, ~ exp_env_factor*gen_factor))
  
  # Gmeans
  G_matrix <- data.frame("G_means" = emm_G$emmean, "gen_factor" = emm_G$gen_factor)
  E_matrix <- data.frame("E_means" = emm_E$emmean, "exp_env_factor" = emm_E$exp_env_factor)
  
  # Covariance 
  Cov_matrix = G_matrix
  Cov_matrix$exp_env_factor = native_df$nat_env_factor[match(G_matrix$gen_factor,native_df$gen_factor)]
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
  
  if(is.perm == FALSE){ # No seed needed here
  
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
                       emm_G$emmean[emm_G$gen_factor == G_levels[i]] - # mean phenotype of ith Genotype
                       emm_E$emmean[emm_E$exp_env_factor == E_levels[j]] + # mean phenotype of jth Environment
                       mean(emm_GxE$emmean)) # Overall mean
      allGE <- c(allGE, loopGxE)
    }
  }

    GxE_emm_loop = mean(allGE)
  
  }else{ # Below generates null distribution for null of G+E means
    
    allGE <- NULL
    for (i in 1:nlevels(input_df$gen_factor)){
      for (j in 1:nlevels(input_df$exp_env_factor)){
        
        G_levels <- levels(input_df$gen_factor)
        E_levels <- levels(input_df$exp_env_factor)
        
        Gi_mean <- mean(sample(input_df$phen_corrected[input_df$gen_factor == G_levels[i]], 
                               size = length(input_df$phen_corrected[input_df$gen_factor == G_levels[i]]),
                               replace = TRUE))
        Ej_mean <- mean(sample(input_df$phen_corrected[input_df$exp_env_factor == E_levels[j]],
                               size = length(input_df$phen_corrected[input_df$exp_env_factor == E_levels[j]]),
                               replace = TRUE))
        GE_sd <- input_df %>%
          filter(gen_factor == G_levels[i]) %>%
          filter(exp_env_factor == E_levels[j]) %>%
          summarize("GEsd" = mean(e))

        # Create a sample of the null expectation for the Gi+Ej
        set.seed = seed
        GiEj_null_samp <- rnorm(1, mean = (Gi_mean + Ej_mean), sd = abs(GE_sd[[1]]))
        
        # Estimate 
        GxE_mean.temp <- abs(GiEj_null_samp - # G+E (Phenotype of ith genotype in jth environment)
                             Gi_mean - # mean phenotype of ith Genotype
                             Ej_mean + # mean phenotype of jth Environment
                             mean(input_df$phen_corrected)) # Overall mean
        allGE <- c(allGE, GxE_mean.temp)
      }
    }
    
    GxE_emm_loop = mean(allGE)
    
  }
  return(list(Cov_matrix, GxE_emm_original, GxE_emm_loop, allGE, w2_GxE, eta_GxE, GxE_SumsSquares, mod_df, delta_E, delta_H, aov.coefs,emm_GxE))
}

mean.GxE <- function(input_df,is.perm, seed){ # input is mean_df
  
  # Clear outputs
  allGEmeans <- c()
  GxE_mean.temp <- c()
  GiEj_mean = Gi_mean = Ej_mean = GiEj_null_samp = NULL
  
  if(is.perm == FALSE){ # No seed needed here
    
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
   
    # hist(allGEmeans)
    GxE_means = mean(allGEmeans)
    
  }else{ # Below generates null distribution for null of G+E means
    
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
        set.seed = seed
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
      new_phen <- rnorm(nrow(cond), mean = cond$avg_phen_corrected, sd = cond$se) # generate replicate means
      #set.seed(bootseed2)
      #new_phen <- sample(new_phen., size = length(new_phen.), replace = TRUE) # shuffle
      
      # Output
      new_mean_temp <- data.frame("gen_factor" = cond$gen_factor,
                                  "exp_env_factor" = cond$exp_env_factor,
                                  "nat_env_factor" = cond$nat_env_factor,
                                  "avg_phen_corrected" = new_phen)
      new_means <- rbind(new_means, new_mean_temp)
    }
  }
  
  # Standardize resampled means
  # new_means$avg_phen_corrected = (new_means$mean_phen - mean(new_means$mean_phen))/sd(new_means$mean_phen) 
  
  return(new_means)
}

permutation_raw <- function(input_df, perm.seed){ # input is model_df
  
  # Clear outputs
  perm_dat = data.frame()
  null_temp <- NULL
  
  # Shuffle raw data
  set.seed(perm.seed)
  null_temp <- sample(input_df$phen_corrected, size=nrow(input_df), replace=FALSE)
  
  perm_dat <- data.frame("gen_factor" = input_df$gen_factor,
                         "exp_env_factor" = input_df$exp_env_factor,
                         "nat_env_factor" = input_df$nat_env_factor,
                         "e" = input_df$e,
                         "phen_corrected" = null_temp)
  return(perm_dat)
}

permutation_means <- function(input_df, perm.seed2){ # means dataframe (mean_df)
  
  # Clear outputs
  perm_means <- data.frame()
  null_gen = null_env = null_means = NULL
  
  # Shuffle means data (same set.seed keeps phen and corresponding se matched)
  set.seed(perm.seed2)
  null_means <- sample(input_df$avg_phen_corrected, size = length(input_df$avg_phen_corrected), replace = FALSE)
  
  set.seed(perm.seed2)
  null_se <- sample(input_df$se, size = length(input_df$se), replace = FALSE)
  
  #null_means <- rnorm(length(null_means.), mean = null_means., sd = null_se) # create replicate mean
  
  perm_means <- data.frame("gen_factor" = input_df$gen_factor,
                           "exp_env_factor" = input_df$exp_env_factor,
                           "nat_env_factor" = input_df$nat_env_factor,
                           "avg_phen_corrected" = null_means,
                           "se" = null_se)
  # Restandardize
  # perm_means$avg_phen_corrected = (perm_means$avg_phen - mean(perm_means$avg_phen))/sd(perm_means$avg_phen)
  
  return(perm_means)
}

permutation_varpar <- function(input_df){
  
  # Clear outputs
  var_dat = data.frame()
  null_temp1 = null_temp2 = null_temp3 = null_temp4 = null_temp5 <- NULL
  
  # Shuffle raw data
  null_temp1 <- sample(input_df$phen_corrected, size=nrow(input_df), replace=FALSE)
  #null_temp2 <- sample(input_df$G_stn_est, size=nrow(input_df), replace=FALSE)
  #null_temp3 <- sample(input_df$E_stn_est, size=nrow(input_df), replace=FALSE)
  #null_temp4 <- sample(input_df$GE_stn_est, size=nrow(input_df), replace=FALSE)
  #null_temp5 <- sample(input_df$mean_stn_est, size=nrow(input_df), replace=FALSE)
  
  
  var_dat <- data.frame("gen_factor" = input_df$gen_factor,
                        "exp_env_factor" = input_df$exp_env_factor,
                        "GE_factor" = input_df$GE_factor,
                        "phen_corrected" = null_temp1,
                        "G_stn_est" = input_df$G_stn_est,
                        "int_stn_est" = input_df$int_stn_est,
                        "E_stn_est" = input_df$E_stn_est,
                        "GE_stn_est" = input_df$GE_stn_est,
                        "mean_stn_est" = input_df$mean_stn_est,
                        "I" = input_df$I)
  return(var_dat)
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

cov.function <- function(input_df, emm_df, is.sample = TRUE){ # input_df = cov matrix
  
  pvar <- function(x) { #population variance - for sample variance can use var() in R
    sum((x - mean(x))**2) / length(x)
  }
  
  N = length(input_df$gen_factor)
  overallmean = mean(emm_df$emmean,na.rm=TRUE) 
  numerator = sum((input_df$G_means - overallmean)*(input_df$E_means - overallmean))
  
  if(is.sample == TRUE){
    sample_correcter = max(var(input_df$E_means),var(input_df$G_means))
    cv = (1/(N-1))*(numerator/sample_correcter)
  }else{
    population_correcter = max(pvar(input_df$E_means),pvar(input_df$G_means))
    cv = (1/(N))*(numerator/population_correcter)
  }
  return(cv)
} 

cov.function_means <- function(input_df, phen_df, is.sample = TRUE){ # input_df = cov_matrix of G_means and E_means
  
  pvar <- function(x) { #population variance - for sample variance can use var() in R
    sum((x - mean(x))**2) / length(x)
  }
  
  N = length(input_df$gen_factor)
  
  overallmean = mean(phen_df$avg_phen_corrected) #not mean of means, mean of overall data
  numerator = sum((input_df$G_means - overallmean)*(input_df$E_means - overallmean))

  if(is.sample == TRUE){
    sample_correcter = max(var(input_df$E_means),var(input_df$G_means))
    cv = (1/(N-1))*(numerator/sample_correcter)
  }else{
    population_correcter = max(pvar(input_df$E_means),pvar(input_df$G_means))
    cv = (1/(N))*(numerator/population_correcter)
  }
  return(cv)
}

########### Simulation Summary/Plotting Functions ####################
is.empty <- function(x, mode=NULL){
  if (is.null(mode)) mode <- class(x)
  identical(vector(mode,1),c(x,vector(class(x),1)))
}

## Confusion Matrix data wrangling ##
fpr.fnr <- function(input_df, divided, scenario){
  
  is.empty <- function(x, mode=NULL){
    if (is.null(mode)) mode <- class(x)
    identical(vector(mode,1),c(x,vector(class(x),1)))
  }
  
  df <- data.frame()
  if(divided == TRUE){
    
  for(i in 1:length(unique(input_df$sample_size))){
    for(j in 1:length(unique(input_df$n_pop))){
      
      fn1 = fn = fp1 = fp = tn1 = tn = tp1 = tp = fnr = fpr = NULL
      ss = unique(input_df$sample_size)[i]
      np = unique(input_df$n_pop)[j]
      
      tempdf <- input_df %>% 
        filter(sample_size == ss) %>% 
        filter(n_pop == np)
      
      if(nrow(tempdf)==0){next}
      
      fn1 = tempdf$n[tempdf$name == "False Negative"]
      fp1 = tempdf$n[tempdf$name == "False Positive"]
      tn1 = tempdf$n[tempdf$name == "True Negative"]
      tp1 = tempdf$n[tempdf$name == "True Positive"]
      
      if(is.empty(fn1) == TRUE){fn = 0}else{fn = fn1}
      if(is.empty(fp1) == TRUE){fp = 0}else{fp = fp1}
      if(is.empty(tn1) == TRUE){tn = 0}else{tn = tn1}
      if(is.empty(tp1) == TRUE){tp = 0}else{tp = tp1}
      
      fnr = fn/(fn+tp)
      fpr = fp/(fp+tn)
      
      n = c(fn, fp, tn, tp)
      rate = c(fnr,fpr,NA,NA)
      if(scenario == 1){n_env = unique(tempdf$n_pop)}else{n_env = 2}
      total = sum(n)
      percent = (n/total)*100
      
      df. <- data.frame("name" = c("False Negative", "False Positive", "True Negative", "True Positive"),
                        "sample_size" = rep(ss,4),
                        "n_pop" = rep(np,4),
                        "totsamp" = rep(unique(tempdf$n_pop) * n_env * unique(tempdf$sample_size),4),
                        "n" = n, 
                        "total" = total, 
                        "percent" = round(percent,2),
                        "rate" = round(rate,2))
      df <- rbind(df,df.)
    }
  }
      
    }else{
      
      fn1 = fn = fp1 = fp = tn1 = tn = tp1 = tp = fnr = fpr = NA
      
      fn1 = input_df$n[input_df$name == "False Negative"]
      fp1 = input_df$n[input_df$name == "False Positive"]
      tn1 = input_df$n[input_df$name == "True Negative"]
      tp1 = input_df$n[input_df$name == "True Positive"]
    
      
      if(is.empty(fn1) == TRUE){fn = 0}else{fn = fn1}
      if(is.empty(fp1) == TRUE){fp = 0}else{fp = fp1}
      if(is.empty(tn1) == TRUE){tn = 0}else{tn = tn1}
      if(is.empty(tp1) == TRUE){tp = 0}else{tp = tp1}
      
      fnr = fn/(fn+tp)
      fpr = fp/(fp+tn)
      
      n = c(fn, fp, tn, tp)
      rate = c(fnr,fpr,NA,NA)
      total = sum(n)
      percent = (n/total)*100
      
      df <- data.frame("name" = c("False Negative", "False Positive", "True Negative", "True Positive"),
                       "n" = n, 
                       "total" = total, 
                       "percent" = round(percent,2),
                       "rate" = round(rate,2))
     }
  return(df)
} # Result goes into heatmap_fun

## Confusion Matrix Heatmaps ## 
heatmap_fun <- function(plot_data, plot_type){ #plot_type is "percent" or "rate"
  
  p <- list()
  
  if(plot_type == "percent"){
  for(i in 1:length(unique(plot_data$name))){
    
    subcat <- filter(plot_data, name == unique(plot_data$name)[i])
    
    p[[i]] <- ggplot(subcat,aes(x = factor(sample_size), y = factor(n_pop), fill = percent)) + 
      geom_tile() + 
      geom_text(aes(label= paste(round(percent,3), totsamp,sep = '\n')), size = 5) +
      theme_classic(base_size = 12, base_family = "Times") + 
      scale_fill_gradient2(#low="#DDDDDD", mid="#99CCEE", high="#000044", #colors in the scale
        #                     midpoint=mean(rng.gxe2),    #same midpoint for plots (mean of the range)
        #                     breaks=seq(0,1,0.25), #breaks in the scale bar
        limits=c(0,1))+
      xlab("Sample Size") + ylab("Number of Populations")+
      #labs(fill = "Percent")+
      theme(axis.text = element_text(colour = "black"))+
      theme(legend.position = "none")+
      ggtitle(unique(subcat$name))+
      theme(plot.title = element_text(size = 16, face = "bold"))
  }
  }else{

    for(i in 1:length(unique(plot_data$name))){
      
      subcat <- filter(plot_data, name == unique(plot_data$name)[i])
      if(unique(is.na(subcat$rate))==TRUE){next}
      
      p[[i]] <- ggplot(subcat,aes(x = factor(sample_size), y = factor(n_pop), fill = rate)) + 
        geom_tile() + 
        geom_text(aes(label= paste(round(rate,3), totsamp,sep = '\n')), size = 5) +
        theme_classic(base_size = 12, base_family = "Times")+ 
        scale_fill_gradient2(limits=c(0,1))+
        xlab("Sample Size") + ylab("Number of Populations")+
        #labs(fill = "Percent")+
        theme(axis.text = element_text(colour = "black"))+
        theme(legend.position = "none")+
        ggtitle(unique(subcat$name))+
        theme(plot.title = element_text(size = 16, face = "bold"))
      
      
    
  }
  }
  return(do.call(grid.arrange, p))
} 

## Calculate Power and False Negative rates
fnr.effsize <- function(x, metric, data.type, analysis, scenario = 1, resolution){ # metric = Cov or GxE; analysis is perm or boot or anova
   # metric = "Cov" or "GxE"
   # data.type = "raw" or "means"
   # analysis = "perm" or "boot"
   # scenario = 1 or 2
   # resolution = "fine" (for heatmaps) or "coarse" (for barplots)
  
  output = data.frame()
  
  if(data.type == "raw"){

    if(metric == "Cov"){
    
      if(resolution == "coarse"){
        x$binCov = "NA"
        for(i in 1:nrow(x)){
        if(abs(x$true_cov[i]) > 0 & abs(x$true_cov[i]) <= 0.25){x$binCov[i] = 0.25
        }else if(abs(x$true_cov[i]) > 0.25 & abs(x$true_cov[i]) <= 0.5){x$binCov[i] = 0.5
        }else if(abs(x$true_cov[i]) > 0.5 & abs(x$true_cov[i]) <= 0.75){x$binCov[i] = 0.75
        }else{x$binCov[i] = 1}
      }
    }else{
      x$binCov = "NA"
      for(i in 1:nrow(x)){
        if(x$true_cov[i] == 0){x$binCov[i] = 0
        }else if(abs(x$true_cov[i]) > 0 & abs(x$true_cov[i]) <= 0.15){x$binCov[i] = 0.1
        }else if(abs(x$true_cov[i]) > 0.15 & abs(x$true_cov[i]) <= 0.25){x$binCov[i] = 0.2
        }else if(abs(x$true_cov[i]) > 0.25 & abs(x$true_cov[i]) <= 0.35){x$binCov[i] = 0.3
        }else if(abs(x$true_cov[i]) > 0.35 & abs(x$true_cov[i]) <= 0.45){x$binCov[i] = 0.4
        }else if(abs(x$true_cov[i]) > 0.45 & abs(x$true_cov[i]) <= 0.55){x$binCov[i] = 0.5
        }else if(abs(x$true_cov[i]) > 0.55 & abs(x$true_cov[i]) <= 0.65){x$binCov[i] = 0.6
        }else if(abs(x$true_cov[i]) > 0.65 & abs(x$true_cov[i]) <= 0.75){x$binCov[i] = 0.7
        }else if(abs(x$true_cov[i]) > 0.75 & abs(x$true_cov[i]) <= 0.85){x$binCov[i] = 0.8
        }else if(abs(x$true_cov[i]) > 0.85 & abs(x$true_cov[i]) <= 0.95){x$binCov[i] = 0.9
        }else{x$binCov[i] = 1}
      }
      
    }
      
      for(i in 1:length(unique(x$sample_size))){
        for(j in 1:length(unique(x$n_pop))){
          for(k in 1:length(unique(x$binCov))){
            
            fn1 = fn = fp1 = fp = tn1 = tn = tp1 = tp = fnr = fpr = NULL
            ss = unique(x$sample_size)[i]
            np = unique(x$n_pop)[j]
            bc = unique(x$binCov)[k]
            
            if(analysis == "perm"){
            
              tempdf <- x %>% 
                filter(sample_size == ss) %>% 
                filter(n_pop == np) %>%
                filter(binCov == bc) %>%
                group_by("name" = Covconfintperm,sample_size,n_pop,binCov)%>%
                summarise(n = n())
              
              if(nrow(tempdf)==0){next}
              
              fn1 = tempdf$n[tempdf$name == "False Negative"]
              fp1 = tempdf$n[tempdf$name == "False Positive"]
              tn1 = tempdf$n[tempdf$name == "True Negative"]
              tp1 = tempdf$n[tempdf$name == "True Positive"]
              
              if(is.empty(fn1) == TRUE){fn = 0}else{fn = fn1}
              if(is.empty(fp1) == TRUE){fp = 0}else{fp = fp1}
              if(is.empty(tn1) == TRUE){tn = 0}else{tn = tn1}
              if(is.empty(tp1) == TRUE){tp = 0}else{tp = tp1}
              
              fnr = fn/(fn+tp)
              
              output1 = data.frame(sample_size = ss, 
                                   n_pop = np, 
                                   bin = bc,
                                   fnr = fnr,
                                   power = 1-fnr)
              output = rbind(output, output1)
            
              }else{
              
                tempdf <- x %>% 
                filter(sample_size == ss) %>% 
                filter(n_pop == np) %>%
                filter(binCov == bc) %>%
                group_by("name" = Covconfintboot,sample_size,n_pop,binCov)%>%
                summarise(n = n())
                
                if(nrow(tempdf)==0){next}
                
                fn1 = tempdf$n[tempdf$name == "False Negative"]
                fp1 = tempdf$n[tempdf$name == "False Positive"]
                tn1 = tempdf$n[tempdf$name == "True Negative"]
                tp1 = tempdf$n[tempdf$name == "True Positive"]
                
                if(is.empty(fn1) == TRUE){fn = 0}else{fn = fn1}
                if(is.empty(fp1) == TRUE){fp = 0}else{fp = fp1}
                if(is.empty(tn1) == TRUE){tn = 0}else{tn = tn1}
                if(is.empty(tp1) == TRUE){tp = 0}else{tp = tp1}
                
                fnr = fn/(fn+tp)
                
                output1 = data.frame(sample_size = ss, 
                                     n_pop = np, 
                                     bin = bc,
                                     fnr = fnr,
                                     power = 1-fnr)
                output = rbind(output, output1)
              }
          }
        }
      }
          }else{
            
            if(resolution == "course"){
            x$binGxE = "NA"
            for(i in 1:nrow(x)){
              if(abs(x$true_GxE_emm[i]) > 0 & abs(x$true_GxE_emm[i]) <= 0.25){x$binGxE[i] = 0.25
              }else if(abs(x$true_GxE_emm[i]) > 0.25 & abs(x$true_GxE_emm[i]) <= 0.5){x$binGxE[i] = 0.5
              }else if(abs(x$true_GxE_emm[i]) > 0.5 & abs(x$true_GxE_emm[i]) <= 0.75){x$binGxE[i] = 0.75
              }else{x$binGxE[i] = 1}
            }
            
            }else{
              x$binGxE = "NA"
              for(i in 1:nrow(x)){
                if(x$true_GxE_emm[i] == 0){x$binGxE[i] = 0
                }else if(abs(x$true_GxE_emm[i]) > 0 & abs(x$true_GxE_emm[i]) <= 0.15){x$binGxE[i] = 0.1
                }else if(abs(x$true_GxE_emm[i]) > 0.15 & abs(x$true_GxE_emm[i]) <= 0.25){x$binGxE[i] = 0.2
                }else if(abs(x$true_GxE_emm[i]) > 0.25 & abs(x$true_GxE_emm[i]) <= 0.35){x$binGxE[i] = 0.3
                }else if(abs(x$true_GxE_emm[i]) > 0.35 & abs(x$true_GxE_emm[i]) <= 0.45){x$binGxE[i] = 0.4
                }else if(abs(x$true_GxE_emm[i]) > 0.45 & abs(x$true_GxE_emm[i]) <= 0.55){x$binGxE[i] = 0.5
                }else if(abs(x$true_GxE_emm[i]) > 0.55 & abs(x$true_GxE_emm[i]) <= 0.65){x$binGxE[i] = 0.6
                }else if(abs(x$true_GxE_emm[i]) > 0.65 & abs(x$true_GxE_emm[i]) <= 0.75){x$binGxE[i] = 0.7
                }else if(abs(x$true_GxE_emm[i]) > 0.75 & abs(x$true_GxE_emm[i]) <= 0.85){x$binGxE[i] = 0.8
                }else if(abs(x$true_GxE_emm[i]) > 0.85 & abs(x$true_GxE_emm[i]) <= 0.95){x$binGxE[i] = 0.9
                }else{x$binGxE[i] = 1}
              }
              
            }
            
            for(i in 1:length(unique(x$sample_size))){
              for(j in 1:length(unique(x$n_pop))){
                for(k in 1:length(unique(x$binGxE))){
                  
                  fn1 = fn = fp1 = fp = tn1 = tn = tp1 = tp = fnr = fpr = NULL
                  ss = unique(x$sample_size)[i]
                  np = unique(x$n_pop)[j]
                  bc = unique(x$binGxE)[k]
                  
                  if(analysis == "perm"){
                    
                    tempdf <- x %>% 
                      filter(sample_size == ss) %>% 
                      filter(n_pop == np) %>%
                      filter(binGxE == bc) %>%
                      group_by("name" = GxEconfintperm,sample_size,n_pop,binGxE)%>%
                      summarise(n = n())
                    
                    if(nrow(tempdf)==0){next}
                    
                    fn1 = tempdf$n[tempdf$name == "False Negative"]
                    fp1 = tempdf$n[tempdf$name == "False Positive"]
                    tn1 = tempdf$n[tempdf$name == "True Negative"]
                    tp1 = tempdf$n[tempdf$name == "True Positive"]
                    
                    if(is.empty(fn1) == TRUE){fn = 0}else{fn = fn1}
                    if(is.empty(fp1) == TRUE){fp = 0}else{fp = fp1}
                    if(is.empty(tn1) == TRUE){tn = 0}else{tn = tn1}
                    if(is.empty(tp1) == TRUE){tp = 0}else{tp = tp1}
                    
                    fnr = fn/(fn+tp)
                    
                    output1 = data.frame(sample_size = ss, 
                                         n_pop = np, 
                                         bin = bc,
                                         fnr = fnr,
                                         power = 1-fnr)
                    output = rbind(output, output1)
                    
                  }else if(analysis == "boot"){
                    
                    tempdf <- x %>% 
                      filter(sample_size == ss) %>% 
                      filter(n_pop == np) %>%
                      filter(binGxE == bc) %>%
                      group_by("name" = GxEconfintboot,sample_size,n_pop,binGxE)%>%
                      summarise(n = n())
                    
                    if(nrow(tempdf)==0){next}
                    
                    fn1 = tempdf$n[tempdf$name == "False Negative"]
                    fp1 = tempdf$n[tempdf$name == "False Positive"]
                    tn1 = tempdf$n[tempdf$name == "True Negative"]
                    tp1 = tempdf$n[tempdf$name == "True Positive"]
                    
                    if(is.empty(fn1) == TRUE){fn = 0}else{fn = fn1}
                    if(is.empty(fp1) == TRUE){fp = 0}else{fp = fp1}
                    if(is.empty(tn1) == TRUE){tn = 0}else{tn = tn1}
                    if(is.empty(tp1) == TRUE){tp = 0}else{tp = tp1}
                    
                    fnr = fn/(fn+tp)
                    
                    output1 = data.frame(sample_size = ss, 
                                         n_pop = np, 
                                         bin = bc,
                                         fnr = fnr,
                                         power = 1-fnr)
                    output = rbind(output, output1)
                    
                  }else{
                    
                    if(scenario == 1){
                    
                      tempdf <- x %>% 
                      filter(sample_size == ss) %>% 
                      filter(n_pop == np) %>%
                      filter(binGxE == bc) %>%
                      group_by("name" = GxEanova_conf,sample_size,n_pop,binGxE)%>%
                      summarise(n = n())
                    }else{
                    
                    tempdf <- x %>% 
                      filter(sample_size == ss) %>% 
                      filter(n_pop == np) %>%
                      filter(binGxE == bc) %>%
                      group_by("name" = GxEanova_conf,sample_size,n_pop,binGxE)%>%
                      summarise(n = n())
                    }
                  
                    
                    if(nrow(tempdf)==0){next}
                    
                    fn1 = tempdf$n[tempdf$name == "False Negative"]
                    fp1 = tempdf$n[tempdf$name == "False Positive"]
                    tn1 = tempdf$n[tempdf$name == "True Negative"]
                    tp1 = tempdf$n[tempdf$name == "True Positive"]
                    
                    if(is.empty(fn1) == TRUE){fn = 0}else{fn = fn1}
                    if(is.empty(fp1) == TRUE){fp = 0}else{fp = fp1}
                    if(is.empty(tn1) == TRUE){tn = 0}else{tn = tn1}
                    if(is.empty(tp1) == TRUE){tp = 0}else{tp = tp1}
                    
                    fnr = fn/(fn+tp)
                    
                    output1 = data.frame(sample_size = ss, 
                                         n_pop = np, 
                                         bin = bc,
                                         fnr = fnr,
                                         power = 1-fnr)
                    output = rbind(output, output1)
                  }
                }
              }
            }
          }
  }else{
    
    if(metric == "Cov"){
      
      if(resolution == "coarse"){
        x$binCov = "NA"
        for(i in 1:nrow(x)){
          if(abs(x$true_cov_means[i]) > 0 & abs(x$true_cov_means[i]) <= 0.25){x$binCov[i] = 0.25
          }else if(abs(x$true_cov_means[i]) > 0.25 & abs(x$true_cov_means[i]) <= 0.5){x$binCov[i] = 0.5
          }else if(abs(x$true_cov_means[i]) > 0.5 & abs(x$true_cov_means[i]) <= 0.75){x$binCov[i] = 0.75
          }else{x$binCov[i] = 1}
        }
      }else{
        x$binCov = "NA"
        for(i in 1:nrow(x)){
          if(x$true_cov_means[i] == 0){x$binCov[i] = 0
          }else if(abs(x$true_cov_means[i]) > 0 & abs(x$true_cov_means[i]) <= 0.15){x$binCov[i] = 0.1
          }else if(abs(x$true_cov_means[i]) > 0.15 & abs(x$true_cov_means[i]) <= 0.25){x$binCov[i] = 0.2
          }else if(abs(x$true_cov_means[i]) > 0.25 & abs(x$true_cov_means[i]) <= 0.35){x$binCov[i] = 0.3
          }else if(abs(x$true_cov_means[i]) > 0.35 & abs(x$true_cov_means[i]) <= 0.45){x$binCov[i] = 0.4
          }else if(abs(x$true_cov_means[i]) > 0.45 & abs(x$true_cov_means[i]) <= 0.55){x$binCov[i] = 0.5
          }else if(abs(x$true_cov_means[i]) > 0.55 & abs(x$true_cov_means[i]) <= 0.65){x$binCov[i] = 0.6
          }else if(abs(x$true_cov_means[i]) > 0.65 & abs(x$true_cov_means[i]) <= 0.75){x$binCov[i] = 0.7
          }else if(abs(x$true_cov_means[i]) > 0.75 & abs(x$true_cov_means[i]) <= 0.85){x$binCov[i] = 0.8
          }else if(abs(x$true_cov_means[i]) > 0.85 & abs(x$true_cov_means[i]) <= 0.95){x$binCov[i] = 0.9
          }else{x$binCov[i] = 1}
        }
        
      }
      
      for(i in 1:length(unique(x$sample_size))){
        for(j in 1:length(unique(x$n_pop))){
          for(k in 1:length(unique(x$binCov))){
            
            fn1 = fn = fp1 = fp = tn1 = tn = tp1 = tp = fnr = fpr = NULL
            ss = unique(x$sample_size)[i]
            np = unique(x$n_pop)[j]
            bc = unique(x$binCov)[k]
            
            if(analysis == "perm"){
              
              tempdf <- x %>% 
                filter(sample_size == ss) %>% 
                filter(n_pop == np) %>%
                filter(binCov == bc) %>%
                group_by("name" = meansCovconfintperm,sample_size,n_pop,binCov)%>%
                summarise(n = n())
              
              if(nrow(tempdf)==0){next}
              
              fn1 = tempdf$n[tempdf$name == "False Negative"]
              fp1 = tempdf$n[tempdf$name == "False Positive"]
              tn1 = tempdf$n[tempdf$name == "True Negative"]
              tp1 = tempdf$n[tempdf$name == "True Positive"]
              
              if(is.empty(fn1) == TRUE){fn = 0}else{fn = fn1}
              if(is.empty(fp1) == TRUE){fp = 0}else{fp = fp1}
              if(is.empty(tn1) == TRUE){tn = 0}else{tn = tn1}
              if(is.empty(tp1) == TRUE){tp = 0}else{tp = tp1}
              
              fnr = fn/(fn+tp)
              
              output1 = data.frame(sample_size = ss, 
                                   n_pop = np, 
                                   bin = bc,
                                   fnr = fnr,
                                   power = 1-fnr)
              output = rbind(output, output1)
              
            }else{
              
              tempdf <- x %>% 
                filter(sample_size == ss) %>% 
                filter(n_pop == np) %>%
                filter(binCov == bc) %>%
                group_by("name" = MeansCovconfintboot,sample_size,n_pop,binCov)%>%
                summarise(n = n())
              
              if(nrow(tempdf)==0){next}
              
              fn1 = tempdf$n[tempdf$name == "False Negative"]
              fp1 = tempdf$n[tempdf$name == "False Positive"]
              tn1 = tempdf$n[tempdf$name == "True Negative"]
              tp1 = tempdf$n[tempdf$name == "True Positive"]
              
              if(is.empty(fn1) == TRUE){fn = 0}else{fn = fn1}
              if(is.empty(fp1) == TRUE){fp = 0}else{fp = fp1}
              if(is.empty(tn1) == TRUE){tn = 0}else{tn = tn1}
              if(is.empty(tp1) == TRUE){tp = 0}else{tp = tp1}
              
              fnr = fn/(fn+tp)
              
              output1 = data.frame(sample_size = ss, 
                                   n_pop = np, 
                                   bin = bc,
                                   fnr = fnr,
                                   power = 1-fnr)
              output = rbind(output, output1)
            }
          }
        }
      }
    }else{
      
      if(resolution == "course"){
        x$binGxE = "NA"
        for(i in 1:nrow(x)){
          if(abs(x$true_GxE_means[i]) > 0 & abs(x$true_GxE_means[i]) <= 0.25){x$binGxE[i] = 0.25
          }else if(abs(x$true_GxE_means[i]) > 0.25 & abs(x$true_GxE_means[i]) <= 0.5){x$binGxE[i] = 0.5
          }else if(abs(x$true_GxE_means[i]) > 0.5 & abs(x$true_GxE_means[i]) <= 0.75){x$binGxE[i] = 0.75
          }else{x$binGxE[i] = 1}
        }
        
      }else{
        x$binGxE = "NA"
        for(i in 1:nrow(x)){
          if(x$true_GxE_means[i] == 0){x$binGxE[i] = 0
          }else if(abs(x$true_GxE_means[i]) > 0 & abs(x$true_GxE_means[i]) <= 0.15){x$binGxE[i] = 0.1
          }else if(abs(x$true_GxE_means[i]) > 0.15 & abs(x$true_GxE_means[i]) <= 0.25){x$binGxE[i] = 0.2
          }else if(abs(x$true_GxE_means[i]) > 0.25 & abs(x$true_GxE_means[i]) <= 0.35){x$binGxE[i] = 0.3
          }else if(abs(x$true_GxE_means[i]) > 0.35 & abs(x$true_GxE_means[i]) <= 0.45){x$binGxE[i] = 0.4
          }else if(abs(x$true_GxE_means[i]) > 0.45 & abs(x$true_GxE_means[i]) <= 0.55){x$binGxE[i] = 0.5
          }else if(abs(x$true_GxE_means[i]) > 0.55 & abs(x$true_GxE_means[i]) <= 0.65){x$binGxE[i] = 0.6
          }else if(abs(x$true_GxE_means[i]) > 0.65 & abs(x$true_GxE_means[i]) <= 0.75){x$binGxE[i] = 0.7
          }else if(abs(x$true_GxE_means[i]) > 0.75 & abs(x$true_GxE_means[i]) <= 0.85){x$binGxE[i] = 0.8
          }else if(abs(x$true_GxE_means[i]) > 0.85 & abs(x$true_GxE_means[i]) <= 0.95){x$binGxE[i] = 0.9
          }else{x$binGxE[i] = 1}
        }
        
      }
      
      for(i in 1:length(unique(x$sample_size))){
        for(j in 1:length(unique(x$n_pop))){
          for(k in 1:length(unique(x$binGxE))){
            
            fn1 = fn = fp1 = fp = tn1 = tn = tp1 = tp = fnr = fpr = NULL
            ss = unique(x$sample_size)[i]
            np = unique(x$n_pop)[j]
            bc = unique(x$binGxE)[k]
            
            if(analysis == "perm"){
              
              tempdf <- x %>% 
                filter(sample_size == ss) %>% 
                filter(n_pop == np) %>%
                filter(binGxE == bc) %>%
                group_by("name" = meansGxEconfintperm,sample_size,n_pop,binGxE)%>%
                summarise(n = n())
              
              if(nrow(tempdf)==0){next}
              
              fn1 = tempdf$n[tempdf$name == "False Negative"]
              fp1 = tempdf$n[tempdf$name == "False Positive"]
              tn1 = tempdf$n[tempdf$name == "True Negative"]
              tp1 = tempdf$n[tempdf$name == "True Positive"]
              
              if(is.empty(fn1) == TRUE){fn = 0}else{fn = fn1}
              if(is.empty(fp1) == TRUE){fp = 0}else{fp = fp1}
              if(is.empty(tn1) == TRUE){tn = 0}else{tn = tn1}
              if(is.empty(tp1) == TRUE){tp = 0}else{tp = tp1}
              
              fnr = fn/(fn+tp)
              
              output1 = data.frame(sample_size = ss, 
                                   n_pop = np, 
                                   bin = bc,
                                   fnr = fnr,
                                   power = 1-fnr)
              output = rbind(output, output1)
              
            }else{
              
              tempdf <- x %>% 
                filter(sample_size == ss) %>% 
                filter(n_pop == np) %>%
                filter(binGxE == bc) %>%
                group_by("name" = MeanGxEconfintboot,sample_size,n_pop,binGxE)%>%
                summarise(n = n())
              
              if(nrow(tempdf)==0){next}
              
              fn1 = tempdf$n[tempdf$name == "False Negative"]
              fp1 = tempdf$n[tempdf$name == "False Positive"]
              tn1 = tempdf$n[tempdf$name == "True Negative"]
              tp1 = tempdf$n[tempdf$name == "True Positive"]
              
              if(is.empty(fn1) == TRUE){fn = 0}else{fn = fn1}
              if(is.empty(fp1) == TRUE){fp = 0}else{fp = fp1}
              if(is.empty(tn1) == TRUE){tn = 0}else{tn = tn1}
              if(is.empty(tp1) == TRUE){tp = 0}else{tp = tp1}
              
              fnr = fn/(fn+tp)
              
              output1 = data.frame(sample_size = ss, 
                                   n_pop = np, 
                                   bin = bc,
                                   fnr = fnr,
                                   power = 1-fnr)
              output = rbind(output, output1)
              
            }
              
              
          }
        }
      }
    }
    
  }
      
  return(output)
}

is.nan.data.frame <- function(x){
  do.call(cbind, lapply(x, is.nan))}


             
      
  
