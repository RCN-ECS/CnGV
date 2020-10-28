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

df.foundations2 <- function(delta_env, delta_gen, sample_size, n_env, std_dev, n_pop, interaction, seed1, seed2, seed3, errpop){
  
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
  model_df$GE_true_temp = delta_env * model_df$env + delta_gen * model_df$gen + model_df$int 
  
  set.seed = seed3
  model_df$e_pop <- rep(rnorm(n_pop, 0, sd = errpop), each = n_env*sample_size)
  
  model_df$GE_true <- model_df$GE_true_temp + model_df$e_pop
  
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

mod.GxE <- function(input_df){ # input is model_df
  
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
  
  # Magnitude of GxE -- EMMs
  GxE_emm_original<- abs(emm_GxE$emmean[emm_GxE$gen_factor == "G_1" & emm_GxE$exp_env_factor == "E_1"] - # GxE (Phenotype of ith genotype in jth environment)
                           emm_G$emmean[emm_G$gen_factor == "G_1"] - # phenotype of ith Genotype
                           emm_E$emmean[emm_E$exp_env_factor == "E_1"] + # phenotype of jth Environment
                           mean(emm_GxE$emmean)) # Overall mean phenotype
  
  # Output based on Stamps/Hadfield approach
  delta_E = ((emm_GxE$emmean[3]-emm_GxE$emmean[4])+(emm_GxE$emmean[1]-emm_GxE$emmean[2]))/2 
  delta_H = (emm_GxE$emmean[1]-emm_GxE$emmean[4])
  aov.coefs = coef(aov.test)
  
  
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
  
  return(list(Cov_matrix, GxE_emm_original, GxE_emm_loop, allGE, w2_GxE, eta_GxE, GxE_SumsSquares, mod_df, delta_E, delta_H, aov.coefs))
}

mean.GxE <- function(input_df,is.perm, seed){ # input is mean_df
  
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
    
  }else{ # Below generates null distribution for GxE means
    
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

bootstrap_means <- function(input_df, seedset){ # input is means_df
  
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
      set.seed = seedset
      new_phen. <- rnorm(nrow(cond), mean = cond$avg_phen, sd = cond$se) # generate replicate means
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

permutation_means <- function(input_df, permseed){ # means dataframe (mean_df)
  
  # Clear outputs
  perm_means <- data.frame()
  null_gen = null_env = null_means = NULL
  
  # Shuffle means data (same set.seed keeps phen and corresponding se matched)
  set.seed(permseed)
  null_means <- sample(input_df$avg_phen, size = length(input_df$avg_phen), replace = FALSE)
  
  set.seed(permseed)
  null_se <- sample(input_df$se, size = length(input_df$se), replace = FALSE)
  
  #null_means <- rnorm(length(null_means.), mean = null_means., sd = null_se) # create replicate mean
  
  perm_means <- data.frame("gen_factor" = input_df$gen_factor,
                           "exp_env_factor" = input_df$exp_env_factor,
                           "nat_env_factor" = input_df$nat_env_factor,
                           "avg_phen" = null_means,
                           "se" = null_se)
  # Restandardize
  perm_means$avg_phen_corrected = (perm_means$avg_phen - mean(perm_means$avg_phen))/sd(perm_means$avg_phen)
  
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

########### Meta-Analysis Functions ####################

amarillo_armadillo <- function(input_df, n_boot, data_type){ # Data, Number of bootstraps, data_type = "raw" or "means"
  
  # Load packages
  library("emmeans")
  library("lme4")
  library("tidyverse")
  
  if(data_type == "raw"){
    
    # Output 
    output = data.frame()
    
    # Standardize data
    input_df$phen_corrected = (input_df$phen_data - mean(input_df$phen_data))/sd(input_df$phen_data)
    
    # Rename Native Environments
    input_df$nat_env_factor = gsub("N_", "E_", input_df$nat_env_factor)
    
    # Sanity Check 
    # ggplot(input_df,aes(x=exp_env_factor,y=phen_corrected, group = gen_factor, colour=nat_env_factor))+geom_point()+geom_smooth(method="glm")+theme_classic()
    # ggplot(input_df,aes(x=exp_env_factor,y=phen_data, group = gen_factor, colour=nat_env_factor))+geom_point()+geom_smooth(method="glm")+theme_classic()
    
    # Anova model fit & GxE estimates
    m1 <- mod.GxE(input_df) # Raw phenotype dataframe
    
    # GxE Estimates
    cov_matrix <- m1[[1]]
    GxE_emm_original <- m1[[2]]
    GxE_emm <- m1[[3]]
    GxE_loop_output <- m1[[4]] # All GxE estimates from loop 
    omega2 <- m1[[5]]
    eta2 <- m1[[6]]
    GxE_SSq <- m1[[7]] 
    aov.df1 <- m1[[8]] # Anova SSq output
    aov_coefs <- m1[[11]]
    
    # Covariance Estimates
    cov_est = cov(cov_matrix$G_means,cov_matrix$E_means)
    cor_est = cor(cov_matrix$G_means,cov_matrix$E_means)
    correction_raw = max(sd(cov_matrix$E_means),sd(cov_matrix$G_means))
    cov_corrected = round(cov(cov_matrix$G_means, cov_matrix$E_means)/(correction_raw^2),2)
    
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
      #GxE_emm_original_boot <- m2[[2]]
      GxE_emm_boot <- m2[[3]]
      GxE_loop_output_boot <- m2[[4]] # GxE output 
      omega2_boot <- m2[[5]]
      eta2_boot <- m2[[6]]
      GxE_SSq_boot <- m2[[7]] 
      
      # Covariance Estimates
      cov_est_boot = cov(cov_matrix_boot$G_means,cov_matrix_boot$E_means)
      cor_est_boot = cor(cov_matrix_boot$G_means,cov_matrix_boot$E_means)
      correction_raw_boot = max(sd(cov_matrix_boot$E_means),sd(cov_matrix_boot$G_means))
      cov_corrected_boot = round(cov(cov_matrix_boot$G_means, cov_matrix_boot$E_means)/(correction_raw_boot^2),2)
      
      # Bootstrap dataframe
      boot_dat_raw <- data.frame("covariance" = cov_est_boot,
                                 "cor_est_boot" = cor_est_boot,
                                 "cov_corrected_boot" = cov_corrected_boot,
                                 #"GxE_emm_original_boot" = GxE_emm_original_boot,
                                 "GxE_emm_boot" = GxE_emm_boot,
                                 "GxE_omega_boot" = omega2_boot,
                                 "GxE_eta_boot" = eta2_boot,
                                 "GxE_SSq_boot" = GxE_SSq_boot)
      boot_df_raw <- rbind(boot_df_raw,boot_dat_raw)
    }
    
    # Covariance Confidence Intervals 
    cov_CI = quantile(boot_df_raw$covariance, probs=c(0.025, 0.975), type=1) 
    cor_CI = quantile(boot_df_raw$cor_est_boot, probs=c(0.025, 0.975), type=1) 
    cov_corrected_CI = quantile(boot_df_raw$cov_corrected_boot, probs=c(0.025, 0.975), type=1) 
    
    # GxE Confidence Intervals
    GxE_orig_CI = quantile(boot_df_raw$GxE_emm_original_boot, probs=c(0.025, 0.975), type=1) 
    GxE_emm_CI = quantile(boot_df_raw$GxE_emm_boot, probs = c(0.025, 0.975), type=1)
    GxE_omega_CI = quantile(boot_df_raw$GxE_omega_boot, probs=c(0.025, 0.975), type=1)
    GxE_eta_CI = quantile(boot_df_raw$GxE_eta_boot, probs=c(0.025,0.975), type = 1)
    GxE_SSq_CI = quantile(boot_df_raw$GxE_SSq_boot, probs = c(0.025,0.975), type = 1)
    
    
    #######################################
    #####   Permutation -- Raw Data   #####
    #######################################
    
    # Output dataframe
    perm_df_raw = perm_dat_raw = data.frame()
    
    for(i in 1:n_boot){
      
      # Resample Data
      perm_dat <- permutation_raw(input_df)
      
      # Anova model fit & GxE estimates
      m3 <- mod.GxE(perm_dat) # Insert shuffled permutation raw data frame
      
      # GxE Estimates
      cov_matrix_perm <- m3[[1]]
      #GxE_emm_original_perm <- m3[[2]]
      GxE_emm_perm <- m3[[3]]
      GxE_loop_output_perm <- m3[[4]] # GxE output 
      omega2_perm <- m3[[5]]
      eta2_perm <- m3[[6]]
      GxE_SSq_perm <- m3[[7]] 
      
      # Covariance Estimates
      cov_est_perm = cov(cov_matrix_perm$G_means,cov_matrix_perm$E_means)
      cor_est_perm = cor(cov_matrix_perm$G_means,cov_matrix_perm$E_means)
      correction_raw_perm = max(sd(cov_matrix_perm$E_means),sd(cov_matrix_perm$G_means))
      cov_corrected_perm = round(cov(cov_matrix_perm$G_means, cov_matrix_perm$E_means)/(correction_raw_perm^2),2)
      
      # Permutation dataframe
      perm_dat_raw <- data.frame("covariance_perm" = cov_est_perm,
                                 "cor_est_perm" = cor_est_perm,
                                 "cov_corrected_perm" = cov_corrected_perm,
                                 #"GxE_emm_original_perm" = GxE_emm_original_perm,
                                 "GxE_emm_perm" = GxE_emm_perm,
                                 "GxE_omega_perm" = omega2_perm,
                                 "GxE_eta_perm" = eta2_perm,
                                 "GxE_SSq_perm" = GxE_SSq_perm)
      perm_df_raw <- rbind(perm_df_raw,perm_dat_raw)
    }
    
    # Check: Permutation histogram
    # hist(round(perm_df_raw$cov_corrected_perm,3))
    
    # Covariance P-values
    cov_original_pvalue <- pvalue_fun(cov_est,perm_df_raw$covariance_perm,"twotail",n_boot)
    #cor_pvalue <- pvalue_fun(cor_est,perm_df_raw$cor_est_perm,"twotail",n_boot)
    cov_corrected_pvalue <- pvalue_fun(cov_corrected,perm_df_raw$cov_corrected_perm,"twotail",n_boot)
    
    # GxE P-values
    #GxE_emm_orig_pvalue <- pvalue_fun(GxE_emm_original,perm_df_raw$GxE_emm_original_perm,"righttail",n_boot)
    GxE_emm_pvalue <- pvalue_fun(GxE_emm,perm_df_raw$GxE_emm_perm,"righttail",n_boot)
    GxE_omega_pvalue <- pvalue_fun(omega2,perm_df_raw$GxE_omega_perm,"righttail",n_boot)
    GxE_eta_pvalue <- pvalue_fun(eta2,perm_df_raw$GxE_eta_perm,"righttail",n_boot)
    
    # Output
    output = data.frame("Covariance Estimate" = cov_corrected,
                        "Covariance Lower CI" = cov_corrected_CI[[1]],
                        "Covariance Upper CI" = cov_corrected_CI[[2]],
                        "Covariance p-value" = cov_corrected_pvalue,
                        "GxE Estimate" = GxE_emm,
                        "GxE Lower CI" = GxE_emm_CI[[1]],
                        "GxE Upper CI" = GxE_emm_CI[[2]],
                        "GxE p-value" = GxE_emm_pvalue,
                        "seed" = NA)
    return(output)
    
  }else{ 
    
    # Establish seeds
    seed = sample(1:15000000,1)
    set.seed(seed)
    sim_seeds <- round(runif(4*n_boot)*1000000) # More than enough
    seed1 = sim_seeds[1] # df.foundations
    seed2 = sim_seeds[2] # df.foundations
    seed3 = sim_seeds[3] # mean_gxe
    seed.set1 = sim_seeds[c(4:(4+n_boot))] # Bootstrap Means seeds
    seed.set2 = sim_seeds[c((5+n_boot):(5+2*n_boot))] # Permutation means set 1
    seed.set3 = sim_seeds[c((6+2*n_boot):(6+3*n_boot))] # Permutation means set 1
    
    # Output 
    output = data.frame()
    
    # Standardize data
    input_df$se = input_df$phen_se
    input_df$avg_phen = (input_df$phen_data - mean(input_df$phen_data))/sd(input_df$phen_data)
    input_df$avg_phen_corrected = (input_df$phen_data - mean(input_df$phen_data))/sd(input_df$phen_data)
    
    # Plot 
    # ggplot(input_df, aes(x = exp_env_factor, y = phen_data, group = gen_factor))+geom_point() + geom_smooth(method = "glm")
    # ggplot(input_df, aes(x = exp_env_factor, y = avg_phen_corrected, group = gen_factor))+geom_point() + geom_smooth(method = "glm")
    
    # GxE estimates
    m4 <- mean.GxE(input_df,is.perm = FALSE, seed = NA) # Insert means data frame (seed not necessary if is.perm is False)
    
    # GxE 
    Cov_mean_matrix <- m4[[1]]
    GxE_means <- m4[[2]]
    GxE_means_loop_output <- m4[[3]]
    
    # Covariance
    cov_est_means = cov(Cov_mean_matrix$G_means,Cov_mean_matrix$E_means)
    #cor_est_means = cor(Cov_mean_matrix$G_means,Cov_mean_matrix$E_means)
    means_correction = max(sd(Cov_mean_matrix$E_means),sd(Cov_mean_matrix$G_means))
    cov_means_corrected = round(cov(Cov_mean_matrix$G_means, Cov_mean_matrix$E_means)/(means_correction^2),2)
    
    ###################################
    ##### BOOTSTRAP -- MEANS DATA #####
    ###################################
    
    # Output Dataframes
    boot_df_means = boot_dat_means = data.frame()
    
    for(i in 1:n_boot){
      
      # Shuffle Data
      boot.seed = seed.set1
      shuffle_means <- bootstrap_means(input_df, boot.seed[i]) # Insert means data, Need n_boot seeds
      
      # GxE :: Covariance Matrix
      m5 <- mean.GxE(shuffle_means, is.perm = FALSE, seed = NA) # Insert shuffled up means data frame
      
      # GxE Estimates
      Cov_mean_matrix_boot <- m5[[1]]
      GxE_means_boot <- m5[[2]]
      
      # Covariance Estimates
      cov_mean_boot = cov(Cov_mean_matrix_boot$G_means,Cov_mean_matrix_boot$E_means)
      cor_mean_boot = cor(Cov_mean_matrix_boot$G_means,Cov_mean_matrix_boot$E_means)
      correction_mean_boot = max(sd(Cov_mean_matrix_boot$E_means),sd(Cov_mean_matrix_boot$G_means))
      cov_corrected_mean_boot = round(cov(Cov_mean_matrix_boot$G_means, Cov_mean_matrix_boot$E_means)/(correction_mean_boot^2),2)
      
      # Bootstrap dataframe
      boot_dat_means <- data.frame("cov_means_boot" = cov_mean_boot,
                                   "cor_mean_boot" = cor_mean_boot,
                                   "cov_corrected_mean_boot" = cov_corrected_mean_boot,
                                   "GxE_means_boot" = GxE_means_boot)
      boot_df_means <- rbind(boot_df_means,boot_dat_means)
    }
    
    # Covariance Confidence Intervals -- Means
    cov_means_CI = quantile(boot_df_means$cov_means_boot, probs=c(0.025, 0.975), na.rm = TRUE, type=1) 
    #cor_means_CI = quantile(boot_df_means$cor_mean_boot, probs=c(0.025, 0.975), na.rm = TRUE, type=1) 
    cov_corrected_means_CI = quantile(boot_df_means$cov_corrected_mean_boot, probs=c(0.025, 0.975), na.rm = TRUE, type=1) 
    
    # GxE Confidence Intervals -- Means
    GxE_means_CI = quantile(boot_df_means$GxE_means_boot, probs=c(0.025, 0.975), na.rm = TRUE, type=1) 
    
    #######################################
    #####  Permutation -- Means Data  #####
    #######################################
    
    # Output
    perm_df_means = perm_dat_means = data.frame()
    
    for(i in 1:n_boot){
      
      # Set seeds for perm_means and mean.GxE
      perm.seeds1 = seed.set2
      perm.seeds2 = seed.set3
      
      # Resample Data
      perm_means <- permutation_means(input_df,perm.seeds1[i])
      
      # GxE :: Covariance Matrix
      m6 <- mean.GxE(perm_means, is.perm = TRUE,perm.seeds2[i]) # Insert resampled mean phenotype dataframe
      
      # GxE Estimates
      Cov_mean_matrix_perm <- m6[[1]]
      GxE_means_perm <- m6[[2]]
      GxE_means_output_perm <- m6[[3]]
      
      # Covariance Estimates
      cov_mean_perm = cov(Cov_mean_matrix_perm$G_means,Cov_mean_matrix_perm$E_means)
      cor_mean_perm = cor(Cov_mean_matrix_perm$G_means,Cov_mean_matrix_perm$E_means)
      correction_mean_perm = max(sd(Cov_mean_matrix_perm$E_means),sd(Cov_mean_matrix_perm$G_means))
      cov_corrected_mean_perm = round(cov(Cov_mean_matrix_perm$G_means, Cov_mean_matrix_perm$E_means)/(correction_mean_perm^2),2)
      
      # Check: GxE Histogram
      # hist(GxE_means_output_perm)
      
      # Permutation dataframe -- Means
      perm_dat_means <- data.frame("cov_means_perm" = cov_mean_perm,
                                   "cor_mean_perm" = cor_mean_perm,
                                   "cov_corrected_mean_perm" = cov_corrected_mean_perm,
                                   "GxE_means_perm" = GxE_means_perm)
      perm_df_means <- rbind(perm_df_means,perm_dat_means)
    }
    
    # Check: Histogram
    perm_df_means$line = GxE_means
    ggplot(perm_df_means, aes(x = GxE_means_perm)) + geom_histogram() + geom_vline(aes(xintercept = line),colour = "red")
    
    # Covariance P-values
    cov_original_mean_pvalue <- pvalue_fun(cov_est_means,perm_df_means$cov_means_perm,"twotail", n_boot)
    #cor_mean_pvalue <- pvalue_fun(cor_est_means,perm_df_means$cor_mean_perm,"twotail", n_boot)
    cov_corrected_mean_pvalue <- pvalue_fun(cov_means_corrected,perm_df_means$cov_corrected_mean_perm,"twotail", n_boot)
    
    # GxE P-values
    GxE_mean_pvalue <- pvalue_fun(GxE_means,perm_df_means$GxE_means_perm,"righttail",n_boot)
    
    # Output
    output = data.frame("Covariance Estimate" = cov_means_corrected,
                        "Covariance Lower CI" = cov_corrected_means_CI[[1]],
                        "Covariance Upper CI" = cov_corrected_means_CI[[2]],
                        "Covariance p-value" = cov_corrected_mean_pvalue,
                        "GxE Estimate" = GxE_means,
                        "GxE Lower CI" = GxE_means_CI[[1]],
                        "GxE Upper CI" = GxE_means_CI[[2]],
                        "GxE p-value" = GxE_mean_pvalue,
                        "seed" = seed)
    return(output)
  } 
}
rm(amarillo_armadillo)

empty_as_na <- function(x){
  if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
  ifelse(as.character(x)!="", x, NA)
}
rm(empty_as_na)

fun <- function(x) {
  df <- read_excel(x)
}
rm(fun)
