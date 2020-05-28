checkfunc <- function(param_table, n_boot){
  
  # Load packages
  library("emmeans")
  library("lme4")
  library("tidyverse")
  library("rlist")
  
  # Output dataframe
  output <- data.frame()
  index = 0
  
  for(i in 1:nrow(param_table)){
    
    # Counter
    cat(i, "\n")
    index = index + 1
    
    # For reproducibility
    set.seed = 999
    
    n_environments = param_table$n_pop[i]
    
    # Approximate Cov(G,E)
    cov_GE_approx = param_table$delta_env[i] * param_table$delta_gen[i]
    
    # Dataframe foundations
    gen <- rep(1:param_table$n_pop[i], each = param_table$sample_size[i]*n_environments)
    env <- rep(1:n_environments, each = param_table$sample_size[i],times = n_environments) 
    
    # Random Noise
    noise <- rnorm(param_table$sample_size[i] * param_table$n_pop[i] * n_environments, 0, sd = param_table$std_dev[i]) 
    
    # Interaction Terms
    int <- rep(rnorm(param_table$n_pop[i] * n_environments, 0, sd = param_table$interaction[i]),each = param_table$sample_size[i]) # interaction term - one for each GE level
    
    # Create the model dataframe 
    model_df <- data.frame(gen, env, noise, int)
    model_df$gen_factor = factor(paste("G", model_df$gen, sep = "_"))
    model_df$exp_env_factor = factor(paste("E", model_df$env, sep = "_"))
    
    # Generate phenotype data using regression equation
    phen = param_table$delta_env[i] * model_df$env + param_table$delta_gen[i] * model_df$gen  + model_df$noise + model_df$int 
    model_df$phen = phen
    
    # Standardize data
    dat_avg <- mean(phen) 
    dat_std <- sd(phen)
    model_df$phen_corrected <- ((phen - dat_avg)/dat_std)
    
    # Compile phenotype data
    model_df$index = rep(unique(index),nrow(model_df))
    model_df$replicate = rep(unique(param_table$replicate[i]),nrow(model_df))
    model_df$delta_env = rep(unique(param_table$delta_env[i]),nrow(model_df))
    model_df$delta_gen = rep(unique(param_table$delta_gen[i]),nrow(model_df))  
    model_df$sample_size = rep(unique(param_table$sample_size[i]),nrow(model_df))
    model_df$n_pop = rep(unique(param_table$n_pop[i]),nrow(model_df))
    model_df$std_dev = rep(unique(param_table$std_dev[i]),nrow(model_df))
    model_df$interaction = rep(unique(param_table$interaction[i]),nrow(model_df))

    output = rbind(output, model_df)
  }
    return(output)
}

# Generates phenotype data only
df3 <- checkfunc(df,1)

ggplot(df3, aes(x = exp_env_factor, y = phen_corrected, colour = factor(gen_factor), group = factor(gen_factor))) + 
  geom_point()+theme_classic()+#geom_boxplot()
  xlab("Environment") + ylab("Phenotype") +
  #scale_color_brewer(palette="Set1") +
  #scale_fill_manual(values = getPalette(16))+
  facet_wrap(~index)

str(df_16gen)
View(df_16gen)
df_16gen$group = paste(df_16gen$sample_size,df_16gen$std_dev,sep="_")
df_16gen$effect = paste(df_16gen$delta_env,df_16gen$delta_gen,sep="__")

## Boxplots showing different distributions for different parameter sets
require(RColorBrewer)
df3$exp_env_factor <- factor(df3$exp_env_factor, levels=c("E_1", "E_2", "E_3","E_4","E_5","E_6","E_7","E_8","E_9","E_10","E_11","E_12","E_13","E_14","E_15","E_16"))
df3$gen_factor <- factor(df3$gen_factor, levels=c("G_1", "G_2", "G_3","G_4","G_5","G_6","G_7","G_8","G_9","G_10","G_11","G_12","G_13","G_14","G_15","G_16"))


ggplot(df4, aes(x = exp_env_factor, y = phen_corrected)) + 
  geom_boxplot() + theme_classic() + 
  xlab("Environment") + ylab("Phenotype") +
  scale_color_brewer(palette="Set1") +
  #scale_fill_manual(values = getPalette(16))+
  facet_wrap(~index,ncol = 4)


checkfunc2 <- function(param_table, n_boot){
  
  # Load packages
  library("emmeans")
  library("lme4")
  library("tidyverse")
  library("rlist")
  
  # Output dataframe
  output <- data.frame()
  perm_df1 <- data.frame()
  index = 0
  
  for(i in 1:nrow(param_table)){
    
    # Counter
    cat(i, "\n")
    index = index + 1
    
    # For reproducibility
    set.seed = 999
    
    n_environments = param_table$n_pop[i]
    
    # Approximate Cov(G,E)
    cov_GE_approx = param_table$delta_env[i] * param_table$delta_gen[i]
    
    # Dataframe foundations
    gen <- rep(1:param_table$n_pop[i], each = param_table$sample_size[i]*n_environments)
    env <- rep(1:n_environments, each = param_table$sample_size[i],times = n_environments) 
    
    # Random Noise
    noise <- rnorm(param_table$sample_size[i] * param_table$n_pop[i] * n_environments, 0, sd = param_table$std_dev[i]) 
    
    # Interaction Terms
    int <- rep(rnorm(param_table$n_pop[i] * n_environments, 0, sd = param_table$interaction[i]),each = param_table$sample_size[i]) # interaction term - one for each GE level
    
    # Create the model dataframe 
    model_df <- data.frame(gen, env, noise, int)
    model_df$gen_factor = factor(paste("G", model_df$gen, sep = "_"))
    model_df$exp_env_factor = factor(paste("E", model_df$env, sep = "_"))
    
    # Generate phenotype data using regression equation
    phen = param_table$delta_env[i] * model_df$env + param_table$delta_gen[i] * model_df$gen  + model_df$noise + model_df$int
    model_df$phen = phen
    
    # Standardize data
    model_df$phen_corrected = (phen-mean(phen))/sd(phen)
    
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
      gmean <- sum(gtemp[,3])/length(unique(emm_GxE$gen_factor))
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
    no_err_phen = param_table$delta_env[i] * model_df$env + param_table$delta_gen[i] * model_df$gen + model_df$int 
    
    # Standardize no error data
    model_df$no_err_phen_corrected <- (no_err_phen - mean(no_err_phen))/sd(no_err_phen)
    
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
      gmean <- sum(gtemp[,3])/length(unique(emm_GxE$gen_factor))
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

      perm_dat <- data.frame("index" = rep(unique(index),nrow(model_df)),
                             "replicate" = rep(unique(param_table$replicate[i]),nrow(model_df)),
                             "delta_env" = rep(unique(param_table$delta_env[i]),nrow(model_df)),
                             "delta_gen" = rep(unique(param_table$delta_gen[i]),nrow(model_df)) , 
                             "sample_size" = rep(unique(param_table$sample_size[i]),nrow(model_df)),
                             "n_pop" = rep(unique(param_table$n_pop[i]),nrow(model_df)),
                             "std_dev" = rep(unique(param_table$std_dev[i]),nrow(model_df)),
                             "interaction" = rep(unique(param_table$interaction[i]),nrow(model_df)),
                             "gen_factor" = model_df$gen_factor,
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
    GxE_pvalue = 1-ptemp1
    #if(ptemp1 < 0.5){GxE_pvalue = ptemp1}else{GxE_pvalue = (1-ptemp1)} # 1-tailed
    
    # Generate Outputs
    temp_out <- data.frame("row" = param_table$row[i],
                           "delta_env" = param_table$delta_env[i],
                           "delta_gen" = param_table$delta_gen[i],
                           "sample_size" = param_table$sample_size[i],
                           "n_pop" = param_table$n_pop[i],
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
  return(list(output,perm_df))
}
# add row with one change
#testset[4,] <- c(11,1,0,-1,5,5,0.1,1.5)
permcheck = checkfunc2(df,50)

outs = permcheck[[1]]
diag = permcheck[[2]]

#Histograms
library(gridExtra)

#Cov Hists
p <- list()
for(i in 1:nrow(outs)){
  
  covdat = diag[,2]
  
  p[[i]] <- ggplot(diag, aes(x = covariance)) + geom_histogram() + theme_classic() + 
    geom_vline(xintercept=outs$cov_estimate[i])+
    geom_vline(xintercept=outs$true_cov[i],colour = "red")+
  #annotate("text", label = paste0("GxE_est=",round(annotdat$GxE_estimate,2),": GxE_pval=",round(annotdat$GxE_pvalue,2)), size = 4, x = 0, y =1500)+
  annotate("text", label = paste0("cov_est=",round(outs$cov_estimate[i],3),": cov_pval=",round(outs$cov_pvalue[i],3)), size = 4, x = 0, y =30)+  
  ggtitle(paste0("Row =",row(outs)[i]))
    
}
do.call(grid.arrange,p)


#GxE Hists
q <- list()
for(i in 1:nrow(outs)){
  
  covdat = diag[,2]
  
  q[[i]] <- ggplot(diag, aes(x = GxE_mag)) + geom_histogram() + theme_classic() + 
    geom_vline(xintercept=outs$GxE_estimate[i])+
    geom_vline(xintercept=outs$true_GxE[i],colour = "red")+
    annotate("text", label = paste0("GxE_est=",round(outs$GxE_estimate[i],3),": GxE_pval=",round(outs$GxE_pvalue[i],3)), size = 4, x = 0.25, y =30)+  
    ggtitle(paste0("Row =",row(outs)[i]))
  
}
do.call(grid.arrange,q)

#write.csv(permcheck, "~/Desktop/permcheck.csv")
permcheck$exp_env_factor <- factor(permcheck$exp_env_factor, levels=c("E_1", "E_2", "E_3","E_4","E_5","E_6","E_7","E_8","E_9","E_10","E_11","E_12","E_13","E_14","E_15","E_16"))
permcheck$gen_factor <- factor(permcheck$gen_factor, levels=c("G_1", "G_2", "G_3","G_4","G_5","G_6","G_7","G_8","G_9","G_10","G_11","G_12","G_13","G_14","G_15","G_16"))

ggplot(permcheck, aes(x = exp_env_factor, y = phen_corrected, group = gen_factor,colour = gen_factor)) + 
  geom_point() + theme_classic() + 
  xlab("Environment") + ylab("Phenotype") +
  #scale_color_brewer(palette="Set1") +
  #scale_fill_manual(values = getPalette(16))+
  facet_grid(delta_env ~ interaction)

  

str(df_16gen)
View(df_16gen)
df_16gen$group = paste(df_16gen$sample_size,df_16gen$std_dev,sep="_")
df_16gen$effect = paste(df_16gen$delta_env,df_16gen$delta_gen,sep="__")

## Boxplots showing different distributions for different parameter sets
require(RColorBrewer)


