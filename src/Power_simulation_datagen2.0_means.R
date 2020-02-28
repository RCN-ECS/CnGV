# Starting list of parameters
param_list <- list(
  reps = c(5),
  delta_env = c(1),#,1), # the amount the phenotype changes across 1 value of the environment (i.e., the slope). This is essentially the amount/degree of phenotypic plasticity that is the same across genotypes.
  delta_gen = c(-1),#,0,1), # the amount the phenotype changes from one genotype to the next. This is essitially the increase intercept from one genotype to the next.
  sample_size = c(5,100), 
  n_genotypes = c(3),
  n_environments = NULL,
  std_dev= c(0.01,0.5,1),#,0.25,0.5,0.75,1,1.5,2), # Random noise, with standard deviation of 1,
  interaction= c(0.5)) # this sd determines the amount of GxE)


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
#write.csv(df, file = "~/Desktop/df.csv")


ring_means <- function(param_table, n_boot){
  
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
    
    # Interaction Terms
    int <- rep(rnorm(param_table$n_genotypes[i] * n_environments, 0, sd = param_table$interaction[i]),each = param_table$sample_size[i]) # interaction term - one for each GE level

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
    
    # Calculate means and SE
    mean_df = model_df %>%
        group_by(gen_factor, exp_env_factor) %>%
        summarise(new_mean = mean(phen_corrected),
                  new_sd = sd(phen_corrected))
    mean_df$se = mean_df$new_sd/sqrt(param_table$sample_size[i])

    overall_mean <- mean(mean_df$new_mean)
    
    # Marginal Means = Each value - overall mean
    mm_df = data.frame()
    for(a in 1:nrow(mean_df)){
      value = mean_df$new_mean[a] - overall_mean
      mm_df. = data.frame("G_mm" = mean_df$gen_factor[a],
                          "E_mm" = mean_df$exp_env_factor[a],
                          value)
      mm_df = rbind(mm_df, mm_df.)
    }
    
    # GxE Magnitude
    G1_mean <- mean(mm_df$value[mm_df$G_mm == "G_1"])
    E1_mean <- mean(mm_df$value[mm_df$E_mm == "E_1"])
    GxEmag = abs(overall_mean - G1_mean - E1_mean + mean_df$new_mean[1])
  
    # Gmeans
    G_matrix = data.frame()
    for(h in 1:length(unique(mean_df$gen_factor))){
      Gname <- unique(mean_df$gen_factor)[h]
      Gmean <- mean(mean_df$new_mean[mean_df$gen_factor == Gname])
      tempdat = data.frame("G_means" = Gmean,
                           "gen_factor" = Gname)
      G_matrix = rbind(G_matrix,tempdat)
    }
    
    # Emeans  
    E_matrix = data.frame()
    for(j in 1:length(unique(mean_df$exp_env_factor))){
      Ename <- unique(mean_df$exp_env_factor)[j]
      Emean <- mean(mean_df$new_mean[mean_df$exp_env_factor == Ename])
      tempdat. = data.frame("E_means" = Emean,
                            "exp_env_factor" = Ename)
      E_matrix = rbind(E_matrix,tempdat.)
    }
    
    # Covariance
    Cov_matrix = data.frame()
    Cov_matrix <- cbind(G_matrix,E_matrix)
    cov_est = cov(Cov_matrix$G_means,Cov_matrix$E_means)
    
    #############################
    ## True Covariance and GxE ##
    #############################
    
    # Generate noiseless phenotype data 
    no_err_phen = param_table$delta_env[i] * model_df$env + param_table$delta_gen[i] * model_df$gen + model_df$int # no error

    # Standardize no error data
    dat_avg_noerror <- mean(no_err_phen) 
    dat_std_noerror <- sd(no_err_phen)
    model_df$no_err_phen_corrected <- ((no_err_phen - dat_avg_noerror)/dat_std_noerror)
    
    # ANOVA
    test_temp_noerror <- aov(no_err_phen_corrected ~ exp_env_factor * gen_factor, data = model_df)
    
    # Estimated Marginal Means
    emm_options(msg.interaction = FALSE)
    emm_E = as.data.frame(emmeans(test_temp_noerror,"exp_env_factor"))
    emm_G = as.data.frame(emmeans(test_temp_noerror, "gen_factor"))
    emm_GxE = as.data.frame(emmeans(test_temp_noerror, ~ exp_env_factor*gen_factor))
    
    # Magnitude of GxE using EMMs
    true_GxE <- abs(mean(model_df$no_err_phen_corrected) - # Overall mean
                      (emm_G$emmean[emm_G$gen_factor=="G_1"])- # G
                      (emm_E$emmean[emm_E$exp_env_factor=="E_1"])+ # E
                      (emm_GxE[1,3])) # GxE
  
    # Gmeans
    G_matrix = data.frame()
    for(y in 1:length(unique(emm_GxE$gen_factor))){
      gtemp <- filter(emm_GxE, gen_factor == unique(emm_GxE$gen_factor)[y])
      gmean <- sum(gtemp[,3])/param_table$n_genotypes[y]
      tempdat = data.frame("G_means" = gmean,
                           "gen_factor" = unique(model_df$gen_factor)[y])
      G_matrix = rbind(G_matrix,tempdat)
    }
    
    # Emeans
    E_matrix = data.frame()
    for(f in 1:length(unique(emm_GxE$exp_env_factor))){
      etemp <- filter(emm_GxE, exp_env_factor == unique(emm_GxE$exp_env_factor)[f])
      emean <- sum(etemp[,3])/n_environments
      tempdat. = data.frame("E_means" = emean,
                            "exp_env_factor" = unique(model_df$exp_env_factor)[f])
      E_matrix = rbind(E_matrix,tempdat.)
    }
    
    # Covariance
    Cov_matrix = data.frame()
    Cov_matrix <- cbind(G_matrix,E_matrix)
    true_cov = cov(Cov_matrix$G_means,Cov_matrix$E_means)
    
    ###############
    ## Bootstrap ##
    ###############
    
    # Output
    boot_df = data.frame()
    
    # Resample data 
    for(a in 1:n_boot){
      replicate_mean <- NULL
      shuffle_dat <- data.frame()
      
      for (n in 1:nlevels(mean_df$gen_factor)){
        for (j in 1:nlevels(mean_df$exp_env_factor)){
          cond_G <- filter(mean_df, gen_factor == unique(mean_df$gen_factor)[n])
          cond_E <- filter(cond_G, exp_env_factor == unique(mean_df$exp_env_factor)[j])
          
          replicate_mean <- rnorm(nrow(cond_E), mean = cond_E$new_mean, sd = cond_E$se)
          
          shuffle_dat_temp <- data.frame(gen_factor = cond_E$gen_factor,
                                         exp_env_factor = cond_E$exp_env_factor,
                                         new_mean = replicate_mean)
          shuffle_dat <- rbind(shuffle_dat, shuffle_dat_temp)
        }
      }
      
      # Bootstrap Marginal Means
      overall_mean_boot <- mean(shuffle_dat$new_mean)
      
      mm_df_boot = data.frame()
      for(a in 1:nrow(shuffle_dat)){
        value = shuffle_dat$new_mean[a] - overall_mean_boot
        mm_df_boot. = data.frame("G_mm" = shuffle_dat$gen_factor[a],
                                 "E_mm" = shuffle_dat$exp_env_factor[a],
                                 value)
        mm_df_boot = rbind(mm_df_boot, mm_df_boot.)
      }
      
      # bootstrap GxE Magnitude
      G1_mean_boot <- mean(mm_df_boot$value[mm_df_boot$G_mm == "G_1"])
      E1_mean_boot <- mean(mm_df_boot$value[mm_df_boot$E_mm == "E_1"])
      GxE_boot = abs(overall_mean_boot - G1_mean_boot - E1_mean_boot + mm_df_boot$value[1])
      
      # Bootstrap Gmeans
      G_matrix_boot = data.frame()
      for(h in 1:length(unique(shuffle_dat$gen_factor))){
        Gname <- unique(shuffle_dat$gen_factor)[h]
        Gmean <- mean(shuffle_dat$new_mean[shuffle_dat$gen_factor == Gname])
        tempdat = data.frame("G_means" = Gmean,
                             "gen_factor" = Gname)
        G_matrix_boot = rbind(G_matrix_boot,tempdat)
      }
      
      # Bootstrap Emeans  
      E_matrix_boot = data.frame()
      for(j in 1:length(unique(shuffle_dat$exp_env_factor))){
        Ename <- unique(shuffle_dat$exp_env_factor)[j]
        Emean <- mean(shuffle_dat$new_mean[shuffle_dat$exp_env_factor == Ename])
        tempdat. = data.frame("E_means" = Emean,
                              "exp_env_factor" = Ename)
        E_matrix_boot = rbind(E_matrix_boot,tempdat.)
      }
      
      # Bootstrap Covariance
      Cov_matrix_boot = data.frame()
      Cov_matrix_boot <- cbind(G_matrix_boot,E_matrix_boot)
      cov_est_boot = cov(Cov_matrix_boot$G_means,Cov_matrix_boot$E_means)
      
      boot_dat. <- data.frame("covariance" = cov_est_boot,
                              "GxE_mag" = GxE_boot)
      
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
    
    new_perm_df <- data.frame()
    
    # Shuffle dataframe
    for(b in 1:n_boot){
      
      shuffled_means = sample(mean_df$new_mean,nrow(mean_df),replace = FALSE)
      perm_df = data.frame("gen_factor" = mean_df$gen_factor, "exp_env_factor" = mean_df$exp_env_factor, shuffled_means)
      perm_df$se = mean_df$se[match(perm_df$shuffled_means,mean_df$new_mean)]
      perm_dat = perm_df
      #perm_dat <- data.frame()
      #for(g in 1:nrow(perm_df)){
       # perm_replicate_mean = rnorm(1, mean = perm_df$sampled_means[g], sd =  perm_df$se[g])
      #  perm_dat_temp <- data.frame("gen_factor" = perm_df$gen_factor[g],
       #                             "exp_env_factor" = perm_df$exp_env_factor[g],
        #                            "new_mean" = perm_replicate_mean)
        #perm_dat <- rbind(perm_dat, perm_dat_temp)
        #}
      
      # Marginal Means
      overall_mean_perm <- mean(perm_dat$shuffled_means)
      
      mm_df_perm = data.frame()
      for(a in 1:nrow(perm_dat)){
        value = perm_dat$shuffled_means[a] - overall_mean_perm
        mm_df_perm. = data.frame("G_mm" = perm_dat$gen_factor[a],
                                 "E_mm" = perm_dat$exp_env_factor[a], 
                                 value)
        mm_df_perm = rbind(mm_df_perm, mm_df_perm.)
      }
      
      # GxE Magnitude
      G1_mean_perm <- mean(mm_df_perm$value[mm_df_perm$G_mm == "G_1"])
      E1_mean_perm <- mean(mm_df_perm$value[mm_df_perm$E_mm == "E_1"])
      GxE_perm = abs(overall_mean_perm - G1_mean_perm - E1_mean_perm + mm_df_perm$value[1])
      
      # Gmeans
      G_matrix_perm = data.frame()
      for(h in 1:length(unique(perm_dat$gen_factor))){
        Gname <- unique(perm_dat$gen_factor)[h]
        Gmean <- mean(perm_dat$shuffled_means[perm_dat$gen_factor == Gname])
        tempdat = data.frame("G_means" = Gmean,
                             "gen_factor" = Gname)
        G_matrix_perm = rbind(G_matrix_perm,tempdat)
      }
      
      # Emeans
      E_matrix_perm = data.frame()
      for(j in 1:length(unique(perm_dat$exp_env_factor))){
        Ename <- unique(perm_dat$exp_env_factor)[j]
        Emean <- mean(perm_dat$shuffled_means[perm_dat$exp_env_factor == Ename])
        tempdat. = data.frame("E_means" = Emean,
                              "exp_env_factor" = Ename)
        E_matrix_perm = rbind(E_matrix_perm,tempdat.)
      }
      
      # Covariance
      Cov_matrix_perm = data.frame()
      Cov_matrix_perm <- cbind(G_matrix_perm,E_matrix_perm)
      cov_est_perm = cov(Cov_matrix_perm$G_means,Cov_matrix_perm$E_means)
      
      perm_dat. <- data.frame("covariance" = cov_est_perm,
                              "GxE_mag" = GxE_perm)
      
      new_perm_df <- rbind(new_perm_df,perm_dat.)
    }
      
    # Covariance p-value
    ptemp = (rank(c(cov_est,new_perm_df[,1]))[1])/(n_boot+1) 
    cov_pvalue = NULL
    if(ptemp < 0.5){cov_pvalue = ptemp}else{cov_pvalue = (1-ptemp)}
    
    # GxE p-value
    ptemp1 = (rank(c(GxEmag,new_perm_df[,2]))[1])/(n_boot+1) 
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
                           "GxE_estimate" = GxEmag,
                           "GxE_lwrCI" = GxE_CI[[1]],
                           "GxE_uprCI" = GxE_CI[[2]],
                           "GxE_pvalue" = GxE_pvalue)
    output = rbind(output, temp_out)
  }
  return(output)
}

test = ring_means(df,25) # Parameter table, then number of bootstraps/perms 

ggplot(test,aes(x=true_cov,y=cov_estimate, group = factor(sample_size),shape = factor(sample_size), colour = factor(std_dev)))+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin = cov_lwrCI,ymax = cov_uprCI))+
  theme_classic()+
  geom_abline(aes(slope = 1,intercept = 0))
ggplot(test,aes(x=true_GxE,y=GxE_estimate,group = factor(sample_size),shape = factor(sample_size), colour = factor(std_dev)))+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin = GxE_lwrCI,ymax = GxE_uprCI))+
  theme_classic()+
  geom_abline(aes(slope = 1,intercept = 0))

ggplot(test,aes(x=std_dev,y=mean(true_cov),group = std_dev,colour = std_dev))+geom_boxplot()#geom_point()#+#geom_abline(aes(slope = 1,intercept = 0))+theme_classic()

