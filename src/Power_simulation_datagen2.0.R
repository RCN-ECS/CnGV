## ForLoop for Cluster

# Starting list of parameters
param_list <- list( 
  reps = c(100), # or more?
  delta_env = c(0.01,0.5,1),
  delta_gen = c(-1,0.01,1),
  sample_size = c(5,10,20), 
  n_pop = c(1), 
  n_environments = c(2,3,5,10,15),
  std_dev= c(0.75,1.5), 
  interaction = 5) # Vector LENGTH not magnitude

# Starting list of parameters
param_list <- list( 
  reps = c(5), # or more?
  delta_env = c(0.01,0.5,1),
  delta_gen = c(-1,0.01,1),
  sample_size = c(5,10),#c(5,10,20),#c(5,10,20), 
  n_pop = c(1),#c(2,3,4,5),#c(2,3,5,10,15), 
  n_environments = c(2,5,10),#c(2,3,5,7,10),
  std_dev= c(0.5,1.5),#c(0.5,1), 
  interaction = NULL) # Vector LENGTH not magnitude

# Table of parameters
parameter_generation <- function(param_list){
  
  # Basic parameters
  param_temp <- expand.grid("delta_env" = param_list$delta_env,
                            "delta_gen" = param_list$delta_gen,
                            "sample_size" = param_list$sample_size,
                            "n_gen" = param_list$n_pop,
                            "n_env" = param_list$n_environments,
                            "std_dev" = param_list$std_dev)
  # Add genotypes
  param_temp$n_pop = param_temp$n_gen*param_temp$n_env
  param_temp = param_temp[,-4] # Delete n_gen

  # Add replicates
  reps <- rep(c(1:param_list$reps), each = nrow(param_temp))
  param_temp2 <- data.frame("replicate" = reps, param_temp)
  
  # Generate Interaction term - N_levels depends on n_pop
  param_temp3 = data.frame()
  for(i in 1:length(unique(param_temp2$replicate))){
    for(j in 1:length(unique(param_temp2$n_pop))){
      sub = dplyr::filter(param_temp2,replicate==unique(param_temp2$replicate)[i])
      subsub = dplyr::filter(sub,n_pop == unique(sub$n_pop)[j])
      interaction_term = seq(from = 0, to = unique(subsub$n_pop), length.out = unique(subsub$n_pop))
      inter_data = merge(subsub,interaction_term)
      colnames(inter_data)[8]<- "interaction"
      param_temp3 = rbind(inter_data,param_temp3)
    }
  }
  
  # Assign Rows to dataframe
  row = seq(1:nrow(param_temp3))
  param_table = data.frame("row" = row, param_temp3)
  
  return(param_table)
  }
  

df = parameter_generation(param_list) 
args = df[108,]
args = df[81,] # Expect no GxE
dim(df)

#write.csv(df,"~/Desktop/df.csv",)

ring <- function(param_table, n_boot){
  start_time <- Sys.time()
  
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
   # set.seed = 999
    model_df = NULL
    
    # Set Conditional Parameters
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
    test_temp <- lm(phen_corrected ~ exp_env_factor * gen_factor, data = model_df)
    
    # Estimated Marginal Means
    emm_options(msg.interaction = FALSE)
    emm_E = emm_G = emm_GxE = NULL
    emm_E = as.data.frame(emmeans(test_temp,"exp_env_factor"))
    emm_G = as.data.frame(emmeans(test_temp, "gen_factor"))
    emm_GxE = as.data.frame(emmeans(test_temp, ~ exp_env_factor*gen_factor))
    
    # Gmeans
    G_matrix = data.frame()
    gtemp = gmean = tempdat = NULL
    for(h in 1:length(unique(emm_GxE$gen_factor))){
      gtemp <- filter(emm_GxE, gen_factor == unique(emm_GxE$gen_factor)[h])
      gmean <- sum(gtemp[,3])/length(unique(emm_GxE$gen_factor))
      tempdat = data.frame("G_means" = gmean,
                           "gen_factor" = unique(emm_GxE$gen_factor)[h])
      G_matrix = rbind(G_matrix,tempdat)
    }
    
    # Emeans
    E_matrix = data.frame()
    etemp = emean = tempdat. = NULL
    for(j in 1:length(unique(emm_GxE$exp_env_factor))){
      etemp <- filter(emm_GxE, exp_env_factor == unique(emm_GxE$exp_env_factor)[j])
      emean <- sum(etemp[,3])/n_environments
      tempdat. = data.frame("E_means" = emean,
                            "exp_env_factor" = unique(emm_GxE$exp_env_factor)[j])
      E_matrix = rbind(E_matrix,tempdat.)
    }
    
    # Covariance
    if(length(G_matrix[,1]) == 2){cov_est = cov(G_matrix$G_means,E_matrix$E_means)
    }else{cov_est = cor(G_matrix$G_means,E_matrix$E_means)}

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
    model_df$no_err_phen_corrected = (no_err_phen-mean(no_err_phen))/sd(no_err_phen)
    
    # Anova with no error
    test_temp_noerror <- lm(no_err_phen_corrected ~ exp_env_factor * gen_factor, data = model_df)
    
    # Estimated Marginal Means with no error
    emm_options(msg.interaction = FALSE)
    emm_options(quietly = TRUE) # Will get perfect fit warning
    emm_E_ne = emm_G_ne = emm_GxE_ne = NULL
    emm_E_ne = as.data.frame(emmeans(test_temp_noerror,"exp_env_factor"))
    emm_G_ne = as.data.frame(emmeans(test_temp_noerror, "gen_factor"))
    emm_GxE_ne = as.data.frame(emmeans(test_temp_noerror, ~ exp_env_factor*gen_factor))
    
    # Gmeans with no error
    G_matrix_ne = data.frame()
    gtemp_ne = gmean_ne  = tempdat_ne = NULL
    for(y in 1:length(unique(emm_GxE_ne$gen_factor))){
      gtemp_ne <- filter(emm_GxE_ne, gen_factor == unique(emm_GxE_ne$gen_factor)[y])
      gmean_ne <- sum(gtemp_ne[,3])/length(unique(emm_GxE_ne$gen_factor))
      tempdat_ne = data.frame("G_means" = gmean_ne,
                           "gen_factor" = unique(emm_GxE_ne$gen_factor)[y])
      G_matrix_ne = rbind(G_matrix_ne,tempdat_ne)
    }
    
    # Emeans with no error
    E_matrix_ne = data.frame()
    etemp_ne = emean_ne = tempdat._ne = NULL
    for(f in 1:length(unique(emm_GxE_ne$exp_env_factor))){
      etemp_ne <- filter(emm_GxE_ne, exp_env_factor == unique(emm_GxE_ne$exp_env_factor)[f])
      emean_ne <- sum(etemp_ne[,3])/n_environments
      tempdat._ne = data.frame("E_means" = emean_ne,
                                "exp_env_factor" = unique(emm_GxE_ne$exp_env_factor)[f])
      E_matrix_ne = rbind(E_matrix_ne,tempdat._ne)
    }
    
    # Covariance
    if(length(G_matrix_ne[,1]) == 2){true_cov = cov(G_matrix_ne$G_means,E_matrix_ne$E_means)
    }else{true_cov = cor(G_matrix_ne$G_means,E_matrix_ne$E_means)}

    # Magnitude of GxE with no error 
    true_GxE <- abs(mean(model_df$no_err_phen_corrected) - # Overall mean
                     (emm_G_ne$emmean[emm_G_ne$gen_factor=="G_1"])- # G
                     (emm_E_ne$emmean[emm_E_ne$exp_env_factor=="E_1"])+ # E
                     (emm_GxE_ne[1,3])) # GxE
    
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
      test_boot <- lm(phen_corrected ~ exp_env_factor * gen_factor, data = shuffle_dat)
      
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
      if(param_table$n_pop[i]==2){cov_est_boot = cov(Cov_matrix_boot$G_means,Cov_matrix_boot$E_means)
      }else{cov_est_boot = cor(Cov_matrix_boot$G_means,Cov_matrix_boot$E_means)}

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
      test_perm <- lm(phen_corrected ~ exp_env_factor * gen_factor, data = perm_dat)
      
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
      if(param_table$n_pop[i]==2){cov_est_perm = cov(Cov_matrix_perm$G_means,Cov_matrix_perm$E_means)
      }else{cov_est_perm = cor(Cov_matrix_perm$G_means,Cov_matrix_perm$E_means)}
      
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
    GxE_pvalue = 1-ptemp1 # Right-tailed

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

    end_time <- Sys.time()
    
    time_difference = end_time - start_time
  }
  return(list(output,time_difference))
}

test2 = ring(df[df$row == 831,],100) # Parameter table, then number of bootstraps/perms  

test2$col = NULL
for(i in 1:nrow(test2)){
  if(test2$cov_pvalue[i] <= 0.05 & test2$GxE_pvalue[i] <= 0.05){test2$col[i] = "red" # Both significant
  }else if(test2$cov_pvalue[i] <= 0.05 & test2$GxE_pvalue[i] > 0.05){test2$col[i] = "darkgreen" # Cov significant
  }else if(test2$cov_pvalue[i] > 0.05 & test2$GxE_pvalue[i] <= 0.05){test2$col[i] = "dodgerblue" # GxE significant
  }else{test2$col[i] = "grey"} # None significant
}



# dat_csv$GxE_estimate = 1.855347e+00

# True GxE by Covariance
require(RColorBrewer)
ggplot(test2, aes(x = true_cov, y = true_GxE, group = factor(n_pop),colour=col)) + 
  geom_point() + theme_classic() + ylim(0,2)+ xlim(-1,1)+
  xlab("True Covariance") + ylab("True GxE") +
  scale_colour_identity()+
  #scale_color_brewer(palette="Set1",name = "Sample Size") +
  facet_grid(~interaction)

p1=ggplot(test2, aes(x = cov_estimate, y = GxE_estimate, group = factor(n_pop),colour = factor(sample_size))) + 
  geom_point() + theme_classic() + 
  xlab("Covariance Estimate") + ylab("GxE Estimate") +
  scale_color_brewer(palette="Set1",name = "Sample Size")+
  ggtitle("Simulation - runif coding")

p2=ggplot(nint, aes(x = cov_estimate, y = GxE_estimate, group = factor(n_pop),colour = factor(sample_size))) + 
  geom_point() + theme_classic() + 
  xlab("Covariance Estimate") + ylab("GxE Estimate") +
  scale_color_brewer(palette="Set1",name = "Sample Size")+
  ggtitle("Full Simulation - original coding")

require(lattice)
require(gridExtra)
# really want to know trade off between n_pop and sample size.
#write.table(test,"Power_data.txt")
#write.csv(test,"Power_data.csv")



##########
## Plot ##
##########

pdata <- read.csv("~/Desktop/power_output/compiled_simdata.csv")


pdata$tickmark <- NULL
for(i in 1:nrow(pdata)){
if(pdata$cov_pvalue[i] > 0.05){pdata$tickmark[i]=0}else{pdata$tickmark[i]=1} # do 100 replicates for this power
}
# How many times is null false? Check true covariance and true GxE when noise = 0 for null expectation. (no sense in looking at power for noise = 0)

powerdata = data.frame()
for(i in 1:length(unique(pdata$index))){
  tempdata <- filter(pdata, index == unique(pdata$index)[i])
  num = sum(tempdata$tickmark)
  denom = length(tempdata$tickmark)
  power = num/denom
  powerdata. = data.frame("index" = unique(tempdata$index),
                          "power" = power)
  powerdata = rbind(powerdata,powerdata.)
  power_df = inner_join(powerdata,df,by = "index")
}

str(power_df)
ggplot(power_df,aes(x = sample_size, y = std_dev, fill = power)) + geom_tile() + scale_fill_gradient(low="white", high="blue") +
  theme_classic()

ggplot(power_df,aes(x = sample_size, y = n_genotypes, fill = power)) + geom_tile() + scale_fill_gradient(low="white", high="blue") +
  theme_classic()

ggplot(power_df,aes(x = n_genotypes, y = std_dev, fill = power)) + geom_tile() + scale_fill_gradient(low="white", high="blue") +
  theme_classic()
