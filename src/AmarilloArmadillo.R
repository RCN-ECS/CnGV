
# Function to analyze data from empirical work/ studies

amarillo_armadillo <- function(input_df, n_boot){ # Data, Number of bootstraps
  
  # Load packages
  library("emmeans")
  library("lme4")
  library("tidyverse")
  
  if(input_df$data_type == "raw"){
    
    # Output 
    output = data.frame()
    
    # Standardize data
    input_df$phen_corrected = (input_df$phen_data - mean(input_df$phen_data))/sd(input_df$phen_data)

    # Sanity Check 
    ggplot(input_df,aes(x=exp_env_factor,y=phen_corrected, group = gen_factor, colour=nat_env_factor))+geom_point()+geom_smooth(method="glm")+theme_classic()
    
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
    G_matrix$exp_env_factor <- input_df$nat_env_factor[match(G_matrix$gen_factor,input_df$gen_factor)]
    G_matrix$E_means <- E_matrix$E_means[match(G_matrix$exp_env_factor,E_matrix$exp_env_factor)]
    
    # Covariances
    correction_ratio = max(sd(G_matrix$E_means),sd(G_matrix$G_means))
    cov_corrected = round(cov(G_matrix$E_means, G_matrix$G_means)/(correction_ratio^2),2)
    
    # Magnitude of GxE -- EMMs
    GxE_emm <- abs(mean(input_df$phen_corrected) - # Overall mean
                     (emm_G$emmean[emm_G$gen_factor == "G_2"])- # G
                     (emm_E$emmean[emm_E$exp_env_factor == "E_2"])+ # E
                     (emm_GxE[4,3])) # GxE
    
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
      G_matrix_boot$exp_env_factor <- input_df$nat_env_factor[match(G_matrix_boot$gen_factor,input_df$gen_factor)]
      G_matrix_boot$E_means <- E_matrix_boot$E_means[match(G_matrix_boot$exp_env_factor,E_matrix_boot$exp_env_factor)]
      
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
      G_matrix_perm$exp_env_factor <- input_df$nat_env_factor[match(G_matrix_perm$gen_factor,input_df$gen_factor)]
      G_matrix_perm$E_means <- E_matrix_perm$E_means[match(G_matrix_perm$exp_env_factor,E_matrix_perm$exp_env_factor)]
      
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
    
  }else{ # MEANS
   
    # Output 
    output = data.frame()
    
    # Standardize data
    input_df$phen_corrected = (input_df$phen_data - mean(input_df$phen_data))/sd(input_df$phen_data)
    
    # Means of Means 
    E_means <- tapply(input_df$phen_corrected, input_df$exp_env_factor, mean)
    G_means <- tapply(input_df$phen_corrected, input_df$gen_factor, mean)
    
    # Match Genotype with Native environment
    G_matrix = data.frame("gen_factor" = names(G_means),
                          "G_means" = G_means)
    E_matrix = data.frame("exp_env_factor" = names(E_means),
                          "E_means" = E_means)
    G_matrix$exp_env_factor <- unique(input_df$nat_env_factor)[match(G_matrix$gen_factor,factor(unique(input_df$gen_factor)))]
    G_matrix$E_means <- E_matrix$E_means[match(E_matrix$exp_env_factor,G_matrix$exp_env_factor)]
    
    # Covariance
    means_correction = max(sd(E_means),sd(G_means))
    cov_means_corrected = round(cov(G_matrix$E_means, G_matrix$G_means)/(means_correction^2),2)
    
    # Magnitude of GxE -- Means
    GxE_means = abs(mean(input_df$phen_corrected) - # Overall mean
                      (mean(input_df$phen_corrected[input_df$gen_factor == "G_1"]))- # G1
                      (mean(input_df$phen_corrected[input_df$exp_env_factor == "E_1"]))+  # E1
                      (mean(input_df$phen_corrected[input_df$exp_env_factor == "E_1" & input_df$gen_factor == "G_1"]))) # G1E1
    
    ###############
    ## Bootstrap ##
    ###############
    
    # Output Dataframe
    boot_df = data.frame()
    
    # Resampling Loop
    for(a in 1:n_boot){

      # Resample means dataframe
      new_means <- data.frame()
    
    for (u in 1:nlevels(input_df$gen_factor)){
      for (r in 1:nlevels(input_df$exp_env_factor)){
        
        cond_G <- filter(input_df, gen_factor == unique(input_df$gen_factor)[u])
        cond_E <- filter(cond_G, exp_env_factor == unique(input_df$exp_env_factor)[r])
        
        # Create new means data 
        new_phen <- rnorm(nrow(cond_E), mean = cond_E$phen_data, sd = cond_E$phen_se)
        
        # Output
        new_mean_temp <- data.frame("gen_factor" = cond_E$gen_factor,
                                    "exp_env_factor" = cond_E$exp_env_factor,
                                    "nat_env_factor" = cond_E$nat_env_factor,
                                    "phen_data" = new_phen)
        new_means <- rbind(new_means, new_mean_temp)
      }
    }
    
    # Standardize resampled means
    new_means$phen_data_corrected = (new_means$phen_data - mean(new_means$phen_data))/sd(new_means$phen_data) 
    
    # Means on Means -- bootstrap
    E_means_shuffle <- tapply(new_means$phen_data_corrected, new_means$exp_env_factor, mean)
    G_means_shuffle <- tapply(new_means$phen_data_corrected, new_means$gen_factor, mean)
    
    # Match Genotype with Native environment
    G_matrix_boot = data.frame("gen_factor" = names(G_means_shuffle),
                               "G_means" = G_means_shuffle)
    E_matrix_boot = data.frame("exp_env_factor" = names(E_means_shuffle),
                               "E_means" = E_means_shuffle)
    G_matrix_boot$exp_env_factor <- unique(input_df$nat_env_factor)[match(G_matrix_boot$gen_factor,factor(unique(input_df$gen_factor)))]
    G_matrix_boot$E_means <- E_matrix_boot$E_means[match(E_matrix_boot$exp_env_factor,G_matrix_boot$exp_env_factor)]
    
    # Covariance - Bootstrap
    means_correction_shuffle = max(sd(E_means_shuffle),sd(G_means_shuffle)) 
    cov_means_boot = round(cov(G_matrix_boot$G_means, G_matrix_boot$E_means)/(means_correction_shuffle^2),2)
    
    # Magnitude of GxE - Bootstrap
    GxE_means_boot = abs(mean(new_means$phen_data_corrected) - # Overall mean
                        (mean(new_means$phen_data_corrected[new_means$gen_factor == "G_1"]))- # G1
                        (mean(new_means$phen_data_corrected[new_means$exp_env_factor == "E_1"]))+  # E1
                        (mean(new_means$phen_data_corrected[new_means$exp_env_factor == "E_1" & new_means$gen_factor == "G_1"]))) # G1E1
    
    # Bootstrap dataframe
    boot_dat. <- data.frame("covariance" = cov_means_boot,
                            "GxE_mag" = GxE_means_boot)
    boot_df <- rbind(boot_df,boot_dat.)
    }
    
    # 95% Confidence Intervals
    cov_CI = quantile(boot_df$covariance, probs=c(0.025, 0.975), type=1) 
    GxE_CI = quantile(boot_df$GxE_mag, probs=c(0.025, 0.975), type=1) 
    
    #################
    ## Permutation ##
    #################
    
    # Output dataframe
    perm_df <- data.frame()
    
    # Shuffling loop
    for(b in 1:n_boot){
      
      # Shuffle means data
      null_means. <- rnorm(nrow(input_df), mean = input_df$phen_data, sd = input_df$phen_se)
      null_means <- sample(null_means., size=length(null_means.), replace=FALSE)
      
      perm_means <- data.frame("gen_factor" = input_df$gen_factor,
                               "exp_env_factor" = input_df$exp_env_factor,
                               "nat_env_factor" = input_df$nat_env_factor,
                               "phen_data" = null_means)
      
      #Re-standardize
      perm_means$phen_data_corrected = (perm_means$phen_data - mean(perm_means$phen_data))/sd(perm_means$phen_data) 
      
      # Means on Means -- Permutation
      E_means_perm <- tapply(perm_means$phen_data_corrected, perm_means$exp_env_factor, mean)
      G_means_perm <- tapply(perm_means$phen_data_corrected, perm_means$gen_factor, mean)
      
      # Match Genotype with Native environment
      G_matrix_perm = data.frame("gen_factor" = names(G_means_perm),
                                 "G_means" = G_means_perm)
      E_matrix_perm = data.frame("exp_env_factor" = names(E_means_perm),
                                 "E_means" = E_means_perm)
      G_matrix_perm$exp_env_factor <- unique(input_df$nat_env_factor)[match(G_matrix_perm$gen_factor,factor(unique(input_df$gen_factor)))]
      G_matrix_perm$E_means <- E_matrix_perm$E_means[match(E_matrix_perm$exp_env_factor,G_matrix_perm$exp_env_factor)]
      
      # Covariance - Bootstrap
      means_correction_perm = max(sd(E_means_perm),sd(G_means_perm)) 
      cov_means_perm = round(cov(G_matrix_perm$G_means, G_matrix_perm$E_means)/(means_correction_perm^2),2)
      
      # Magnitude of GxE - Bootstrap
      GxE_means_perm = abs(mean(perm_means$phen_data_corrected) - # Overall mean
                             (mean(perm_means$phen_data_corrected[perm_means$gen_factor == "G_1"]))- # G1
                             (mean(perm_means$phen_data_corrected[perm_means$exp_env_factor == "E_1"]))+  # E1
                             (mean(perm_means$phen_data_corrected[perm_means$exp_env_factor == "E_1" & perm_means$gen_factor == "G_1"]))) # G1E1
      
      # Bootstrap dataframe
      perm_df. <- data.frame("covariance" = cov_means_perm,
                             "GxE_mag" = GxE_means_perm)
      perm_df <- rbind(perm_df,perm_df.)
       
    }
    
    #Sanity check the histograms
    #ggplot(perm_df,aes(x = covariance))+geom_histogram()
    
    # Covariance p-value
    ptemp1 = (rank(c(cov_means_corrected,perm_df$covariance))[1])/(n_boot+1) 
    cov_pvalue = NULL
    if(ptemp1 < 0.5){cov_pvalue = ptemp1}else{cov_pvalue = (1-ptemp1)} # 2-tailed
    
    # GxE p-value
    ptemp7 = (rank(c(GxE_means,perm_df$GxE_mag))[1])/(n_boot+1) 
    GxE_pvalue = 1-ptemp7 # Right-tailed
    
    # Output
    output = data.frame("Covariance Estimate" = cov_means_corrected,
                        "Covariance Lower CI" = cov_CI[[1]],
                        "Covariance Upper CI" = cov_CI[[2]],
                        "Covariance p-value" = cov_pvalue,
                        "GxE Estimate" = GxE_means,
                        "GxE Lower CI" = GxE_CI[[1]],
                        "GxE Upper CI" = GxE_CI[[2]],
                        "GxE p-value" = GxE_pvalue)
    return(output)
    
    } 
    }
  
MollyTest = amarillo_armadillo(ma, 1000) # Molly's raw data
GeoffTest = amarillo_armadillo(gt, 1000)

# Plots for paper
colors = c("G_1" = "#0066BB", "G_2" = "#FF6633")
molly_labels = c("G_1" = "Coastal populations", "G_2" = "Inland populations")
geoff_labels = c("G_1" = "Sheltered populations", "G_2" = "Wave-exposed populations")

mollyplot = ggplot(ma, aes(x = exp_env_factor, y = phen_data, group = gen_factor, colour = gen_factor)) + 
  geom_point(position = position_dodge(width = 0.1))+geom_smooth(aes(fill=gen_factor),method = "glm")+
  xlab("Environment")+ylab("Phenotype (age at metamorphosis in days)")+
  scale_x_discrete(breaks=c("E_1","E_2"),labels=c("Freshwater", "Saltwater"))+
  scale_colour_manual(values = colors,labels = molly_labels)+
  scale_fill_manual(values = colors,labels = molly_labels)+
  labs(col=" ",fill=" ")+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(axis.text.x = element_text(size=14,colour = "black"),
        axis.title.x = element_text(size=16,face="bold")) +
  theme(axis.text.y = element_text(size=14,colour = "black"),
        axis.title.y = element_text(size=16,face="bold")) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2)) +
          annotate("text", x = "E_1", y = 60, 
                   label = paste0("Cov = ", MollyTest$Covariance.Estimate,", p = ", round(MollyTest$Covariance.p.value,2)), size = 6, hjust = 0)+
          annotate("text", x = "E_1", y = 57, 
                   label = paste0("GxE = ",round(MollyTest$GxE.Estimate,2),", p = ",round(MollyTest$GxE.p.value,2)), size = 6, hjust = 0)
mollyplot

geoffplot = ggplot(gt, aes(x = exp_env_factor, y = phen_data, group = gen_factor, colour = gen_factor)) + 
  geom_point()+geom_line()+
  geom_errorbar(aes(ymin = (phen_data-phen_se),ymax = (phen_data+phen_se)),width = 0.1)+
  xlab("Environment")+ylab("Phenotype (Shell mass growth (mg))")+
  scale_x_discrete(breaks=c("E_1","E_2"),labels=c("Low Flow", "High Flow"))+
  scale_colour_manual(values = colors,labels = geoff_labels)+
  scale_fill_manual(values = colors,labels = geoff_labels)+
  labs(col=" ")+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(axis.text.x = element_text(size=14,colour = "black"),
        axis.title.x = element_text(size=16,face="bold")) +
  theme(axis.text.y = element_text(size=14,colour = "black"),
        axis.title.y = element_text(size=16,face="bold")) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2)) +
  annotate("text", x = "E_1", y = 14, 
           label = paste0("Cov = ", GeoffTest$Covariance.Estimate,", p = ", round(GeoffTest$Covariance.p.value,2)), size = 6, hjust = 0)+
  annotate("text", x = "E_1", y = 13, 
           label = paste0("GxE = ",round(GeoffTest$GxE.Estimate,2),", p = ",round(GeoffTest$GxE.p.value,2)), size = 6, hjust = 0)
geoffplot
                   