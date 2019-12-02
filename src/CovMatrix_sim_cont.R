
Cov_matrix_sim_cont <- function(genphenenv_df){
  # Input:  raw, continuous, simulated data (NOT means/SE)
  # Output: lm model summaries, covariance estimates, and magnitude of GxE 
  
  # Output Dataframes
  mod_dat = data.frame()
  Cov_matrix = data.frame()
  model_specs = data.frame()
  
  # Splits into different dataframes
  for(x in 1:length(unique(genphenenv_df$index))){
    id = unique(genphenenv_df$index)[x]
    ind_dat = genphenenv_df[genphenenv_df$index == id,]
    
    # Standardize data
    dat_avg = mean(ind_dat$phen) 
    dat_std = sd(ind_dat$phen)
    env_avg = mean(ind_dat$env) 
    env_std = sd(ind_dat$env)
    ind_dat$phen_corrected = ((ind_dat$phen-dat_avg)/dat_std)
    ind_dat$env_corrected = ((ind_dat$env-env_avg)/env_std)
    
    # Model Comparison # run more on unstandardized too
    test_temp_a = lm(phen_corrected ~ env_corrected + gen, data = ind_dat)
    test_temp_b = lm(phen_corrected ~ env_corrected * gen, data = ind_dat)
    result = anova(test_temp_a,test_temp_b)
    
    # Model Outputs
    if(result[[2,6]] > 0.05){
      
      test_temp = lm(phen_corrected ~ factor(env_corrected) + gen, data = ind_dat)
      
      # Extract model outputs
      lm_result = "No_GxE"
      GxE_pval = result[[2,6]]
      emm_E = as.data.frame(emmeans(test_temp, ~ env_corrected))
      emm_G = as.data.frame(emmeans(test_temp, ~ gen))
      emm_GxE = 0
      G1_slope_predicted = test_temp[[1]][2]
      G2_slope_predicted = (test_temp[[1]][[2]]+test_temp[[1]][[3]])
      E_R2 = summary(aov(test_temp))[[1]][1,2]/sum(summary(aov(test_temp))[[1]][,2])
      G_R2 = summary(aov(test_temp))[[1]][2,2]/sum(summary(aov(test_temp))[[1]][,2])
      GxE_R2 = 0
      w2_env <- (summary(aov(test_temp))[[1]][1,2]-(summary(aov(test_temp))[[1]][1,1]*summary(aov(test_temp))[[1]][3,3]))/
        (sum(summary(aov(test_temp))[[1]][,2])+summary(aov(test_temp))[[1]][3,3])
      w2_gen <- (summary(aov(test_temp))[[1]][2,2]-summary(aov(test_temp))[[1]][2,1]*summary(aov(test_temp))[[1]][3,3])/
        (sum(summary(aov(test_temp))[[1]][,2])+summary(aov(test_temp))[[1]][3,3])
      w2_GxE <- 0
      
      mod_dat. <- data.frame() 
      for(i in 1:length(test_temp[[1]])){
        id = names(test_temp[[1]][i]) 
        coef = test_temp[[1]][[i]]
        lwr_CI = confint(test_temp)[i,1]
        upr_CI = confint(test_temp)[i,2]
        pval = summary(test_temp)[[4]][i,4]
        mod_dat.. = data.frame("index" = unique(ind_dat$index),
                              "id" = id,
                              "coef" = coef,
                              "lwr_CI" = lwr_CI,
                              "upr_CI" = upr_CI,
                              "pval" = pval)
        mod_dat. = rbind(mod_dat.,mod_dat..)  
      }
      mod_dat = rbind(mod_dat,mod_dat.)
      
      # Generate Matrices for Covariance and Permutations
      cov_temp <- data.frame(gen = rep(unique(ind_dat$gen),
                                       each = length(seq(from = min(ind_dat$env_corrected), 
                                                         to = (max(ind_dat$env_corrected)+0.1),by=0.1))),
                             env_corrected = seq(from = min(ind_dat$env_corrected), 
                                                 to = (max(ind_dat$env_corrected)+0.1),by=0.1))
      cov_temp$phen_predicted = predict(test_temp_a ,cov_temp)
    }else{
      
      test_temp = lm(phen_corrected ~ factor(env_corrected) * gen, data = ind_dat)
      
      # Extract model outputs                
      lm_result = "Yes_GxE"
      GxE_pval = result[[2,6]]
      slope_predicted = simple_slopes(test_temp_b)
      G1_slope_predicted = slope_predicted[4,3]
      G2_slope_predicted = slope_predicted[5,3]
      emm_E = as.data.frame(emmeans(test_temp, ~ env_corrected))
      emm_G = as.data.frame(emmeans(test_temp, ~ gen))
      emm_GxE = as.data.frame(emmeans(test_temp, ~ env_corrected*gen))
      E_R2 = summary(aov(test_temp))[[1]][1,2]/sum(summary(aov(test_temp))[[1]][,2])
      G_R2 = summary(aov(test_temp))[[1]][2,2]/sum(summary(aov(test_temp))[[1]][,2])
      GxE_R2 = summary(aov(test_temp))[[1]][3,2]/sum(summary(aov(test_temp))[[1]][,2])
      w2_env <- (summary(aov(test_temp))[[1]][1,2]-summary(aov(test_temp))[[1]][1,1]*summary(aov(test_temp))[[1]][4,3])/
        (sum(summary(aov(test_temp))[[1]][,2])+summary(aov(test_temp))[[1]][4,3])
      w2_gen <- (summary(aov(test_temp))[[1]][2,2]-summary(aov(test_temp))[[1]][2,1]*summary(aov(test_temp))[[1]][4,3])/
        (sum(summary(aov(test_temp))[[1]][,2])+summary(aov(test_temp))[[1]][4,3])
      w2_GxE <- (summary(aov(test_temp))[[1]][3,2]-summary(aov(test_temp))[[1]][3,1]*summary(aov(test_temp))[[1]][4,3])/
        (sum(summary(aov(test_temp))[[1]][,2])+summary(aov(test_temp))[[1]][4,3])
      
      mod_dat. <- data.frame() 
      for(i in 1:length(test_temp[[1]])){
        id = names(test_temp[[1]][i]) 
        coef = test_temp[[1]][[i]]
        lwr_CI = confint(test_temp)[i,1]
        upr_CI = confint(test_temp)[i,2]
        pval = summary(test_temp)[[4]][i,4]
        mod_dat.. = data.frame("index" = unique(ind_dat$index),
                               "id" = id,
                               "coef" = coef,
                               "lwr_CI" = lwr_CI,
                               "upr_CI" = upr_CI,
                               "pval" = pval)
        mod_dat. = rbind(mod_dat.,mod_dat..)
      }
      mod_dat = rbind(mod_dat,mod_dat.)
      
      # Generate Matrices for Covariance and Permutations
      cov_temp <- data.frame(gen = rep(unique(ind_dat$gen),
                                       each = length(seq(from = min(ind_dat$env_corrected), 
                                                         to = (max(ind_dat$env_corrected)+0.1),by=0.1))),
                             env_corrected = seq(from = min(ind_dat$env_corrected), 
                                                 to = (max(ind_dat$env_corrected)+0.1),by=0.1))
      cov_temp$phen_predicted = predict(test_temp_b,cov_temp)
    }                 
    
    # Re-assign environmental variables
    cov_temp$env_corrected <- round(cov_temp$env_corrected,digits = 1) 
    cov_temp$native_env <- ind_dat$native_env[match(cov_temp$gen,ind_dat$gen)]
    ind_dat$env_corrected <- round(ind_dat$env_corrected, digits = 1)
    cov_temp$env <- ind_dat$env[match(cov_temp$env_corrected, ind_dat$env_corrected)]
    ngen <- length(unique(ind_dat$gen))
    
    # G_means
    E_hat = mean(ind_dat$env)
    G_means = cov_temp[cov_temp$env_corrected == E_hat,]
    
    # E_means
    E_mean_temp = data.frame()
    
    for(m in 1:length(unique(cov_temp$native_env))){
      E = unique(cov_temp$native_env)[m]
      E_temp = cov_temp[which(cov_temp$env == E),]
      E_mean = mean(E_temp$phen_predicted)
      E_mean_temp1 = data.frame("E_mean" = E_mean,
                                "env" = unique(E_temp$env))
      E_mean_temp = rbind(E_mean_temp, E_mean_temp1)
    }
    
    Cov_matrix. = data.frame("Index" = rep(unique(ind_dat$index),ngen),
                             "gen" = unique(cov_temp$gen),
                             "native_env" = as.factor(E_mean_temp[,2]),
                             "G_means" = unique(G_means$phen_predicted),
                             "E_means" = E_mean_temp[,1])
    Cov_matrix <- rbind(Cov_matrix,Cov_matrix.)
    
    cov_est = cov(Cov_matrix.$G_means,Cov_matrix.$E_means)
    
    # Assign covariance type 
    cov_type = NULL
    if(cov_est < -0.01){cov_type = "CnGV"
    }else if{cov_est > 0.01){cov_type = "CoGV"
    }else{cov_est = "pure_GxE"}
    
    # Estimated Marginal Means (Manual)
    overall_mean <- mean(ind_dat$phen_corrected)
    if(lm_result == "Yes_GxE"){
    G1M = abs(overall_mean -
      (mean(ind_dat$phen_corrected[ind_dat$gen=="G1"]))- # G1
      (mean(ind_dat$phen_corrected[ind_dat$env=="-2"]))+  # E1
      (mean(ind_dat$phen_corrected[ind_dat$env=="-2" & ind_dat$gen=="G1"]))) #G1E1
    G2M = abs(overall_mean - 
      (mean(ind_dat$phen_corrected[ind_dat$gen=="G2"]))- # G2
      (mean(ind_dat$phen_corrected[ind_dat$env=="-2"]))+ # E1
      (mean(ind_dat$phen_corrected[ind_dat$env=="-2" & ind_dat$gen=="G2"]))) #G2E1
    E1M = abs(overall_mean - 
      (mean(ind_dat$phen_corrected[ind_dat$gen=="G1"]))- # G1
      (mean(ind_dat$phen_corrected[ind_dat$env=="2"]))+ # E2
      (mean(ind_dat$phen_corrected[ind_dat$env=="2" & ind_dat$gen=="G1"]))) #G1E2
    E2M = abs(overall_mean - 
      (mean(ind_dat$phen_corrected[ind_dat$gen=="G2"]))- # G2
      (mean(ind_dat$phen_corrected[ind_dat$env=="2"]))+ # E2
      (mean(ind_dat$phen_corrected[ind_dat$env=="2" & ind_dat$gen=="G2"]))) #G2E2
    }else{
      G1M <- 0
      G2M <- 0
      E1M <- 0
      E2M <- 0
    }
    
    # Estimated Marginal Means (Emmeans)
    if(lm_result == "Yes_GxE"){
      G1E1_emm = abs(overall_mean -
                       (emm_G$emmean[emm_G$gen=="G1"])- # G1
                       (emm_E[1,2])+ # E1
                       (emm_GxE[1,3])) # G1E1
      G2E1_emm = abs(overall_mean -
                       (emm_G$emmean[emm_G$gen=="G2"])- # G2
                       (emm_E[1,2])+ # E1
                       (emm_GxE[3,3])) # G2E1
      G1E2_emm = abs(overall_mean -
                       (emm_G$emmean[emm_G$gen=="G1"])- # G1
                       (emm_E[2,2])+ # E2
                       (emm_GxE[2,3])) # G1E2
      G2E2_emm = abs(overall_mean -
                       (emm_G$emmean[emm_G$gen=="G2"])- # G2
                       (emm_E[2,2])+ # E2
                       (emm_GxE[4,3])) # G2E2
    }else{
      G1E1_emm <- 0
      G2E1_emm <- 0
      G1E2_emm <- 0
      G2E2_emm <- 0
    }
    
    model_specs. = data.frame("Index" = unique(ind_dat$index),
                              "cov_type" = cov_type,
                              "Covariance_est" = cov_est,
                              "G1_slope" = unique(ind_dat$slope[ind_dat$gen == "G1"]),
                              "G1_slope_predicted" = G1_slope_predicted,
                              "G2_slope" = unique(ind_dat$slope[ind_dat$gen == "G2"]),
                              "G2_slope_predicted" = G2_slope_predicted,
                              "slope_diff" = unique(ind_dat$slope_diff),
                              "error" = unique(ind_dat$stdev),
                              "lm_result" = lm_result,
                              "GxE_pval" = GxE_pval,
                              "G_eta" = G_R2,
                              "E_eta" = E_R2,
                              "GxE_eta" = GxE_R2,
                              "GxE_emm"= G2E2_emm,
                              "GxE_lot" = G1M,
                              "E_omega" = w2_env,
                              "G_omega" = w2_gen,
                              "GxE_omega" = w2_GxE)
    model_specs <- rbind(model_specs,model_specs.)
    
  }
  return(list(mod_dat,Cov_matrix,model_specs))
}

