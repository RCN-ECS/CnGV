
Cov_matrix_sim <- function(genphenenv_df){
  
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
      
      test_temp = lm(phen_corrected ~ env_corrected + gen, data = ind_dat)
      
      # Extract model outputs
      lm_result = "No_GxE"
      GxE_pval = result[[2,6]]
      emm_E = emmeans(test_temp, ~ env_corrected)
      emm_G = emmeans(test_temp, ~ gen)
      emm_GxE = NA
      G1_slope_predicted = test_temp[[1]][2]
      G2_slope_predicted = (test_temp[[1]][[2]]+test_temp[[1]][[3]])
      E_R2 = summary(aov(test_temp))[[1]][1,2]/sum(summary(aov(test_temp))[[1]][,2])
      G_R2 = summary(aov(test_temp))[[1]][2,2]/sum(summary(aov(test_temp))[[1]][,2])
      GxE_R2 = NA
      w2_env <- (summary(aov(test_temp))[[1]][1,2]-(summary(aov(test_temp))[[1]][1,1]*summary(aov(test_temp))[[1]][3,3]))/
        (sum(summary(aov(test_temp))[[1]][,2])+summary(aov(test_temp))[[1]][3,3])
      w2_gen <- (summary(aov(test_temp))[[1]][2,2]-summary(aov(test_temp))[[1]][2,1]*summary(aov(test_temp))[[1]][3,3])/
        (sum(summary(aov(test_temp))[[1]][,2])+summary(aov(test_temp))[[1]][3,3])
      w2_GxE <- NA
      
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
      cov_temp$phen_predicted = predict(test_temp,cov_temp)
    }else{
      
      test_temp = lm(phen_corrected ~ env_corrected * gen, data = ind_dat)
      
      # Extract model outputs                
      lm_result = "Yes_GxE"
      GxE_pval = result[[2,6]]
      slope_predicted = simple_slopes(test_temp)
      G1_slope_predicted = slope_predicted[4,3]
      G2_slope_predicted = slope_predicted[5,3]
      emm_E = emmeans(test_temp, ~ env_corrected)
      emm_G = emmeans(test_temp, ~ gen)
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
      cov_temp$phen_predicted = predict(test_temp,cov_temp)
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
    
    # Estimated Marginal Means (Manual)
    overall_mean <- mean(ind_dat$phen_corrected)
    if(lm_result == "Yes_GxE"){
    G1M = abs(overall_mean -
      (mean(ind_dat$phen_corrected[ind_dat$gen=="G1"]))- # G1
      (mean(ind_dat$phen_corrected[ind_dat$env=="-2"]))+  # E1
      (mean(ind_dat$phen_corrected[ind_dat$env=="-2" & ind_dat$gen=="G1"])))
    G2M = abs(overall_mean - 
      (mean(ind_dat$phen_corrected[ind_dat$gen=="G2"]))- # G2
      (mean(ind_dat$phen_corrected[ind_dat$env=="-2"]))+ # E1
      (mean(ind_dat$phen_corrected[ind_dat$env=="-2" & ind_dat$gen=="G2"])))
    E1M = abs(overall_mean - 
      (mean(ind_dat$phen_corrected[ind_dat$gen=="G1"]))- # G1
      (mean(ind_dat$phen_corrected[ind_dat$env=="2"]))+ # E2
      (mean(ind_dat$phen_corrected[ind_dat$env=="2" & ind_dat$gen=="G1"])))
    E2M = abs(overall_mean - 
      (mean(ind_dat$phen_corrected[ind_dat$gen=="G2"]))- # G2
      (mean(ind_dat$phen_corrected[ind_dat$env=="2"]))+ # E2
      (mean(ind_dat$phen_corrected[ind_dat$env=="2" & ind_dat$gen=="G2"])))
    }else{
      G1M <- NA
      G2M <- NA
      E1M <- NA
      E2M <- NA
    }
    
    # Estimated Marginal Means (Emmeans)
    if(lm_result == "Yes_GxE"){
    G1E1_emm = abs(overall_mean -
                (emm_GxE$emmean[emm_GxE$gen=="G1"])- # G1
                (emm_GxE$env_corrected[1])+ # E1
                (mean(emm_GxE$env_corrected[1] & emm_GxE$emmean[emm_GxE$gen=="G1"])))
    G2E1_emm = abs(overall_mean -
                    (emm_GxE$emmean[emm_GxE$gen=="G2"])- # G2
                    (emm_GxE$env_corrected[1])+ # E1
                    (mean(emm_GxE$env_corrected[1] & emm_GxE$emmean[emm_GxE$gen=="G2"])))
    G1E2_emm = abs(overall_mean -
                    (emm_GxE$emmean[emm_GxE$gen=="G1"])- # G1
                    (emm_GxE$env_corrected[2])+ # E2
                    (mean(emm_GxE$env_corrected[2] & emm_GxE$emmean[emm_GxE$gen=="G1"])))
    G2E2_emm = abs(overall_mean -
                    (emm_GxE$emmean[emm_GxE$gen=="G2"])- # G2
                    (emm_GxE$env_corrected[2])+ # E2
                    (mean(emm_GxE$env_corrected[2] & emm_GxE$emmean[emm_GxE$gen=="G2"])))
    }else{
      G1E1_emm <- NA
      G2E1_emm <- NA
      G1E2_emm <- NA
      G2E2_emm <- NA
    }
    
    model_specs. = data.frame("Index" = unique(ind_dat$index),
                              "G1_slope" = unique(ind_dat$slope[ind_dat$gen == "G1"]),
                              "G1_slope_predicted" = G1_slope_predicted,
                              "G2_slope" = unique(ind_dat$slope[ind_dat$gen == "G2"]),
                              "G2_slope_predicted" = G2_slope_predicted,
                              "slope_diff" = unique(ind_dat$slope_diff),
                              "error" = unique(ind_dat$stdev),
                              "lm_result" = lm_result,
                              "GxE_pval" = GxE_pval,
                              "Covariance_est" = cov_est,
                              "eta_G" = G_R2,
                              "eta_E" = E_R2,
                              "eta_GxE" = GxE_R2,
                              "G11_emm" = G1E1_emm,
                              "G12_emm" = G2E1_emm,
                              "G21_emm" = G2E1_emm,
                              "G22_emm"= G2E2_emm,
                              "G11_lot" = G1M,
                              "G12_lot" = G2M,
                              "G21_lot" = E1M,
                              "G22_lot" = E2M,
                              "w2_env" = w2_env,
                              "w2_gen" = w2_gen,
                              "w2_GxE" = w2_GxE)
    model_specs <- rbind(model_specs,model_specs.)
    
  }
  return(list(mod_dat,Cov_matrix,model_specs))
}

