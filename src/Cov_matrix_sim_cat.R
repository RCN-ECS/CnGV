Cov_matrix_sim_cat <- function(genphenenv_df){
  # Input:  raw, categorical, simulated data (NOT means/SE)
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
    ind_dat$phen_corrected = ((ind_dat$phen-dat_avg)/dat_std)

    # Model Comparison 
    test_temp_a = aov(phen_corrected ~ env + gen, data = ind_dat)
    test_temp_b = aov(phen_corrected ~ env * gen, data = ind_dat)
    result = anova(test_temp_a,test_temp_b)
    
    # Model Outputs
    if(result[[2,6]] > 0.05){
      
      test_temp = aov(phen_corrected ~ env + gen, data = ind_dat)
      
      # Extract model outputs
      lm_result = "No_GxE"
      GxE_pval = result[[2,6]]
      emm_E = emmeans(test_temp,"env")
      emm_G = emmeans(test_temp, "gen")
      emm_GxE = as.data.frame(emmeans(test_temp, ~ env*gen))
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
        pval = summary(test_temp)[[1]][i,5]
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
      cov_temp <- expand.grid(gen = unique(ind_dat$gen),
                              env = unique(ind_dat$env))
      cov_temp$phen_predicted = predict(test_temp,cov_temp)
    
      }else{
      
      test_temp = aov(phen_corrected ~ env * gen, data = ind_dat)
      
      # Extract model outputs                
      lm_result = "Yes_GxE"
      GxE_pval = result[[2,6]]
      G1_slope_predicted = test_temp[[1]][2]
      G2_slope_predicted = (test_temp[[1]][[2]]+test_temp[[1]][[4]])
      emm_E = as.data.frame(emmeans(test_temp, ~env))
      emm_G = as.data.frame(emmeans(test_temp, ~gen))
      emm_GxE = as.data.frame(emmeans(test_temp, ~ env*gen))
      E_R2 = summary(aov(test_temp))[[1]][1,2]/sum(summary(aov(test_temp))[[1]][,2])
      G_R2 = summary(aov(test_temp))[[1]][2,2]/sum(summary(aov(test_temp))[[1]][,2])
      GxE_R2 = summary(aov(test_temp))[[1]][3,2]/sum(summary(aov(test_temp))[[1]][,2])
      w2_env = (summary(aov(test_temp))[[1]][1,2]-(summary(aov(test_temp))[[1]][1,1]*summary(aov(test_temp))[[1]][4,3]))/
        (sum(summary(aov(test_temp))[[1]][,2])+summary(aov(test_temp))[[1]][4,3])
      w2_gen = (summary(aov(test_temp))[[1]][2,2]-(summary(aov(test_temp))[[1]][2,1]*summary(aov(test_temp))[[1]][4,3]))/
        (sum(summary(aov(test_temp))[[1]][,2])+summary(aov(test_temp))[[1]][4,3])
      w2_GxE = (summary(aov(test_temp))[[1]][3,2]-(summary(aov(test_temp))[[1]][3,1]*summary(aov(test_temp))[[1]][4,3]))/
        (sum(summary(aov(test_temp))[[1]][,2])+summary(aov(test_temp))[[1]][4,3])
      
      mod_dat. <- data.frame() 
      for(i in 1:length(test_temp[[1]])){
        id = names(test_temp[[1]][i]) 
        coef = test_temp[[1]][[i]]
        lwr_CI = confint(test_temp)[i,1]
        upr_CI = confint(test_temp)[i,2]
        pval = summary(test_temp)[[1]][i,4]
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
      cov_temp <- expand.grid(gen = unique(ind_dat$gen),
                              env = unique(ind_dat$env)) 
      cov_temp$phen_predicted = predict(test_temp,cov_temp)
    }                 
    
    # Re-assign environmental variables
    cov_temp$native_env <- ind_dat$native_env[match(cov_temp$gen,ind_dat$gen)]
    ngen <- length(unique(ind_dat$gen))
    nenv <- length(unique(ind_dat$env))
    
    # Gmeans and Emeans
    G1mean = sum(cov_temp[cov_temp$gen == "G1",-c(1,2,4)])/ngen
    G2mean = sum(cov_temp[cov_temp$gen == "G2",-c(1,2,4)])/ngen
    E1mean = sum(cov_temp[cov_temp$env == "E1",-c(1,2,4)])/nenv
    E2mean = sum(cov_temp[cov_temp$env == "E2",-c(1,2,4)])/nenv
    
    Cov_matrix = data.frame("gen" = unique(cov_temp$gen),
                            "native_env" = unique(cov_temp$env),
                            "G_means" = c(G1mean,G2mean),
                            "E_means" = c(E1mean,E2mean))
    # Covariance
    cov_est = cov(Cov_matrix$G_means,Cov_matrix$E_means)
    
    # Estimated Marginal Means (Manual)
    overall_mean <- mean(ind_dat$phen_corrected)
    if(lm_result == "Yes_GxE"){
      G1M = abs(overall_mean -
                  (mean(ind_dat$phen_corrected[ind_dat$gen=="G1"]))- # G1
                  (mean(ind_dat$phen_corrected[ind_dat$env=="E1"]))+  # E1
                  (mean(ind_dat$phen_corrected[ind_dat$env=="E1" & ind_dat$gen=="G1"]))) # G1E1
      G2M = abs(overall_mean - 
                  (mean(ind_dat$phen_corrected[ind_dat$gen=="G2"]))- # G2
                  (mean(ind_dat$phen_corrected[ind_dat$env=="E1"]))+ # E1
                  (mean(ind_dat$phen_corrected[ind_dat$env=="E1" & ind_dat$gen=="G2"]))) # G2E1
      E1M = abs(overall_mean - 
                  (mean(ind_dat$phen_corrected[ind_dat$gen=="G1"]))- # G1
                  (mean(ind_dat$phen_corrected[ind_dat$env=="E2"]))+ # E2
                  (mean(ind_dat$phen_corrected[ind_dat$env=="E2" & ind_dat$gen=="G1"]))) # G1E2
      E2M = abs(overall_mean - 
                  (mean(ind_dat$phen_corrected[ind_dat$gen=="G2"]))- # G2
                  (mean(ind_dat$phen_corrected[ind_dat$env=="E2"]))+ # E2
                  (mean(ind_dat$phen_corrected[ind_dat$env=="E2" & ind_dat$gen=="G2"]))) # G2E2
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
                       (emm_E$emmean[emm_E$env=="E1"])+ # E1
                       (emm_GxE[1,3])) # G1E1
      G2E1_emm = abs(overall_mean -
                       (emm_G$emmean[emm_G$gen=="G2"])- # G2
                       (emm_E$emmean[emm_E$env=="E1"])+ # E1
                       (emm_GxE[3,3])) # G2E1
      G1E2_emm = abs(overall_mean -
                       (emm_G$emmean[emm_G$gen=="G1"])- # G1
                       (emm_E$emmean[emm_E$env=="E2"])+ # E2
                       (emm_GxE[2,3])) # G1E2
      G2E2_emm = abs(overall_mean -
                       (emm_G$emmean[emm_G$gen=="G2"])- # G2
                       (emm_E$emmean[emm_E$env=="E2"])+ # E2
                       (emm_GxE[4,3])) # G2E2
    }else{
      G1E1_emm <- 0
      G2E1_emm <- 0
      G1E2_emm <- 0
      G2E2_emm <- 0
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

