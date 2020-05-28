
Cov_matrix_meanSE_categorical <- function(genphenenv_df){
  
  # Generate data based on normal distribution of data with MEAN/STDERROR
  newdat = data.frame()
  for(i in 1:nrow(genphenenv_df)){
    n = genphenenv_df$phen_n[i]
    u = genphenenv_df$phen_data[i]
    SE = genphenenv_df$phen_mean_SE[i]
    dat = rnorm(n,u,SE)
    newdat_temp = data.frame(gen_factor = rep(genphenenv_df$gen_factor[i],n),
                             #Native_env_cat = rep(genphenenv_df$Native_env_cat[i],n),
                             nat_env_factor = rep(genphenenv_df$nat_env_factor[i],n),
                             #exp_env_cont = rep(genphenenv_df$exp_env_cont[i],n),
                             exp_env_cat = rep(genphenenv_df$exp_env_cat[i],n),
                             exp_env_factor = rep(genphenenv_df$exp_env_factor[i],n),
                             phen_n = rep(n,n),
                             phen_data = dat)
    newdat = rbind(newdat_temp,newdat)
  }
  
  # Standardize data
  phen_mean = mean(newdat$phen_data) 
  phen_sd = sd(newdat$phen_data) 
  newdat$phen_corrected = ((newdat$phen_data-phen_mean)/phen_sd)
  
  # Model Comparison
  test_temp_a = lm(phen_corrected ~ exp_env_factor + gen_factor, data = newdat)
  test_temp_b = lm(phen_corrected ~ exp_env_factor * gen_factor, data = newdat)
  result = anova(test_temp_a,test_temp_b)
  
  # Model Outputs
  if(result[[2,6]] > 0.05){
    
    test_temp = lm(phen_corrected ~ exp_env_factor + gen_factor, data = newdat)
    
    # Model diagnostics
    res = residuals(test_temp)
    shap_wilkes = shapiro.test(res)
    is.normal <- NULL
    if(shap_wilkes[[2]] > 0.05){is.normal = "Yes"}else{is.normal = "No"}
    if(length(unique(genphenenv_df$exp_env_factor))>2){ # Wont work if only 2 environments
      test_nonlinear = lm(phen_corrected ~ (exp_env_factor^2) + gen_factor, data = genphenenv_df) #Polynomial 
      lin_test = lrtest(test_temp,test_nonlinear) # If p > 0.05, then data likely non-linear
      is.linear <- NULL
      if(lin_test[[2,5]] > 0.05){is.linear = "No"}else{is.linear = "Yes"}
    }else{is.linear <- "NA"}
    
    # Extract model outputs
    lm_result = "No_GxE"
    GxE_pval = result[[2,6]]
    emm_E = emmeans(test_temp, ~ exp_env_factor)
    emm_G = emmeans(test_temp, ~ gen_factor)
    emm_GxE = NA
    E_R2 = summary(aov(test_temp))[[1]][1,2]/sum(summary(aov(test_temp))[[1]][,2])
    G_R2 = summary(aov(test_temp))[[1]][2,2]/sum(summary(aov(test_temp))[[1]][,2])
    GxE_R2 = NA
    w2_env <- (summary(aov(test_temp))[[1]][1,2]-summary(aov(test_temp))[[1]][1,1]*summary(aov(test_temp))[[1]][3,3])/
      (sum(summary(aov(test_temp))[[1]][,2])+summary(aov(test_temp))[[1]][3,3])
    w2_gen <- (summary(aov(test_temp))[[1]][2,2]-summary(aov(test_temp))[[1]][2,1]*summary(aov(test_temp))[[1]][3,3])/
      (sum(summary(aov(test_temp))[[1]][,2])+summary(aov(test_temp))[[1]][3,3])
    w2_GxE <- NA
    
    mod_dat = data.frame()
    for(i in 1:length(test_temp[[1]])){
      id = names(test_temp[[1]][i]) 
      coef = test_temp[[1]][[i]]
      lwr_CI = confint(test_temp)[i,1]
      upr_CI = confint(test_temp)[i,2]
      pval = summary(test_temp)[[4]][i,4]
      mod_dat. = data.frame("id" = id,
                            "coef" = coef,
                            "lwr_CI" = lwr_CI,
                            "upr_CI" = upr_CI,
                            "pval" = pval)
      mod_dat = rbind(mod_dat,mod_dat.)
    }
    
    # Generate Matrices for Covariance and Permutations
    cov_temp <- expand.grid(gen_factor = unique(newdat$gen_factor),
                            exp_env_factor = unique(newdat$exp_env_factor))
    cov_temp$phen_predicted = predict(test_temp,cov_temp)
  }else{
    
    test_temp = lm(phen_corrected ~ exp_env_factor * gen_factor, data = newdat)
    
    # Model diagnostics
    res = residuals(test_temp)
    shap_wilkes = shapiro.test(res)
    is.normal <- NULL
    if(shap_wilkes[[2]] > 0.05){is.normal = "Yes"}else{is.normal = "No"}
    if(length(unique(genphenenv_df$exp_env_factor))>2){ # Wont work if only 2 environments
      test_nonlinear = lm(phen_corrected ~ (exp_env_factor^2) + gen_factor, data = genphenenv_df) #Polynomial 
      lin_test = lrtest(test_temp,test_nonlinear) # If p > 0.05, then data likely non-linear
      is.linear <- NULL
      if(lin_test[[2,5]] > 0.05){is.linear = "No"}else{is.linear = "Yes"}
    }else{is.linear <- "NA"}
    
    # Extract model outputs                
    lm_result = "Yes_GxE"
    GxE_pval = result[[2,6]]
    emm_E = emmeans(test_temp, ~ exp_env_factor)
    emm_G = emmeans(test_temp, ~ gen_factor)
    emm_GxE = emmeans(test_temp, ~ exp_env_factor*gen_factor)
    E_R2 = summary(aov(test_temp))[[1]][1,2]/sum(summary(aov(test_temp))[[1]][,2])
    G_R2 = summary(aov(test_temp))[[1]][2,2]/sum(summary(aov(test_temp))[[1]][,2])
    GxE_R2 = summary(aov(test_temp))[[1]][3,2]/sum(summary(aov(test_temp))[[1]][,2])
    w2_env <- (summary(aov(test_temp))[[1]][1,2]-summary(aov(test_temp))[[1]][1,1]*summary(aov(test_temp))[[1]][4,3])/
      (sum(summary(aov(test_temp))[[1]][,2])+summary(aov(test_temp))[[1]][4,3])
    w2_gen <- (summary(aov(test_temp))[[1]][2,2]-summary(aov(test_temp))[[1]][2,1]*summary(aov(test_temp))[[1]][4,3])/
      (sum(summary(aov(test_temp))[[1]][,2])+summary(aov(test_temp))[[1]][4,3])
    w2_GxE <- (summary(aov(test_temp))[[1]][3,2]-summary(aov(test_temp))[[1]][3,1]*summary(aov(test_temp))[[1]][4,3])/
      (sum(summary(aov(test_temp))[[1]][,2])+summary(aov(test_temp))[[1]][4,3])
    
    mod_dat = data.frame()
    for(i in 1:length(test_temp[[1]])){
      id = names(test_temp[[1]][i]) 
      coef = test_temp[[1]][[i]]
      lwr_CI = confint(test_temp)[i,1]
      upr_CI = confint(test_temp)[i,2]
      pval = summary(test_temp)[[4]][i,4]
      mod_dat. = data.frame("id" = id,
                            "coef" = coef,
                            "lwr_CI" = lwr_CI,
                            "upr_CI" = upr_CI,
                            "pval" = pval)
      mod_dat = rbind(mod_dat,mod_dat.)
    }
    
    # Generate Matrices for Covariance and Permutations
    cov_temp <- expand.grid(gen_factor = unique(newdat$gen_factor),
                            exp_env_factor = unique(newdat$exp_env_factor))
    cov_temp$phen_predicted = predict(test_temp,cov_temp)
  }                 
  
  # Re-assign environmental variables
  cov_temp$nat_env_factor <- newdat$nat_env_factor[match(cov_temp$gen_factor,newdat$gen_factor)]

  # Covariance Emeans and Gmeans
  Cov_matrix = data.frame()
  
  # G_means
  G_mean_temp = data.frame()
  for(l in 1:length(unique(cov_temp$gen_factor))){
    G = unique(cov_temp$gen_factor)[l]
    G_temp = cov_temp[which(cov_temp$gen_factor == G),]
    G_mean = mean(unique(G_temp$phen_predicted))
    G_mean_temp1 = data.frame("G_mean" = G_mean,
                              "gen_factor" = G)
    G_mean_temp = rbind(G_mean_temp, G_mean_temp1)
  }
  
  # E_means
  E_mean_temp = data.frame()
  
  for(m in 1:length(unique(cov_temp$nat_env_factor))){
    E = unique(cov_temp$nat_env_factor)[m]
    E_temp = cov_temp[which(cov_temp$nat_env_factor == E),]
    E_mean = mean(E_temp$phen_predicted)
    E_mean_temp1 = data.frame("E_mean" = E_mean,
                              "native_env" = E)
    E_mean_temp = rbind(E_mean_temp, E_mean_temp1)
  }
  
  Cov_matrix = data.frame("gen" = G_mean_temp[,2],
                          "native_env" = as.factor(E_mean_temp[,2]),
                          "G_means" = G_mean_temp[,1],
                          "E_means" = E_mean_temp[,1])
  cov_est = cov(Cov_matrix$G_means,Cov_matrix$E_means) 

  model_specs = data.frame("Covariance_est" = cov_est,
                           "lm_result" = lm_result,
                           "GxE_pval" = GxE_pval,
                           "is.normal" = is.normal,
                           "is.linear" = is.linear,
                           "eta_G" = G_R2,
                           "eta_E" = E_R2,
                           "eta_GxE" = GxE_R2,
                           "w2_env" = w2_env,
                           "w2_gen" = w2_gen,
                           "w2_GxE" = w2_GxE)
  
  return(list(newdat,mod_dat,Cov_matrix,model_specs))
}
