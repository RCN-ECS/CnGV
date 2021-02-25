Cov_matrix_meanSE_continuous <- function(genphenenv_df){
  
  # Estimate slopes
  newdat = data.frame()
  newdat = meansSE_boot_cont(genphenenv_df, 20)
  
  # Data generation via multiple Regression
  phen1 = c(intercept_G1 + slope_G1[d] * env)
  plot_dat = (unique(newdat$intercept)) + # Intercept
    ((newdat[1,2])*(unique(genphenenv_df$exp_env_cat))) + #G1 estimates
    ((newdat[2,2])*(unique(genphenenv_df$exp_env_cat))) #G2 estimates
  
  for(i in 1:nrow(genphenenv_df)){
    n = genphenenv_df$phen_n[i]
    u = genphenenv_df$phen_data[i]
    SE = genphenenv_df$phen_mean_SE[i]
    dat = rnorm(n,u,SE)
    newdat_temp = data.frame(gen_factor = rep(genphenenv_df$gen_factor[i],n),
                             nat_env_mean = rep(genphenenv_df$nat_env_mean[i],n),
                             Native_env_cat = rep(genphenenv_df$Native_env_cat[i],n),
                             nat_env_factor = rep(genphenenv_df$nat_env_factor[i],n),
                             exp_env_cont = rep(genphenenv_df$exp_env_cont[i],n),
                             exp_env_cat = rep(genphenenv_df$exp_env_cat[i],n),
                             exp_env_factor = rep(genphenenv_df$exp_env_factor[i],n),
                             phen_n = rep(n,n),
                             phen_data = dat)
    newdat = rbind(newdat_temp,newdat)
  }
  
  # Standardize data
  phen_mean = mean(newdat$phen_data) 
  phen_sd = sd(newdat$phen_data) 
  nat_env_mean = mean(newdat$nat_env_mean) 
  nat_env_sd = sd(newdat$nat_env_mean)
  env_avg = mean(newdat$exp_env_cont) 
  env_std = sd(newdat$exp_env_cont)
  newdat$phen_corrected = ((newdat$phen_data-phen_mean)/phen_sd)
  newdat$env_corrected = ((newdat$exp_env_cont-env_avg)/env_std) 
  newdat$native_env_corrected = ((newdat$nat_env_mean-env_avg)/env_std) 
  
  # Model Comparison
  test_temp_a = lm(phen_corrected ~ env_corrected + gen_factor, data = newdat)
  test_temp_b = lm(phen_corrected ~ env_corrected * gen_factor, data = newdat)
  result = anova(test_temp_a,test_temp_b)
  
  # Model Outputs
  if(result[[2,6]] > 0.05){
    
    test_temp = lm(phen_corrected ~ env_corrected + gen_factor, data = newdat)
    
    # Model diagnostics
    res = residuals(test_temp)
    shap_wilkes = shapiro.test(res)
    is.normal <- NULL
    if(shap_wilkes[[2]] > 0.05){is.normal = "Yes"}else{is.normal = "No"}
    if(length(unique(genphenenv_df$env_corrected))>2){ # Wont work if only 2 environments
      test_nonlinear = lm(phen_corrected ~ (exp_env_factor^2) + gen_factor, data = genphenenv_df) #Polynomial 
      lin_test = lrtest(test_temp,test_nonlinear) # If p > 0.05, then data likely non-linear
      is.linear <- NULL
      if(lin_test[[2,5]] > 0.05){is.linear = "No"}else{is.linear = "Yes"}
    }else{is.linear <- "NA"}
    
    # Extract model outputs
    lm_result = "No_GxE"
    GxE_pval = result[[2,6]]
    emm_E = emmeans(test_temp, ~ env_corrected)
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
    cov_temp <- data.frame(gen_factor = rep(genphenenv_df$gen_factor,
                                            each = length(seq(from = min(genphenenv_df$env_corrected), 
                                                              to = max(genphenenv_df$env_corrected),  by = 0.1))),
                           env_corrected = seq(from = min(genphenenv_df$env_corrected), 
                                               to = max(genphenenv_df$env_corrected), by = 0.1))
    cov_temp$phen_predicted = predict(test_temp,cov_temp)
  }else{
    
    test_temp = lm(phen_corrected ~ env_corrected * gen_factor, data = newdat)
    
    # Model diagnostics
    res = residuals(test_temp)
    shap_wilkes = shapiro.test(res)
    is.normal <- NULL
    if(shap_wilkes[[2]] > 0.05){is.normal = "Yes"}else{is.normal = "No"}
    if(length(unique(genphenenv_df$env_corrected))>2){ # Wont work if only 2 environments
      test_nonlinear = lm(phen_corrected ~ (exp_env_factor^2) + gen_factor, data = genphenenv_df) #Polynomial 
      lin_test = lrtest(test_temp,test_nonlinear) # If p > 0.05, then data likely non-linear
      is.linear <- NULL
      if(lin_test[[2,5]] > 0.05){is.linear = "No"}else{is.linear = "Yes"}
    }else{is.linear <- "NA"}
    
    # Extract model outputs                
    lm_result = "Yes_GxE"
    GxE_pval = result[[2,6]]
    emm_E = emmeans(test_temp, ~ env_corrected)
    emm_G = emmeans(test_temp, ~ gen_factor)
    emm_GxE = emmeans(test_temp, ~ env_corrected*gen_factor)
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
    cov_temp <- data.frame(gen_factor = rep(newdat$gen_factor,
                                            each = length(seq(from = min(newdat$env_corrected), 
                                                              to = max(newdat$env_corrected),  by = 0.1))),
                           env_corrected = seq(from = min(newdat$env_corrected), 
                                               to = max(newdat$env_corrected), by = 0.1))
    cov_temp$phen_predicted = predict(test_temp,cov_temp)
  }                 
  
  # Re-assign environmental variables to model predicted data 
  cov_temp$native_env_corrected <- newdat$native_env_corrected[match(cov_temp$gen_factor,newdat$gen_factor)]
  cov_temp$env_corrected <- round(cov_temp$env_corrected, digits = 1)
  cov_temp$native_env_corrected <- round(newdat$native_env_corrected[match(cov_temp$gen_factor,newdat$gen_factor)],digits =1)
  
  # Covariance Emeans and Gmeans
  Cov_matrix = data.frame()
  
  # G_means
  E_hat = round(mean(cov_temp$env_corrected),digits = 1)
  G_means = filter(cov_temp, env_corrected == E_hat)
  
  # E_means
  E_mean_temp = data.frame()
  
  for(m in 1:length(unique(cov_temp$native_env_corrected))){
    E = unique(cov_temp$native_env_corrected)[m]
    E_temp = cov_temp[which(cov_temp$env_corrected == E),]
    E_mean = mean(unique(E_temp$phen_predicted))
    E_mean_temp1 = data.frame("E_mean" = E_mean,
                              "env" = unique(E_temp$env))
    E_mean_temp = rbind(E_mean_temp, E_mean_temp1)
  }
  
  Cov_matrix = data.frame("gen" = unique(cov_temp$gen),
                          "native_env" = as.factor(E_mean_temp[,2]),
                          "G_means" = unique(G_means$phen_predicted),
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
