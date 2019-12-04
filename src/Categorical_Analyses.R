
cat_means = read.csv("~/Desktop/cat_means.csv")
Extraction_Initialize = read.csv("~/Desktop/Extraction_Initialize.csv")


# Categorical Covariance and GxE from raw data
fun1 <- function(input_df){ 
  
  # Load packages
  require("emmeans","lme4","tidyverse")
   
  # Standardize data
  dat_avg <- mean(input_df$phen_data) 
  dat_std <- sd(input_df$phen_data)
  input_df$phen_corrected <- ((input_df$phen_data - dat_avg)/dat_std)
  
  # Anovas
  test_temp_a <- aov(phen_corrected ~ exp_env_factor + gen_factor, data = input_df)
  test_temp_b <- aov(phen_corrected ~ exp_env_factor * gen_factor, data = input_df)
  result <- anova(test_temp_a,test_temp_b)
  
  # Estimated Marginal Means
  GxE_pval = result[[2,6]]
  emm_E = as.data.frame(emmeans(test_temp_b,"exp_env_factor"))
  emm_G = as.data.frame(emmeans(test_temp_b, "gen_factor"))
  emm_GxE = as.data.frame(emmeans(test_temp_b, ~ exp_env_factor*gen_factor))
  
  ngen <- length(unique(input_df$gen_factor))
  nenv <- length(unique(input_df$exp_env_factor))
  
  # Gmeans and Emeans
  G1mean <- sum(emm_GxE[emm_GxE$gen_factor == "G1",3])/ngen
  G2mean <- sum(emm_GxE[emm_GxE$gen_factor == "G2",3])/ngen
  E1mean <- sum(emm_GxE[emm_GxE$exp_env_factor == "E1",3])/nenv
  E2mean <- sum(emm_GxE[emm_GxE$exp_env_factor == "E2",3])/nenv
  
  Cov_matrix = data.frame("G_means" = c(G1mean,G2mean),
                          "E_means" = c(E1mean,E2mean))
  # Covariance
  cov_est = cov(Cov_matrix$G_means,Cov_matrix$E_means)

  # Magnitude of GxE using EMMs
  overall_mean = mean(input_df$phen_corrected)
  GxE_emm <- abs(overall_mean -
                (emm_G$emmean[emm_G$gen_factor=="G1"])- # G1
                (emm_E$emmean[emm_E$exp_env_factor=="E1"])+ # E1
                (emm_GxE[1,3])) # G1E1

  
  output <- data.frame("GxE_magnitude" = GxE_emm, 
                       "Covariance" = cov_est)
  return(output)
}

# Categorical Covariance and GxE from Means and SE
fun2 <- function(input_df){
  
  # Load packages
  require("emmeans","lme4","tidyverse")
  
  # Standardize data
  dat_avg <- mean(input_df$phen_data) 
  dat_std <- sd(input_df$phen_data)
  input_df$phen_corrected <- ((input_df$phen_data - dat_avg)/dat_std)
  
  # Emeans
  nenv <- length(unique(input_df$exp_env_factor))
  E_vec <- data.frame()
  
  for(i in 1:length(unique(input_df$exp_env_factor))){
    e <- filter(input_df, exp_env_factor == unique(exp_env_factor)[i])
    emean <- sum(e$phen_corrected)/nenv
    E_vec. <- data.frame("emean" = emean,
                         "exp_env_factor" = unique(e$exp_env_factor))
    E_vec <- rbind(E_vec,E_vec.)
  }
  
  # Gmeans
  ngen <- length(unique(input_df$gen_factor))  
  G_vec <- data.frame()
  
  for(i in 1:length(unique(input_df$gen_factor))){
    g <- filter(input_df, gen_factor == unique(gen_factor)[i])
    gmean <- sum(g$phen_corrected)/ngen
    G_vec. <- data.frame("gmean" = gmean,
                         "gen_factor" = unique(g$gen_factor))
    G_vec <- rbind(G_vec,G_vec.)
    }

  # Arrange matrix
  input_df$nat_env_factor_corrected <- gsub("N_", "E_",as.character(input_df$nat_env_factor))
  G_vec$native <- input_df$nat_env_factor_corrected[match(G_vec$gen_factor,input_df$gen_factor)]
  G_vec$emean <- E_vec$emean[match(E_vec$exp_env_factor,G_vec$native)]
  
  # Covariance
  cov_est = cov(G_vec$gmean,G_vec$emean)
  
  # Magnitude of GxE 
  overall_mean = mean(input_df$phen_corrected)
  GxE_emm = abs(overall_mean -
              (mean(input_df$phen_corrected[input_df$gen_factor == "G_1"]))- # G1
              (mean(input_df$phen_corrected[input_df$exp_env_factor=="E_1"]))+  # E1
              (mean(input_df$phen_corrected[input_df$exp_env_factor=="E_1" & input_df$gen_factor == "G_1"]))) #G1E1
  
  
  output <- data.frame("GxE_magnitude" = GxE_emm, 
                       "Covariance" = cov_est)
  return(output)
}

bootstrap_categorical <- function(input_df){
  
  # Outputs
  new_phen <- NULL
  shuffle_dat <- data.frame()
  
  # Each genotype and environment
  for (i in 1:nlevels(input_df$gen_factor)){
    for (j in 1:nlevels(input_df$exp_env_factor)){
      G = unique(input_df$gen_factor)[i]
      E = unique(input_df$exp_env_factor)[j]
      cond_G <- filter(input_df, gen_factor == G)
      cond_E <- filter(cond_G, exp_env_factor == E)
      
      # Shuffle data 
      new_phen <- sample(cond_E$phen_data, size=nrow(cond_E), replace=TRUE)
          
      shuffle_dat_temp <- data.frame(gen_factor = cond_E$gen_factor,
                                    exp_env_factor = cond_E$exp_env_factor,
                                    nat_env_factor = cond_E$nat_env_factor,
                                    phen_data = new_phen)
      shuffle_dat <- rbind(shuffle_dat, shuffle_dat_temp)
    }
  }
  
  new_data = fun1(shuffle_dat)
  return(list(new_data))
}

bootstrap_mean <- function(input_df){
  
  # Outputs
  new_phen <- NULL
  shuffle_dat <- data.frame()
  
  # Each genotype and environment
  for (i in 1:nlevels(input_df$gen_factor)){
    for (j in 1:nlevels(input_df$exp_env_factor)){
      G = unique(input_df$gen_factor)[i]
      E = unique(input_df$exp_env_factor)[j]
      cond_G <- filter(input_df, gen_factor == G)
      cond_E <- filter(cond_G, exp_env_factor == E)
      
      # Shuffle data 
      new_phen <- rnorm(nrow(cond_E), mean = cond_E$phen_data, sd =  cond_E$phen_mean_SE)
      
      shuffle_dat_temp <- data.frame(gen_factor = cond_E$gen_factor,
                                     exp_env_factor = cond_E$exp_env_factor,
                                     nat_env_factor = cond_E$nat_env_factor,
                                     phen_data = new_phen)
      shuffle_dat <- rbind(shuffle_dat, shuffle_dat_temp)
    }
  }
  new_data = fun2(shuffle_dat)
  return(list(new_data))
}

permutation_categorical <- function(input_df){
  
  # Outputs
  null_temp <- NULL
  perm_dat = data.frame()
  
  # Shuffle data
  null_temp <- sample(input_df$phen_data, size=nrow(input_df), replace=FALSE)
  
  perm_dat <- data.frame("gen_factor" = input_df$gen_factor,
                         "exp_env_factor" = input_df$exp_env_factor,
                         "nat_env_factor" = input_df$nat_env_factor,
                         "phen_data" = null_temp)
  
  new_data = fun1(perm_dat)
  return(list(new_data))
}

permutation_mean <- function(input_df){
  
  # Outputs
  null_temp <- NULL
  perm_dat = data.frame()
  
  # Shuffle data
  null_temp <- sample(input_df$phen_data, size=nrow(input_df), replace=FALSE)
  
  perm_dat <- data.frame("gen_factor" = input_df$gen_factor,
                         "exp_env_factor" = input_df$exp_env_factor,
                         "nat_env_factor" = input_df$nat_env_factor,
                         "phen_data" = null_temp)
  
  new_data = fun2(perm_dat)
  return(list(new_data))
}

Categorical_sim <- function(input_df,iterations){
  
  # Number of bootstraps
  iterations = iterations
  
  # Output dataframe
  results = data.frame()
  
  # Indexing
  for(i in 1:length(unique(input_df$index))){
  df_temp <- filter(input_df, index == unique(index)[i])
                      
  # Estimate GxE and Covariance
  new_data <- fun1(df_temp)
  true_covariance <- new_data["Covariance"]
  GxE_magnitude <- new_data["GxE_magnitude"]
  
  # Bootstrap
  bootdat <- replicate(iterations, bootstrap_categorical(df_temp), simplify=TRUE)
  
  # Confidence Intervals
  GxE_CI = quantile(unlist(lapply(bootdat, `[`, 1)), probs=c(0.025, 0.975), type=1) 
  GxE_avg = mean(as.numeric(unlist(lapply(bootdat, `[`, 1)))) 
  cov_CI = quantile(unlist(lapply(bootdat, `[`, 2)), probs=c(0.025, 0.975), type=1) 
  cov_avg = mean(as.numeric(unlist(lapply(bootdat, `[`, 2)))) 
  
  # Permutation 
  permdat <- replicate(iterations, permutation_categorical(df_temp), simplify=TRUE) 
  
  # Null distributions
  alpha = 0.05
  cov_nulldist_lwrCI = quantile(unlist(lapply(permdat, `[`, 2)), alpha/2)
  cov_nulldist_uprCI = quantile(unlist(lapply(permdat, `[`, 2)), 1-alpha/2)
  GxE_nulldist_lwrCI = quantile(unlist(lapply(permdat, `[`, 1)), alpha/2)
  GxE_nulldist_uprCI = quantile(unlist(lapply(permdat, `[`, 1)), 1-alpha/2)
  
  is.cov.sig = NA
  if((true_covariance[[1]] < cov_nulldist_lwrCI[[1]]) | (true_covariance[[1]] > cov_nulldist_uprCI[[1]]) ) {is.cov.sig = "yes"}else{is.cov.sig = "no"}
  is.GxE.sig = NA
  if((GxE_magnitude[[1]] < GxE_nulldist_lwrCI[[1]]) | (GxE_magnitude[[1]] > GxE_nulldist_uprCI[[1]]) ) {is.GxE.sig = "yes"}else{is.GxE.sig = "no"}
  
  # Covariance p-value
  obs = true_covariance[[1]]
  new_dat = data.frame(unlist(lapply(permdat, `[`, 2)))
  ptemp = (rank(c(obs,new_dat[,1]))[1])/(iterations+1) 
  cov_pvalue = NULL
  if(ptemp < 0.5){cov_pvalue = ptemp}else{cov_pvalue = (1-ptemp)}
  
  # GxE p-value
  obs1 = GxE_magnitude[[1]]
  new_dat1 = data.frame(unlist(lapply(permdat, `[`, 1)))
  ptemp1 = (rank(c(obs1,new_dat1[,1]))[1])/(iterations+1) 
  GxE_pvalue = NULL
  if(ptemp1 < 0.5){GxE_pvalue = ptemp1}else{GxE_pvalue = (1-ptemp1)}
  
  # Output
  results. = data.frame("Index" = unique(df_temp$index),
                       "true_covariance" = true_covariance[[1]], 
                       "Cov_lowCI" = cov_CI[[1]], 
                       "Cov_highCI" = cov_CI[[2]], 
                       "covariance_p_value" = cov_pvalue,
                       "GxE_magnitude" = GxE_magnitude[[1]],
                       "GxE_lowCI" = GxE_CI[[1]],
                       "GxE_highCI" = GxE_CI[[2]],
                       "GxE_pvalue" = GxE_pvalue)
  results <- rbind(results,results.)
  }
  return(results)
}

Categorical_meta <- function(input_df,meta_data,iterations){
  
    # Number of bootstraps
    iterations = iterations
    
    # Output dataframe
    results = data.frame()
    
    # Establish data type (means or raw)
    input_df$data_type <- unique(meta_data$Raw_data_available[match(input_df$Study_ID_phenotype,meta_data$Study_ID_phenotype)])
    
    # Run Appropriate Function  
    new_data <- NULL
    if(input_df$data_type == "raw"){
      new_data = fun1(input_df)
      bootdat = replicate(iterations, bootstrap_categorical(input_df), simplify=TRUE)
      permdat = replicate(iterations, permutation_categorical(input_df), simplify=TRUE) 
    }else{
      new_data = fun2(input_df)
      bootdat = replicate(iterations, bootstrap_mean(input_df), simplify=TRUE)
      permdat = replicate(iterations, permutation_mean(input_df), simplify=TRUE) 
      }
    
    # Estimate GxE and Covariance
    true_covariance <- new_data["Covariance"]
    GxE_magnitude <- new_data["GxE_magnitude"]
    
    # Confidence Intervals
    GxE_CI = quantile(unlist(lapply(bootdat, `[`, 1)), probs=c(0.025, 0.975), type=1) 
    GxE_avg = mean(as.numeric(unlist(lapply(bootdat, `[`, 1)))) 
    cov_CI = quantile(unlist(lapply(bootdat, `[`, 2)), probs=c(0.025, 0.975), type=1) 
    cov_avg = mean(as.numeric(unlist(lapply(bootdat, `[`, 2)))) 
    
    # Null distributions
    alpha = 0.05
    cov_nulldist_lwrCI = quantile(unlist(lapply(permdat, `[`, 2)), alpha/2)
    cov_nulldist_uprCI = quantile(unlist(lapply(permdat, `[`, 2)), 1-alpha/2)
    GxE_nulldist_lwrCI = quantile(unlist(lapply(permdat, `[`, 1)), alpha/2)
    GxE_nulldist_uprCI = quantile(unlist(lapply(permdat, `[`, 1)), 1-alpha/2)
    
    is.cov.sig = NA
    if((true_covariance[[1]] < cov_nulldist_lwrCI[[1]]) | (true_covariance[[1]] > cov_nulldist_uprCI[[1]]) ) {is.cov.sig = "yes"}else{is.cov.sig = "no"}
    is.GxE.sig = NA
    if((GxE_magnitude[[1]] < GxE_nulldist_lwrCI[[1]]) | (GxE_magnitude[[1]] > GxE_nulldist_uprCI[[1]]) ) {is.GxE.sig = "yes"}else{is.GxE.sig = "no"}
    
    # Covariance p-value
    obs = true_covariance[[1]]
    new_dat = data.frame(unlist(lapply(permdat, `[`, 2)))
    ptemp = (rank(c(obs,new_dat[,1]))[1])/(iterations+1) 
    cov_pvalue = NULL
    if(ptemp < 0.5){cov_pvalue = ptemp}else{cov_pvalue = (1-ptemp)}
    
    # GxE p-value
    obs1 = GxE_magnitude[[1]]
    new_dat1 = data.frame(unlist(lapply(permdat, `[`, 1)))
    ptemp1 = (rank(c(obs1,new_dat1[,1]))[1])/(iterations+1) 
    GxE_pvalue = NULL
    if(ptemp1 < 0.5){GxE_pvalue = ptemp1}else{GxE_pvalue = (1-ptemp1)}
    
    # Output
    results. = data.frame("Study_ID_phenotype" = unique(input_df$Study_ID_phenotype),
                          "true_covariance" = true_covariance[[1]], 
                          "Cov_lowCI" = cov_CI[[1]], 
                          "Cov_highCI" = cov_CI[[2]], 
                          "covariance_p_value" = cov_pvalue,
                          "GxE_magnitude" = GxE_magnitude[[1]],
                          "GxE_lowCI" = GxE_CI[[1]],
                          "GxE_highCI" = GxE_CI[[2]],
                          "GxE_pvalue" = GxE_pvalue)
    results <- rbind(results,results.)
    return(results)
}

# Test Simulated Data
outdat = Categorical_sim(cat_raw,50)

# Test Meta-Analysis Data
outdat2 = Categorical_meta(cat_means,Extraction_Initialize,50)

## Next Steps: 
#3. automate for meta analysis data 
#4. plot significant results
