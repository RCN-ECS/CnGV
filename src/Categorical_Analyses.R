

# Categorical Covariance and GxE from raw data
fun1 <- function(input_df){ 
  
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
  if(ngen == nenv){
    
    Cov_matrix = data.frame()
    Cov_matrix. = data.frame()
    Cov_matrix.. = data.frame()
    
    for(i in 1:length(unique(input_df$gen_factor))){
      gtemp <- sum(emm_GxE[emm_GxE$gen_factor == unique(input_df$gen_factor)[i],3])/ngen
      tempdat = data.frame("G_means" = gtemp,
                           "gen_factor" = unique(input_df$gen_factor)[i])
      Cov_matrix. = rbind(Cov_matrix.,tempdat)
    }
    
    for(j in 1:length(unique(input_df$exp_env_factor))){
      etemp <- sum(emm_GxE[emm_GxE$exp_env_factor == unique(input_df$exp_env_factor)[j] ,3])/nenv
      tempdat. = data.frame("E_means_unmatched" = etemp,
                            "exp_env_factor" = unique(input_df$exp_env_factor)[j])
      Cov_matrix.. = rbind(Cov_matrix..,tempdat.)
    }
    
    # Match genotypes to native environments
    Cov_matrix.$E_means = Cov_matrix..$E_means_unmatched[match(Cov_matrix..$exp_env_factor,unique(input_df$nat_env_factor_corrected))]
    Cov_matrix.$exp_env_factor = Cov_matrix..$exp_env_factor[match(Cov_matrix.$E_means,Cov_matrix..$E_means_unmatched)]
    
    # Final Output
    Cov_matrix = Cov_matrix.
    
  }else if(ngen == (2*nenv)){
    
    Cov_matrix = data.frame()
      tempdat_g2 = data.frame()
    tempdat_e = data.frame()
      tempdat_e2 = data.frame()
    
    for(k in 1:length(unique(input_df$Native_env_cat2))){
      tempdat_g = data.frame()
      subdata <- filter(input_df,Native_env_cat2 == unique(input_df$Native_env_cat2)[k])
      
      for(i in 1:length(unique(subdata$gen_factor))){
        gtemp <- sum(emm_GxE[emm_GxE$gen_factor == unique(subdata$gen_factor)[i],3])/(ngen/2)
        tempdat = data.frame("G_means" = gtemp,
                             "gen_factor" = unique(subdata$gen_factor)[i])
        tempdat_g = rbind(tempdat_g,tempdat)
      }
      tempdat_g2 = rbind(tempdat_g2,tempdat_g)
    }
      
      for(j in 1:length(unique(input_df$exp_env_factor))){
        etemp <- sum(emm_GxE[emm_GxE$exp_env_factor == unique(input_df$exp_env_factor)[j] ,3])/nenv
        tempdat. = data.frame("E_means_unmatched" = etemp,
                              "exp_env_factor" = unique(input_df$exp_env_factor)[j])
        tempdat_e = rbind(tempdat_e,tempdat.)
      }
      
      # Match genotypes to native environments
      tempdat_g2$E_means = tempdat_e$E_means_unmatched[match(tempdat_e$exp_env_factor,unique(input_df$nat_env_factor))]
      tempdat_g2$exp_env_factor = tempdat_e$exp_env_factor[match(tempdat_g2$E_means,tempdat_e$E_means_unmatched)]
      
      # Output dataframe
      Cov_matrix = tempdat_g2
  }else{
    stop("Number of genotypes does not equal number of environments")
  }
  
  # Covariance
  cov_est = cov(Cov_matrix$G_means,Cov_matrix$E_means)

  # Magnitude of GxE using EMMs
  GxE_emm <- abs(mean(input_df$phen_corrected) - # Overall mean
                (emm_G$emmean[emm_G$gen_factor=="G_1"])- # G1
                (emm_E$emmean[emm_E$exp_env_factor=="E_1"])+ # E1
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

  # Covariance
  cov_est = cov(G_vec$gmean,E_vec$emean)
  
  # Magnitude of GxE 
  GxE_emm = abs(mean(input_df$phen_corrected) - # Overall mean
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
  if(length(unique(input_df$gen_factor)) == length(unique(input_df$exp_env_factor))){
    for (i in 1:nlevels(input_df$gen_factor)){
      for (j in 1:nlevels(input_df$exp_env_factor)){
        cond_G <- filter(input_df, gen_factor == unique(input_df$gen_factor)[i])
        cond_E <- filter(cond_G, exp_env_factor == unique(input_df$exp_env_factor)[j])
        
        # Shuffle data 
        new_phen <- sample(cond_E$phen_data, size=nrow(cond_E), replace=TRUE)
            
        shuffle_dat_temp <- data.frame("gen_factor" = cond_E$gen_factor,
                                      "exp_env_factor" = cond_E$exp_env_factor,
                                      "nat_env_factor" = cond_E$nat_env_factor,
                                      "nat_env_factor_corrected" = cond_E$nat_env_factor_corrected,
                                      "phen_data" = new_phen)
        shuffle_dat <- rbind(shuffle_dat, shuffle_dat_temp)
      }
    }
  }else if(length(unique(input_df$gen_factor)) == (2*length(unique(input_df$exp_env_factor)))){
      for (i in 1:nlevels(input_df$Native_env_cat2)){
        for (j in 1:nlevels(input_df$exp_env_factor)){
          for (k in 1:nlevels(input_df$gen_factor)){
            cond_k = filter(input_df, Native_env_cat2 == unique(input_df$Native_env_cat2)[i])
            cond_E <- filter(cond_k, exp_env_factor == unique(input_df$exp_env_factor)[j])
            cond_G <- filter(cond_E, gen_factor == unique(input_df$gen_factor)[k])
          
          # Shuffle data 
          new_phen <- sample(cond_G$phen_data, size=nrow(cond_G), replace=TRUE)
          
          shuffle_dat_temp <- data.frame("gen_factor" = cond_G$gen_factor,
                                         "exp_env_factor" = cond_G$exp_env_factor,
                                         "nat_env_factor" = cond_G$nat_env_factor,
                                         "Native_env_cat" = cond_G$Native_env_cat,
                                         "Native_env_cat2" = cond_G$Native_env_cat2,
                                         "phen_data" = new_phen)
          shuffle_dat <- rbind(shuffle_dat, shuffle_dat_temp)
          }
        }
      }
  }else{
    stop("Number of genotypes does not match number of environments")
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
                                     #nat_env_factor = cond_E$nat_env_factor,
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
  
  # Each genotype and environment
  if(length(unique(input_df$gen_factor)) == length(unique(input_df$exp_env_factor))){
    # Shuffle data
    null_temp <- sample(input_df$phen_data, size=nrow(input_df), replace=FALSE)
    
    perm_dat <- data.frame("gen_factor" = input_df$gen_factor,
                           "exp_env_factor" = input_df$exp_env_factor,
                           "nat_env_factor_corrected" = input_df$nat_env_factor_corrected,
                           "nat_env_factor" = input_df$nat_env_factor,
                           "phen_data" = null_temp)
    
  }else if(length(unique(input_df$gen_factor)) == (2*length(unique(input_df$exp_env_factor)))){
    # Shuffle data
    null_temp <- sample(input_df$phen_data, size=nrow(input_df), replace=FALSE)
    
    perm_dat <- data.frame("gen_factor" = input_df$gen_factor,
                           "exp_env_factor" = input_df$exp_env_factor,
                           "nat_env_factor_corrected" = input_df$nat_env_factor_corrected,
                           "nat_env_factor" = input_df$nat_env_factor,
                           "Native_env_cat" = input_df$Native_env_cat,
                           "Native_env_cat2" = input_df$Native_env_cat2,
                           "phen_data" = null_temp)
  }else{
    stop("Number of genotypes does not match number of environments")
  }
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
                         #"nat_env_factor" = input_df$nat_env_factor,
                         "phen_data" = null_temp)
  
  new_data = fun2(perm_dat)
  return(list(new_data))
}

Categorical_sim <- function(input_df,iterations){
  
  # Load packages
  library("emmeans","lme4","tidyverse")
  
  # Number of bootstraps
  iterations = iterations
  
  # Output dataframe
  results = data.frame()
  
  # Add native environment information
  input_df$nat_env_factor_corrected <- NULL
  for(i in 1:nrow(input_df)){
  if(input_df$gen_factor[i] == "G_1") {input_df$nat_env_factor_corrected[i] = "E_1"}else{input_df$nat_env_factor_corrected[i] = "E_2"}
  }
  
  # Indexing
  for(i in 1:length(unique(input_df$index))){
    df_temp <- filter(input_df, index == unique(index)[i])
                      
  if(input_df$source == "sim"){
    new_data <- fun1(df_temp) # Estimate GxE and Covariance
    bootdat <- replicate(iterations, bootstrap_categorical(df_temp), simplify=TRUE) # Bootstrap
    permdat <- replicate(iterations, permutation_categorical(df_temp), simplify=TRUE)   # Permutation 
  }else{
    new_data <- fun2(df_temp) # Estimate GxE and Covariance
    bootdat <- replicate(iterations, bootstrap_mean(df_temp), simplify=TRUE) # Bootstrap
    permdat <- replicate(iterations, permutation_mean(df_temp), simplify=TRUE)   # Permutation 
  }
  
  #new_data <- fun1(df_temp)
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
    input_df$data_type <- unique(meta_data$Raw_data_available[match(input_df$Study_ID_phenotype,meta_data$study_ID_phenotype)])
    
    # Arrange native environments 
    input_df$nat_env_factor_corrected <- gsub("N_", "E_",as.character(input_df$nat_env_factor))
    
    # Run Appropriate Function  
    new_data <- NULL
    if(input_df$data_type == "raw data"){
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

##############
# Test  Data #
##############

# Categorical Starting parameters
Catdat <- list(
  "cat_cont" = c("categorical"), 
  "intercept_G1" = 0,
  "slope_G1" = 0.5,
  "intercept_G2" = c(1,2,4),#seq(from = -5, to = 5, by = 3),
  "slope_G2" = 0.5,#seq(from = -5, to = 5, by = 3),
  "sd" = 0.5,#seq(from = 0, to = 1, by = 0.5),
  "sample_size" = 5) #seq(from = 3, to = 12, by = 4))
source("~/Documents/GitHub/CnGV/src/data_generation_function.R") # Generate data (either cat. or cont.)
source("~/Documents/GitHub/CnGV/src/sim_means_se.R") # Generate means and SE from raw sim. data (either cat. or cont.)

# Generate categorical data
cat_raw <- data.frame(data_generation(Catdat)) # raw
outdat <- replicate(10,Categorical_sim(cat_raw,10)) # Covariance and GxE on Raw data # need to pull out values for power analysis.
cat_mean <- sim_means_se(cat_raw) # generate means 
outdat2 <- Categorical_sim(cat_mean,20) # Covariance and GxE on means 
cat_new <- mean_generation_cat(cat_mean) # generate raw data from means (rnorm)
cat_new$source <- rep("sim",nrow(cat_new)) # label data type
outdat3 <- Categorical_sim(cat_new,20) # Covariance and GxE on raw data generated by means

outdat2$orig = rep("means",nrow(outdat2)) # identifier
outdat$orig = rep("raw",nrow(outdat)) # identifier
outdat3$orig = rep("new_from_means",nrow(outdat3)) # identifier

big_data = rbind(outdat,outdat2)
ggplot(big_data,aes(x=Index,y=true_covariance,group = Index,colour=orig))+geom_point()+
  theme_classic()+geom_errorbar(ymin=big_data$Cov_lowCI,ymax=big_data$Cov_highCI)
ggplot(big_data,aes(x=Index,y=GxE_magnitude,group = Index,colour=orig))+geom_point()+
  theme_classic()+geom_errorbar(ymin=big_data$GxE_lowCI,ymax=big_data$GxE_highCI)

ggplot(big_data, aes(x = true_covariance,y = GxE_magnitude, group = Index, colour = orig))+
  geom_point()+theme_classic() #+geom_errorbar(ymin=big_data$GxE_lowCI,ymax=big_data$GxE_highCI)+geom_errorbar(ymin=big_data$Cov_lowCI,ymax=big_data$Cov_highCI)

# Test Meta-Analysis Data
test1a = Categorical_meta(test4,Extraction_Initialize,50) # Works! #630_male_wing_length
test2a = Categorical_meta(test5,Extraction_Initialize,50) # Works! #401_growth_coefficient
test3a = Categorical_meta(test1,Extraction_Initialize,50) # Works! #654_second_litter_size

# Test datasets
Extraction_Initialize = read.csv("~/Desktop/Extraction_Initialize.csv")
test1 = read.csv("~/Desktop/test1.csv") #Raw but ngen!= nenv #654_second_litter_size
test2 = read.csv("~/Desktop/test2.csv") #Means with small dataset
test3 = read.csv("~/Desktop/test3.csv") #Means
test4 = read.csv("~/Desktop/test4.csv") #Raw #401_growth_coefficient
test5 = read.csv("~/Desktop/test5.csv") # Raw but ngen!= nenv #630_male_wing_length

## Next Steps: 
#3. automate for meta analysis data 
#4. plot significant results
power_data = read.csv("~/Desktop/power_analysis_data_20200102.csv")
power_analysis <- Categorical_sim(power_data,1) # Covariance and GxE on Raw data
