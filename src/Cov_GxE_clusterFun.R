
# Load packages 
library("emmeans")
emm_options(msg.interaction = FALSE)
library("lme4")
library("rlist")
library("dplyr")
library("ggplot2")

# Load functions
source("~/Documents/GitHub/CnGV/src/Cov_GxE_functions.R")

start.time <- Sys.time()
set.seed(86)
#args = commandArgs(trailingOnly = TRUE)

# Load Parameters
row <- as.numeric(args[1])
replicate <- as.numeric(args[2])
delta_env <- as.numeric(args[3])
delta_gen <- as.numeric(args[4])
sample_size <- as.numeric(args[5])
n_env <- as.numeric(args[6])
std_dev <- as.numeric(args[7])
n_pop <- as.numeric(args[8])
interaction <- as.numeric(args[9])
n_boot <- 19

# Output dataframes
output <- data.frame()
phen_out <- data.frame()
GEmeans_out <- data.frame()
model_out <- data.frame()
boot_df_raw <- data.frame()
boot_df_means <- data.frame()
perm_df_raw <- data.frame()
perm_df_means <- data.frame()

# Simulate data
dfs <- df.foundations(delta_env, delta_gen, sample_size, n_env, std_dev, n_pop, interaction)

# Working dataframes
model_df <- dfs[[1]]     # Raw data
mean_df <- dfs[[2]]      # Means data
model_df.ne <- dfs[[3]]  # Raw data with no error
mean_df.ne <- dfs[[4]]   # Means data with no error

# Tracking: Phenotype output
phen_out. <- data.frame("row" = rep(unique(row),nrow(model_df)), 
                        "replicate" = rep(unique(replicate),nrow(model_df)), 
                        "delta_env" = rep(unique(delta_env),nrow(model_df)), 
                        "delta_gen" = rep(unique(delta_gen),nrow(model_df)), 
                        "sample_size" = rep(unique(sample_size),nrow(model_df)), 
                        "n_pop" = rep(unique(n_pop),nrow(model_df)), 
                        "n_env" = rep(unique(n_env),nrow(model_df)),
                        "std_dev" = rep(unique(std_dev),nrow(model_df)), 
                        "interaction" = rep(unique(interaction),nrow(model_df))) 
phen_out <- cbind(phen_out.,model_df)

# Check: Raw Phenotype 
ggplot(model_df, aes(x = exp_env_factor, y = phen_corrected, group = gen_factor, fill = nat_env_factor,colour = nat_env_factor)) + geom_point() + geom_smooth() + theme_classic()

# Check: Mean Phenotype 
ggplot(mean_df, aes(x = exp_env_factor, y = avg_phen_corrected, group = gen_factor, fill = nat_env_factor,colour = nat_env_factor)) + geom_point() + geom_smooth() + theme_classic()

# Check: Raw Phenotype with no error
ggplot(model_df.ne, aes(x = exp_env_factor, y = phen_corrected, group = gen_factor, fill = nat_env_factor,colour = nat_env_factor)) + geom_point() + geom_line() + theme_classic()

# Check: Mean Phenotype with no error
ggplot(mean_df.ne, aes(x = exp_env_factor, y = avg_phen_corrected, group = gen_factor, fill = nat_env_factor,colour = nat_env_factor)) + geom_point() + geom_line() + theme_classic()

###########################
##   RAW DATA ANALYSES   ##
###########################

# Anova model fit & GxE estimates
m1 <- mod.GxE(model_df) # Raw phenotype dataframe

# GxE Estimates
cov_matrix <- m1[[1]]
GxE_emm_original <- m1[[2]]
GxE_emm <- m1[[3]]
GxE_loop_output <- m1[[4]] # GxE output 
omega2 <- m1[[5]]
eta2 <- m1[[6]]
GxE_SSq <- m1[[7]] 
aov.df1 <- m1[[8]] # Model output

# Covariance Estimates
cov_est = cov(cov_matrix$G_means,cov_matrix$E_means)
cor_est = cor(cov_matrix$G_means,cov_matrix$E_means)
correction_raw = max(sd(cov_matrix$E_means),sd(cov_matrix$G_means))
cov_corrected = round(cov(cov_matrix$G_means, cov_matrix$E_means)/(correction_raw^2),2)

# Check: GxE Loop output
#hist(GxE_loop_output)

# Tracking: Anova Output
aov.df1$data_class <- rep("Raw_anova", nrow(aov.df1))

# Tracking: Covariance Matrix Output
cov_matrix$data.type = rep("raw",nrow(cov_matrix))

###################################
#####  Bootstrap -- Raw Data  #####
###################################

for(i in 1:n_boot){
  
  # Shuffle Data
  shuffle_dat <- bootstrap_raw(model_df) 
  
  # Anova model fit & GxE estimates
  m2 <- mod.GxE(shuffle_dat) # Insert shuffled raw phenotype dataframe
  
  # GxE Estimates
  cov_matrix_boot <- m2[[1]]
  GxE_emm_original_boot <- m2[[2]]
  GxE_emm_boot <- m2[[3]]
  GxE_loop_output_boot <- m2[[4]] # GxE output 
  omega2_boot <- m2[[5]]
  eta2_boot <- m2[[6]]
  GxE_SSq_boot <- m2[[7]] 

  # Covariance Estimates
  cov_est_boot = cov(cov_matrix_boot$G_means,cov_matrix_boot$E_means)
  cor_est_boot = cor(cov_matrix_boot$G_means,cov_matrix_boot$E_means)
  correction_raw_boot = max(sd(cov_matrix_boot$E_means),sd(cov_matrix_boot$G_means))
  cov_corrected_boot = round(cov(cov_matrix_boot$G_means, cov_matrix_boot$E_means)/(correction_raw_boot^2),2)
  
  # Bootstrap dataframe
  boot_dat_raw <- data.frame("covariance" = cov_est_boot,
                             "cor_est_boot" = cor_est_boot,
                             "cov_corrected_boot" = cov_corrected_boot,
                             "GxE_emm_original_boot" = GxE_emm_original_boot,
                             "GxE_emm_boot" = GxE_emm_boot,
                             "GxE_omega_boot" = omega2_boot,
                             "GxE_eta_boot" = eta2_boot,
                             "GxE_SSq_boot" = GxE_SSq_boot)
  boot_df_raw <- rbind(boot_df_raw,boot_dat_raw)
}

# Check: Histograms of Bootstrap
hist(boot_df_raw$GxE_emm_boot)

# Covariance Confidence Intervals 
cov_CI = quantile(boot_df_raw$covariance, probs=c(0.025, 0.975), type=1) 
cor_CI = quantile(boot_df_raw$cor_est_boot, probs=c(0.025, 0.975), type=1) 
cov_corrected_CI = quantile(boot_df_raw$cov_corrected_boot, probs=c(0.025, 0.975), type=1) 

# GxE Confidence Intervals
GxE_orig_CI = quantile(boot_df_raw$GxE_emm_original_boot, probs=c(0.025, 0.975), type=1) 
GxE_emm_CI = quantile(boot_df_raw$GxE_emm_boot, probs = c(0.025, 0.975), type=1)
GxE_omega_CI = quantile(boot_df_raw$GxE_omega_boot, probs=c(0.025, 0.975), type=1)
GxE_eta_CI = quantile(boot_df_raw$GxE_eta_boot, probs=c(0.025,0.975), type = 1)
GxE_SSq_CI = quantile(boot_df_raw$GxE_SSq_boot, probs = c(0.025,0.975), type = 1)

#######################################
#####   Permutation -- Raw Data   #####
#######################################

for(i in 1:n_boot){
  
  # Resample Data
  perm_dat <- permutation_raw(model_df)

  # Anova model fit & GxE estimates
  m3 <- mod.GxE(perm_dat) # Insert raw phenotype data frame
  
  # GxE Estimates
  cov_matrix_perm <- m3[[1]]
  GxE_emm_original_perm <- m3[[2]]
  GxE_emm_perm <- m3[[3]]
  GxE_loop_output_perm <- m3[[4]] # GxE output 
  omega2_perm <- m3[[5]]
  eta2_perm <- m3[[6]]
  GxE_SSq_perm <- m3[[7]] 
  
  # Covariance Estimates
  cov_est_perm = cov(cov_matrix_perm$G_means,cov_matrix_perm$E_means)
  cor_est_perm = cor(cov_matrix_perm$G_means,cov_matrix_perm$E_means)
  correction_raw_perm = max(sd(cov_matrix_perm$E_means),sd(cov_matrix_perm$G_means))
  cov_corrected_perm = round(cov(cov_matrix_perm$G_means, cov_matrix_perm$E_means)/(correction_raw_perm^2),2)
  
  # Permutation dataframe
  perm_dat_raw <- data.frame("covariance_perm" = cov_est_perm,
                             "cor_est_perm" = cor_est_perm,
                             "cov_corrected_perm" = cov_corrected_perm,
                             "GxE_emm_original_perm" = GxE_emm_original_perm,
                             "GxE_emm_perm" = GxE_emm_perm,
                             "GxE_omega_perm" = omega2_perm,
                             "GxE_eta_perm" = eta2_perm,
                             "GxE_SSq_perm" = GxE_SSq_perm)
  perm_df_raw <- rbind(perm_df_raw,perm_dat_raw)
}

# Check: Permutation histogram
hist(perm_df_raw$cov_corrected_perm)

# Covariance P-values
cov_original_pvalue <- pvalue_fun(cov_est,perm_df_raw$cov_est_perm,"twotail")
cor_pvalue <- pvalue_fun(cor_est,perm_df_raw$cor_est_perm,"twotail")
cov_corrected_pvalue <- pvalue_fun(cov_corrected,perm_df_raw$cov_corrected_perm,"twotail")

# GxE P-values
GxE_emm_orig_pvalue <- pvalue_fun(GxE_emm_original,perm_df_raw$GxE_emm_original_perm,"righttail")
GxE_emm_pvalue <- pvalue_fun(GxE_emm,perm_df_raw$GxE_emm_perm,"righttail")
GxE_omega_pvalue <- pvalue_fun(omega2,perm_df_raw$GxE_omega_perm,"righttail")
GxE_eta_pvalue <- pvalue_fun(eta2,perm_df_raw$GxE_eta_perm,"righttail")

############################
##   MEAN DATA ANALYSES   ##
############################

# GxE estimates
m4 <- mean.GxE(mean_df) # Insert means data frame

# GxE 
Cov_mean_matrix <- m4[[1]]
GxE_means <- m4[[2]]
GxE_means_loop_output <- m4[[3]]

# Covariance
cov_est_means = cov(Cov_mean_matrix$G_means,Cov_mean_matrix$E_means)
cor_est_means = cor(Cov_mean_matrix$G_means,Cov_mean_matrix$E_means)
means_correction = max(sd(Cov_mean_matrix$E_means),sd(Cov_mean_matrix$G_means))
cov_means_corrected = round(cov(Cov_mean_matrix$G_means, Cov_mean_matrix$E_means)/(means_correction^2),2)

# Check: GxE means Loop output
#hist(GxE_means_loop_output)

# Tracking: Covariance Matrix Output
Cov_mean_matrix$data.type = rep("means",nrow(Cov_mean_matrix))

###################################
##### BOOTSTRAP -- MEANS DATA #####
###################################

for(i in 1:n_boot){
  
  # Shuffle Data
  shuffle_means <- bootstrap_means(mean_df) # Insert means data
  
  # GxE :: Covariance Matrix
  m5 <- mean.GxE(shuffle_means) # Insert shuffled up means data frame
  
  # GxE Estimates
  Cov_mean_matrix_boot <- m5[[1]]
  GxE_means_boot <- m5[[2]]

  # Covariance Estimates
  cov_mean_boot = cov(Cov_mean_matrix_boot$G_means,Cov_mean_matrix_boot$E_means)
  cor_mean_boot = cor(Cov_mean_matrix_boot$G_means,Cov_mean_matrix_boot$E_means)
  correction_mean_boot = max(sd(Cov_mean_matrix_boot$E_means),sd(Cov_mean_matrix_boot$G_means))
  cov_corrected_mean_boot = round(cov(Cov_mean_matrix_boot$G_means, Cov_mean_matrix_boot$E_means)/(correction_mean_boot^2),2)
  
  # Bootstrap dataframe
  boot_dat_means <- data.frame("cov_means_boot" = cov_mean_boot,
                               "cor_mean_boot" = cor_mean_boot,
                               "cov_corrected_mean_boot" = cov_corrected_mean_boot,
                               "GxE_means_boot" = GxE_means_boot)
  boot_df_means <- rbind(boot_df_means,boot_dat_means)
}

# Covariance Confidence Intervals -- Means
cov_means_CI = quantile(boot_df_means$cov_means_boot, probs=c(0.025, 0.975), type=1) 
cor_means_CI = quantile(boot_df_means$cor_mean_boot, probs=c(0.025, 0.975), type=1) 
cov_corrected_means_CI = quantile(boot_df_means$cov_corrected_mean_boot, probs=c(0.025, 0.975), type=1) 

# GxE Confidence Intervals -- Means
GxE_means_CI = quantile(boot_df_means$GxE_means_boot, probs=c(0.025, 0.975), type=1) 

#######################################
#####  Permutation -- Means Data  #####
#######################################

for(i in 1:n_boot){
  
  # Resample Data
  perm_means <- permutation_means(mean_df)
  
  # GxE :: Covariance Matrix
  m6 <- mean.GxE(perm_means) # Insert resampled mean phenotype dataframe
  
  # GxE Estimates
  Cov_mean_matrix_perm <- m6[[1]]
  GxE_means_perm <- m6[[2]]
  GxE_means_output_perm <- m6[[3]]
  
  # Covariance Estimates
  cov_mean_perm = cov(Cov_mean_matrix_perm$G_means,Cov_mean_matrix_perm$E_means)
  cor_mean_perm = cor(Cov_mean_matrix_perm$G_means,Cov_mean_matrix_perm$E_means)
  correction_mean_perm = max(sd(Cov_mean_matrix_perm$E_means),sd(Cov_mean_matrix_perm$G_means))
  cov_corrected_mean_perm = round(cov(Cov_mean_matrix_perm$G_means, Cov_mean_matrix_perm$E_means)/(correction_mean_perm^2),2)
  
  # Check: GxE Histogram
  hist(GxE_means_output_perm)
  
  # Permutation dataframe -- Means
  perm_dat_means <- data.frame("cov_means_perm" = cov_mean_perm,
                               "cor_mean_perm" = cor_mean_perm,
                               "cov_corrected_mean_perm" = cov_corrected_mean_perm,
                               "GxE_means_perm" = GxE_means_perm)
  perm_df_means <- rbind(perm_df_means,perm_dat_means)
}

# Check: Histogram
hist(perm_df_means$GxE_means_perm)

# Covariance P-values
cov_original_mean_pvalue <- pvalue_fun(cov_est_means,perm_df_means$cov_means_perm,"twotail")
cor_mean_pvalue <- pvalue_fun(cor_est_means,perm_df_means$cor_mean_perm,"twotail")
cov_corrected_mean_pvalue <- pvalue_fun(cov_means_corrected,perm_df_means$cov_corrected_mean_perm,"twotail")

# GxE P-values
GxE_mean_pvalue <- pvalue_fun(GxE_means,perm_df_means$GxE_means_perm,"righttail")

####################################
##   NO ERROR RAW DATA ANALYSES   ##
####################################

# Anova model fit & GxE estimates
m7 <- mod.GxE(model_df.ne) # Insert raw phenotype data frame with no error

# GxE 
cov_matrix.ne <- m7[[1]]
GxE_emm_original.ne <- m7[[2]]
GxE_emm.ne <- m7[[3]]
GxE_loop_output.ne <- m7[[4]] # GxE output 
omega2.ne <- m7[[5]]
eta2.ne <- m7[[6]]
GxE_SSq.ne <- m7[[7]] 
aov.df1.ne <- m7[[8]] # Model output

# Covariance
cov_est.ne = cov(cov_matrix.ne$G_means,cov_matrix.ne$E_means)
cor_est.ne = cor(cov_matrix.ne$G_means,cov_matrix.ne$E_means)
correction_raw.ne = max(sd(cov_matrix.ne$E_means),sd(cov_matrix.ne$G_means))
cov_corrected.ne = round(cov(cov_matrix.ne$E_means, cov_matrix.ne$G_means)/(correction_raw.ne^2),2)

# Check: GxE Loop output
#hist(GxE_loop_output.ne)

# Tracking: Anova Output
aov.df1.ne$data_class <- rep("NoError_anova", nrow(aov.df1.ne))

# Tracking: Covariance Matrix
cov_matrix.ne$data.type <- rep("raw.ne",nrow(cov_matrix.ne))

######################################
##   NO ERROR MEANS DATA ANALYSES   ##
######################################

# GxE estimates
m8 <- mean.GxE(mean_df.ne) # Insert means data frame

# GxE 
Cov_mean_matrix.ne <- m8[[1]]
GxE_means.ne <- m8[[2]]
GxE_means_loop_output.ne <- m8[[3]]

# Covariance
cov_est_means.ne = cov(Cov_mean_matrix.ne$G_means,Cov_mean_matrix.ne$E_means)
cor_est_means.ne = cor(Cov_mean_matrix.ne$G_means,Cov_mean_matrix.ne$E_means)
means_correction.ne = max(sd(Cov_mean_matrix.ne$E_means),sd(Cov_mean_matrix.ne$G_means))
cov_means_corrected.ne = round(cov(Cov_mean_matrix.ne$G_means, Cov_mean_matrix.ne$E_means)/(means_correction.ne^2),2)

# Check: GxE means Loop output
#hist(GxE_means_loop_output.ne)

# Tracking: Covariance Matrix Output
Cov_mean_matrix.ne$data.type <- rep("mean.ne" , nrow(Cov_mean_matrix.ne))

################################
##   GENERATE FINAL OUTPUTS   ##
################################

# Time 
end.time <- Sys.time()
time.taken <- end.time - start.time

# Model Output
model_info = rbind(aov.df1,aov.df1.ne)

# Covariance Matrix Output
Cov_Matrix_Output = rbind(cov_matrix, Cov_mean_matrix, cov_matrix.ne, Cov_mean_matrix.ne)

# Bootstrap Output
boot_df <- cbind(boot_df_raw, boot_df_means)

# Permutation Output
perm_df <- cbind(perm_df_raw,perm_df_means)

# Generate Outputs
Parameters <- data.frame("row" = row, # Original Parameters
                         "replicate" = replicate,
                         "delta_env" = delta_env,
                         "delta_gen" = delta_gen,
                         "sample_size" = sample_size,
                         "n_pop" = n_pop,
                         "std_dev" = std_dev,
                         "interaction" = interaction,
                         "Sim_time" = time.taken)

Covariance <- data.frame("row" = row,
                         "true_cov" = cov_means_corrected.ne, #Corrected Covariance Estimates
                         "covariance" = cov_corrected, 
                         "covariance_lwrCI" = cov_corrected_CI[[1]],
                         "covariance_uprCI" = cov_corrected_CI[[2]],
                         "covariance_pvalue" = cov_corrected_pvalue,
                         
                         "true_cov_uncorrected" = cov_est.ne, #Covariance 
                         "cov_uncorrected" = cov_est,
                         "cov_uncorrected_lwrCI" = cov_CI[[1]],
                         "cov_uncorrected_uprCI" = cov_CI[[2]],
                         "cov_uncorrected_pvalue" = cov_original_pvalue,  
                         
                         "true_cor" = cor_est.ne, #Correlation 
                         "cor" = cor_est,
                         "cor_lwrCI" = cor_CI[[1]],
                         "cor_uprCI" = cor_CI[[2]],
                         "cor_pvalue" = cor_pvalue,
                         
                         "true_cov_means_correct" = cov_means_corrected.ne, # Corrected Covariance -- means
                         "cov_means_correct" = cov_means_corrected,
                         "cov_means_correct_lwrCI" = cov_corrected_means_CI[[1]],
                         "cov_means_correct_uprCI" = cov_corrected_means_CI[[2]],
                         "cov_means_correct_pvalue" = cov_corrected_mean_pvalue,
                         
                         "true_cor_means" = cor_est_means.ne, # Correlation -- means 
                         "cor_means" = cor_est_means,
                         "cor_means_lwrCI" = cor_means_CI[[1]],
                         "cor_means_uprCI" = cor_means_CI[[2]],
                         "cor_means_pvalue" = cor_mean_pvalue)

GxE <- data.frame("row" = row,
                  "true_GxE_emm" = round(GxE_emm.ne,2),# GxE Emmeans from loop
                  "GxE_emm" = round(GxE_emm,2), 
                  "GxE_emm_lwrCI" = round(GxE_emm_CI[[1]],2),
                  "GxE_emm_uprCI" = round(GxE_emm_CI[[2]],2),
                  "GxE_emm_pvalue" = round(GxE_emm_pvalue,2),
                  
                  "true_GxE_old" = round(GxE_emm_original.ne,2), # Emmeans GxE using original method (unused)
                  "GxE_emm_old" = round(GxE_emm_original,2),
                  "GxE_emm_old_lwrCI" = round(GxE_orig_CI[[1]],2),
                  "GxE_emm_old_uprCI" = round(GxE_orig_CI[[2]],2),
                  "GxE_emm_old_pvalue" = round(GxE_emm_orig_pvalue,2),
                  
                  "true_GxE_omega" = round(omega2.ne,2), # Omega^2 
                  "GxE_omega" = round(omega2,2), 
                  "GxE_omega_lwrCI" = round(GxE_omega_CI[[1]],2),
                  "GxE_omega_uprCI" = round(GxE_omega_CI[[2]],2),
                  "GxE_omega_pvalue" = round(GxE_omega_pvalue,2), 
                  
                  "true_GxE_eta" = round(eta2.ne,2), # Eta^2 
                  "GxE_eta" = round(eta2,2),
                  "GxE_eta_lwrCI" = round(GxE_eta_CI[[1]],2),
                  "GxE_eta_uprCI" = round(GxE_eta_CI[[2]],2),
                  "GxE_eta_pvalue" = round(GxE_eta_pvalue,2),
                  
                  "true_GxE_SSq" = round(GxE_SSq.ne,2), # Sums of Squares
                  "GxE_SSq" = round(GxE_SSq,2),
                  "GxE_SSq_lwrCI" = round(GxE_SSq_CI[[1]],2),
                  "GxE_SSq_uprCI" = round(GxE_SSq_CI[[2]],2),

                  "true_GxE_means" = round(GxE_means.ne,2), # GxE on means
                  "GxE_means" = round(GxE_means,2),
                  "GxE_means_lwrCI" = round(GxE_means_CI[[1]],2),
                  "GxE_means_uprCI" = round(GxE_means_CI[[2]],2),
                  "GxE_means_pvalue" = round(GxE_mean_pvalue,2)) 

# Write Files
#write.csv(GxE,paste0("/scratch/albecker/Power_analysis/power_output/GxE_",row,"_output.csv"))
#write.csv(Covariance,paste0("/scratch/albecker/Power_analysis/power_output/Covariance_",row,"_output.csv"))
#write.csv(Parameters,paste0("/scratch/albecker/Power_analysis/power_output/Parameters_",row,"_output.csv"))

#write.csv(phen_out,paste0("/scratch/albecker/Power_analysis/phenotype_output/Phenotype_data",row,"_output.csv"))
#write.csv(perm_df,paste0("/scratch/albecker/Power_analysis/permutation_output/Permutation_data",row,"_output.csv"))
#write.csv(boot_df,paste0("/scratch/albecker/Power_analysis/bootstrap_output/Bootstrap_data",row,"_output.csv"))
#write.csv(Cov_Matrix_Output,paste0("/scratch/albecker/Power_analysis/GEmeans_output/covmatrix_",row,"_output.csv"))
#write.csv(model_info,paste0("/scratch/albecker/Power_analysis/Anova_output/model_info_data",row,"_output.csv"))


