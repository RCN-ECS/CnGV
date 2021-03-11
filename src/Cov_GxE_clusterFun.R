
###############################################
#####   Data Input and Simulation         #####
###############################################

# Load packages 
library("emmeans")
emm_options(msg.interaction = FALSE)
library("lme4")
library("rlist")
library("dplyr")
library("ggplot2")
library("tibble")

# Load functions
source("Cov_GxE_functions.R")

start.time <- Sys.time()
args = commandArgs(trailingOnly = TRUE)

# Load Parameters
row = as.numeric(args[1])
n_pop = as.numeric(args[2])
sample_size = as.numeric(args[3])
std_dev = as.numeric(args[4])
n_env = as.numeric(args[5])
delta_env = as.numeric(args[6])
delta_gen = as.numeric(args[7])
interaction = as.numeric(args[8])
replicate = as.numeric(args[9])
env_scenario = as.numeric(args[10])
seed = as.numeric(args[11])

# Bootstraps
n_boot <- 999

# Clear Dataframes
Variance.partition <- var_perm <- PL_df <- phen_out <- 
Cov_Matrix_Output <- output_data <- output <- GEmeans_out <- 
model_out <- boot_df_raw <- boot_df_means <- perm_df_raw <-
perm_df_means <- data.frame()

# Establish seeds
set.seed(seed)
sim_seeds <- round(runif(10*n_boot)*1000000) # More than enough
seed1 = sim_seeds[1] # df.foundations
seed2 = sim_seeds[2] # df.foundations
seed3 = sim_seeds[3] # df.foundations2 when necessary
seed.set1 = sim_seeds[c(4:(4+n_boot))] # Bootstrap Means 
seed.set2 = sim_seeds[c((5+n_boot):(5+2*n_boot))] # Permutation means set 1
seed.set3 = sim_seeds[c((6+2*n_boot):(6+3*n_boot))] # Permutation means set 2
seed.set4 = sim_seeds[c((7+3*n_boot):(7+4*n_boot))] # Bootstrap raw
seed.set5 = sim_seeds[c((8+4*n_boot):(8+5*n_boot))] # Permutation raw
seed.set6 = sim_seeds[c((9+5*n_boot):(9+6*n_boot))] # Permutation raw

# Simulate data
if(env_scenario == 1){
  dfs <- df.foundations(delta_env, delta_gen, sample_size, n_env, std_dev, n_pop, interaction, seed1, seed2)
}else{
  dfs <- df.foundations2(delta_env, delta_gen, sample_size, n_env, std_dev, n_pop, interaction, seed1, seed2, seed3)
}

# Working dataframes
model_df <- dfs[[1]]     # Raw data
mean_df <- dfs[[2]]      # Means data
model_df.ne <- dfs[[3]]  # Raw data with no error
mean_df.ne <- dfs[[4]]   # Means data with no error
varpar.df <- dfs[[5]]    # Variance Partitioning

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

# Check: Raw Phenotype (All 4 plots should look similar)
# ggplot(model_df, aes(x = exp_env_factor, y = phen_corrected, group = gen_factor, fill = nat_env_factor,colour = nat_env_factor)) + geom_point() + geom_smooth(se = F) + theme_classic()

# Check: Mean Phenotype 
# ggplot(mean_df, aes(x = exp_env_factor, y = avg_phen_corrected, group = gen_factor, fill = nat_env_factor,colour = nat_env_factor)) + geom_point() + geom_line() + theme_classic()

# Check: Raw Phenotype with no error
# ggplot(model_df.ne, aes(x = exp_env_factor, y = phen_corrected, group = gen_factor, fill = nat_env_factor,colour = nat_env_factor)) + geom_point() + geom_line() + theme_classic()

# Check: Mean Phenotype with no error
# ggplot(mean_df.ne, aes(x = exp_env_factor, y = avg_phen_corrected, group = gen_factor, fill = nat_env_factor,colour = nat_env_factor)) + geom_point() + geom_line() + theme_classic()

# Check: Variance Partitioning
# ggplot(varpar.df, aes(x = exp_env_factor, y = phen_corrected, group = gen_factor, fill = nat_env_factor,colour = nat_env_factor)) + geom_smooth(se = FALSE) + geom_point() + theme_classic()

###############################################
#####           Raw Data Analysis         #####
###############################################

# Anova model fit & GxE estimates
m1 <- mod.GxE(model_df, is.perm = FALSE, seed = NA) # Raw phenotype dataframe

# GxE Estimates
cov_matrix <- m1[[1]]
GxE_emm_original <- m1[[2]]
GxE_emm <- m1[[3]]
GxE_loop_output <- m1[[4]] # All GxE estimates from loop 
omega2 <- m1[[5]]
eta2 <- m1[[6]]
GxE_SSq <- m1[[7]] 
aov.df1 <- m1[[8]] # Anova SSq output
aov_coefs <- m1[[11]]

# Covariance Estimates
cov_est = cov(cov_matrix$G_means,cov_matrix$E_means)
cor_est = cor(cov_matrix$G_means,cov_matrix$E_means)
cov_corrected = round(cov.function(cov_matrix, phen_df = model_df, is.sample = TRUE),9)

# Stamps and Hadfield Info
if (n_pop == 2) {
  
  delta_E <- m1[[9]]
  delta_H <- m1[[10]]
  Pl <- (delta_E/delta_H)
  
  G1Emean <- cov_matrix[1,4]
  G2Emean <- cov_matrix[2,4]
  E1Gmean <- cov_matrix[1,1]
  E2Gmean <- cov_matrix[2,1]
  
  # PL Output (for Stamps + Hadfield response)
  PL_df <- data.frame(row =  row, PL = Pl, delta_E = delta_E, delta_H = delta_H, 
                      G1Emean = G1Emean, G2Emean = G2Emean, E1Gmean = E1Gmean, E2Gmean = E2Gmean, 
                      GxE_emm = GxE_emm, covariance = cov_corrected)
  
  write.csv(PL_df,paste0("/scratch/albecker/Power_analysis/PL_output/PL_",row,".csv"))
}

# Variance Partitioning
Variance.partition1 = var.partition(varpar.df)
Variance.partition1 = data.frame(row, Variance.partition1)

# Tracking: Anova Output
aov.df1$data_class <- rep("Raw_anova", nrow(aov.df1))
aov_coefs1 = data.frame(aov_coefs, "row" = row)

# Tracking: Covariance Matrix Output
cov_matrix$data.type = rep("raw",nrow(cov_matrix))


###############################################
#####        Bootstrap - Raw Data         #####
###############################################

# Output
boot_df_raw = data.frame()

for(i in 1:n_boot){
  
  # Sampling seed
  bootseed <- seed.set4[i]
  
  # Shuffle Data
  shuffle_dat <- bootstrap_raw(model_df) 
  
  # Anova model fit & GxE estimates
  m2 <- mod.GxE(shuffle_dat, is.perm = FALSE, seed = NA) # Insert shuffled raw phenotype dataframe
  
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
  cov_corrected_boot = round(cov.function(cov_matrix_boot, phen_df = shuffle_dat, is.sample = TRUE),3)

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
# ggplot(boot_df_raw, aes(x = cov_corrected_boot), alpha = 0.5)+
#   geom_histogram()+ geom_vline(aes(xintercept = GxE_emm))+theme_classic()+ ggtitle("Bootstrap: Raw Data")

# Covariance Confidence Intervals 
cov_CI = quantile(boot_df_raw$covariance, probs=c(0.025, 0.975), type=8) 
cor_CI = quantile(boot_df_raw$cor_est_boot, probs=c(0.025, 0.975), type=8) 
cov_corrected_CI = quantile(boot_df_raw$cov_corrected_boot, probs=c(0.025, 0.975), type=8) 

# GxE Confidence Intervals
GxE_orig_CI = quantile(boot_df_raw$GxE_emm_original_boot, probs=c(0.025, 0.975), type=8) 
GxE_emm_CI = quantile(boot_df_raw$GxE_emm_boot, probs = c(0.025, 0.975), type = 8)
GxE_omega_CI = quantile(boot_df_raw$GxE_omega_boot, probs=c(0.025, 0.975), type = 8)
GxE_eta_CI = quantile(boot_df_raw$GxE_eta_boot, probs=c(0.025,0.975), type = 8)
GxE_SSq_CI = quantile(boot_df_raw$GxE_SSq_boot, probs = c(0.025,0.975), type = 8)


###############################################
#####      Permutation -- Raw Data        #####
###############################################

# Output
perm_df_raw = data.frame()

for(i in 1:n_boot){
  
  # Sampling Seed
  perm.seed = seed.set5[i]
  perm.seed2 = seed.set6[i]
  
  # Resample Data
  perm_dat <- permutation_raw(model_df, perm.seed)
  var_dat <- permutation_varpar(varpar.df)
  
  # Anova model fit & GxE estimates
  m3 <- mod.GxE(perm_dat, is.perm = TRUE, seed = perm.seed2) # Insert shuffled permutation raw data frame
  
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
  cov_corrected_perm = round(cov.function(cov_matrix_perm, phen_df = perm_dat, is.sample = TRUE),3)
  
  # Variance Partitioning
  v3 <- var.partition(var_dat)
  v3$permid = i
  var_perm = rbind(v3, var_perm) 
  
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
# ggplot(perm_df_raw, aes(x = GxE_emm_perm), alpha = 0.5)+geom_histogram()+xlim(0,(GxE_emm+0.1))+geom_vline(aes(xintercept = GxE_emm))+theme_classic()+ ggtitle("Null Histogram: Raw")
 
# Covariance P-values
cov_original_pvalue <- pvalue_fun(cov_est,perm_df_raw$covariance_perm,"twotail",n_boot)
cor_pvalue <- pvalue_fun(cor_est,perm_df_raw$cor_est_perm,"twotail",n_boot)
cov_corrected_pvalue <- pvalue_fun(cov_corrected,perm_df_raw$cov_corrected_perm,"twotail",n_boot)

# GxE P-values
GxE_emm_orig_pvalue <- pvalue_fun(GxE_emm_original,perm_df_raw$GxE_emm_original_perm,"righttail",n_boot)
GxE_emm_pvalue <- pvalue_fun(GxE_emm,perm_df_raw$GxE_emm_perm,"righttail",n_boot)
GxE_omega_pvalue <- pvalue_fun(omega2,perm_df_raw$GxE_omega_perm,"righttail",n_boot)
GxE_eta_pvalue <- pvalue_fun(eta2,perm_df_raw$GxE_eta_perm,"righttail",n_boot)

# Variance Partitioning P-values
Vg = Ve = Vgxe = Vcov = Verror = NULL

for(i in 1:length(unique(var_perm$permid))){
  id = unique(var_perm$permid)[i]
  sub = filter(var_perm, permid == id)
  Vg. = sub$omega2[1]
  Ve. = sub$omega2[2]
  Vgxe. = sub$omega2[3]
  Vcov. = sub$omega2[4]
  Verror. = sub$omega2[5]
  
  Vg = rbind(Vg, Vg.)
  Ve = rbind(Ve, Ve.)
  Vgxe = rbind(Vgxe, Vgxe.)
  Vcov = rbind(Vcov, Vcov.)
  Verror = rbind(Verror, Verror.)
}

Vg.pval <- pvalue_fun(Variance.partition1$omega2[1],Vg,"righttail",n_boot)
Ve.pval <- pvalue_fun(Variance.partition1$omega2[2],Ve[,1],"righttail",n_boot)
Vgxe.pval <- pvalue_fun(Variance.partition1$omega2[3],Vgxe[,1],"righttail",n_boot)
Vcov.pval <- pvalue_fun(Variance.partition1$omega2[4],Vcov[,1],"righttail",n_boot)

Varpar.p.values = data.frame("Variance_Component" = c("V_g", "V_e", "V_gxe", "V_cov", "V_error"),
                             "Pvalue" = c(Vg.pval,Ve.pval,Vgxe.pval,Vcov.pval,NA))


###############################################
#####         Mean Data Analysis          ##### 
###############################################

# GxE estimates
m4 <- mean.GxE(mean_df,is.perm = FALSE, seed = NA) # Insert means data frame (seed not necessary if is.perm is False)

# Covariance and GxE for Means 
Cov_mean_matrix <- m4[[1]]

# Fix ordering issue for N_pop > 10 (How R orders factors)
if(length(Cov_mean_matrix$gen_factor) > 9){
  Cov_mean_matrix[,2] = cov_matrix$gen_factor
  Cov_mean_matrix[,3] = cov_matrix$exp_env_factor
  colnames(Cov_mean_matrix)[4] <- "JustKidding"
  Cov_mean_matrix$E_means = cov_matrix$E_means[match(cov_matrix$exp_env_factor, Cov_mean_matrix$exp_env_factor)]
}

GxE_means <- m4[[2]]
GxE_means_loop_output <- m4[[3]]

# Covariance
cov_est_means = cov(Cov_mean_matrix$G_means,Cov_mean_matrix$E_means)
cor_est_means = cor(Cov_mean_matrix$G_means,Cov_mean_matrix$E_means)
cov_means_corrected = round(cov.function_means(Cov_mean_matrix, phen_df = mean_df, is.sample = TRUE),3)

# Tracking: Covariance Matrix Output
Cov_mean_matrix$data.type = rep("means",nrow(Cov_mean_matrix))


###############################################
#####       Bootstrap -- Means Data       #####
###############################################

for(i in 1:n_boot){
  
  # Shuffle Data
  boot.seed = seed.set1
  shuffle_means <- bootstrap_means(mean_df) # Insert means data, Need n_boot seeds
  
  # GxE :: Covariance Matrix
  m5 <- mean.GxE(shuffle_means, is.perm = FALSE, seed = NA) # Insert shuffled up means data frame
  
  # GxE Estimates
  Cov_mean_matrix_boot <- m5[[1]]
  
  if(length(Cov_mean_matrix_boot$gen_factor) > 9){
    Cov_mean_matrix_boot[,2] = cov_matrix_boot$gen_factor
    Cov_mean_matrix_boot[,3] = cov_matrix_boot$exp_env_factor
    colnames(Cov_mean_matrix_boot)[4] <- "JustKidding"
    Cov_mean_matrix_boot$E_means = cov_matrix_boot$E_means[match(cov_matrix_boot$exp_env_factor, Cov_mean_matrix_boot$exp_env_factor)]
  }
  
  GxE_means_boot <- m5[[2]]
  
  # Covariance Estimates
  cov_mean_boot = cov(Cov_mean_matrix_boot$G_means,Cov_mean_matrix_boot$E_means)
  cor_mean_boot = cor(Cov_mean_matrix_boot$G_means,Cov_mean_matrix_boot$E_means)
  cov_corrected_mean_boot = round(cov.function_means(Cov_mean_matrix_boot,phen_df =shuffle_means, is.sample = TRUE),3)
  
  # Bootstrap dataframe
  boot_dat_means <- data.frame("cov_means_boot" = cov_mean_boot,
                               "cor_mean_boot" = cor_mean_boot,
                               "cov_corrected_mean_boot" = cov_corrected_mean_boot,
                               "GxE_means_boot" = GxE_means_boot)
  boot_df_means <- rbind(boot_df_means,boot_dat_means)
}

# Histograms 
# ggplot(boot_df_means, aes(x = GxE_means_boot), alpha = 0.5)+geom_histogram()+xlim(0,(GxE_means+0.1))+geom_vline(aes(xintercept = GxE_means))+theme_classic()+ ggtitle("Bootstrap: Means Data")

# Covariance Confidence Intervals -- Means
cov_means_CI = quantile(boot_df_means$cov_means_boot, probs=c(0.025, 0.975), na.rm = TRUE, type=8) 
cor_means_CI = quantile(boot_df_means$cor_mean_boot, probs=c(0.025, 0.975), na.rm = TRUE, type=8) 
cov_corrected_means_CI = quantile(boot_df_means$cov_corrected_mean_boot, probs=c(0.025, 0.975), na.rm = TRUE, type=8) 

# GxE Confidence Intervals -- Means
GxE_means_CI = quantile(boot_df_means$GxE_means_boot, probs=c(0.025, 0.975), na.rm = TRUE, type=8) 

###############################################
#####      Permutation -- Means Data      #####
###############################################


for(i in 1:n_boot){
  
  # Set seeds for perm_means and mean.GxE
  perm.seeds1 = seed.set2
  perm.seeds2 = seed.set3
  
  # Resample Data
  perm_means <- permutation_means(mean_df,perm.seeds1[i])
  
  # GxE :: Covariance Matrix
  m6 <- mean.GxE(perm_means, is.perm = TRUE,perm.seeds2[i]) # Insert resampled mean phenotype dataframe
  
  # GxE Estimates
  Cov_mean_matrix_perm <- m6[[1]]
  if(length(Cov_mean_matrix_perm$gen_factor) > 9){
    Cov_mean_matrix_perm[,2] = cov_matrix_perm$gen_factor
    Cov_mean_matrix_perm[,3] = cov_matrix_perm$exp_env_factor
    colnames(Cov_mean_matrix_perm)[4] <- "JustKidding"
    Cov_mean_matrix_perm$E_means = cov_matrix_perm$E_means[match(cov_matrix_perm$exp_env_factor, Cov_mean_matrix_perm$exp_env_factor)]
  }
  GxE_means_perm <- m6[[2]]
  GxE_means_output_perm <- m6[[3]]
  
  # Covariance Estimates
  cov_mean_perm = cov(Cov_mean_matrix_perm$G_means,Cov_mean_matrix_perm$E_means)
  cor_mean_perm = cor(Cov_mean_matrix_perm$G_means,Cov_mean_matrix_perm$E_means)
  cov_corrected_mean_perm = round(cov.function_means(Cov_mean_matrix_perm, phen_df = perm_means, is.sample = TRUE),3)
  
  # Permutation dataframe -- Means
  perm_dat_means <- data.frame("cov_means_perm" = cov_mean_perm,
                               "cor_mean_perm" = cor_mean_perm,
                               "cov_corrected_mean_perm" = cov_corrected_mean_perm,
                               "GxE_means_perm" = GxE_means_perm)
  perm_df_means <- rbind(perm_df_means,perm_dat_means)
}

# Check: Histogram
# ggplot(perm_df_means, aes(x = GxE_means_perm), alpha = 0.5)+geom_histogram()+ geom_vline(aes(xintercept = GxE_means))+theme_classic()+ggtitle("Null Histogram: Means")

# Covariance P-values
cov_original_mean_pvalue <- pvalue_fun(cov_est_means,perm_df_means$cov_means_perm,"twotail", n_boot)
cor_mean_pvalue <- pvalue_fun(cor_est_means,perm_df_means$cor_mean_perm,"twotail", n_boot)
cov_corrected_mean_pvalue <- pvalue_fun(cov_means_corrected,perm_df_means$cov_corrected_mean_perm,"twotail", n_boot)

# GxE P-values
GxE_mean_pvalue <- pvalue_fun(GxE_means,perm_df_means$GxE_means_perm,"righttail",n_boot)


###############################################
#####      Population -- Raw Data         #####
###############################################

# Anova model fit & GxE estimates
m7 <- mod.GxE(model_df.ne, is.perm = FALSE, seed = NA) # Insert raw phenotype data frame with no error

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
cov_corrected.ne = round(cov.function(cov_matrix.ne, phen_df = model_df.ne, is.sample = FALSE),3)

# Check: GxE Loop output
# hist(GxE_loop_output.ne)

# Tracking: Anova Output
aov.df1.ne$data_class <- rep("NoError_anova", nrow(aov.df1.ne))

# Tracking: Covariance Matrix
cov_matrix.ne$data.type <- rep("raw.ne",nrow(cov_matrix.ne))

###############################################
#####       Population -- Means Data      #####
###############################################

# GxE estimates
m8 <- mean.GxE(mean_df.ne, is.perm = FALSE, seed = NA) # Insert means data frame

# GxE 
Cov_mean_matrix.ne <- m8[[1]]
if(length(Cov_mean_matrix.ne$gen_factor) > 9){
  Cov_mean_matrix.ne[,2] = cov_matrix.ne$gen_factor
  Cov_mean_matrix.ne[,3] = cov_matrix.ne$exp_env_factor
  colnames(Cov_mean_matrix.ne)[4] <- "JustKidding"
  Cov_mean_matrix.ne$E_means = cov_matrix.ne$E_means[match(cov_matrix.ne$exp_env_factor, Cov_mean_matrix.ne$exp_env_factor)]
}
GxE_means.ne <- m8[[2]]
GxE_means_loop_output.ne <- m8[[3]]

# Covariance
cov_est_means.ne = cov(Cov_mean_matrix.ne$G_means,Cov_mean_matrix.ne$E_means)
cor_est_means.ne = cor(Cov_mean_matrix.ne$G_means,Cov_mean_matrix.ne$E_means)
cov_means_corrected.ne = round(cov.function_means(Cov_mean_matrix.ne, phen_df = mean_df.ne, is.sample = FALSE),3)

# Tracking: Covariance Matrix Output
Cov_mean_matrix.ne$data.type <- rep("mean.ne" , nrow(Cov_mean_matrix.ne))


###############################################
#####     Bookkeeping -- Final Output     #####
###############################################

# Time 
end.time <- Sys.time()
time.taken <- end.time - start.time

# Variance Partitioning
Variance.partition1$row = row
Variance.partition = full_join(Variance.partition1,Varpar.p.values, by = "Variance_Component")

# Phenotype Data
phen_out$row <- row

# Anova SSQ
model_info = rbind(aov.df1,aov.df1.ne)
model_info$row = row

# Covariance Matrix Output
if(n_pop > 10){
  Cov_Matrix_Output = rbind(cov_matrix, Cov_mean_matrix[,-4], cov_matrix.ne, Cov_mean_matrix.ne[,-4])
}else{
  Cov_Matrix_Output = rbind(cov_matrix, Cov_mean_matrix, cov_matrix.ne, Cov_mean_matrix.ne)}
Cov_Matrix_Output$row = row

# Bootstrap Output
boot_df <- cbind(boot_df_raw, boot_df_means)
boot_df$row = row

# Permutation Output
perm_df <- cbind(perm_df_raw, perm_df_means)
perm_df$row = row

# Generate Outputs
output_data <- data.frame("row" = row, # Original Parameters
                          "replicate" = replicate,
                          "env_scenario" = env_scenario,
                          "delta_env" = delta_env,
                          "delta_gen" = delta_gen,
                          "sample_size" = sample_size,
                          "total_samples" = sample_size*n_pop*n_env,
                          "n_env" = n_env,
                          "n_pop" = n_pop,
                          "std_dev" = std_dev,
                          "interaction" = interaction,
                          "Sim_time" = time.taken,
                          
                          #"pop_cov" = cov_corrected.ne.pop, # Corrected Covariance Estimates
                          "true_cov" = cov_corrected.ne,
                          #"cov_R" = cov_corrected_check,
                          "covariance" = cov_corrected, 
                          "covariance_lwrCI" = cov_corrected_CI[[1]],
                          "covariance_uprCI" = cov_corrected_CI[[2]],
                          "covariance_pvalue" = cov_corrected_pvalue,
                          
                          #"true_cov_uncorrected" = cov_est.ne, #Covariance 
                          #"cov_uncorrected" = cov_est,
                          #"cov_uncorrected_lwrCI" = cov_CI[[1]],
                          #"cov_uncorrected_uprCI" = cov_CI[[2]],
                          #"cov_uncorrected_pvalue" = cov_original_pvalue,  
                          
                          "true_cov_means" = cov_means_corrected.ne, # Corrected Covariance -- means
                          #"cov_means_R" = cov_means_corrected_check,
                          "cov_means" = cov_means_corrected,
                          "cov_means_lwrCI" = cov_corrected_means_CI[[1]],
                          "cov_means_uprCI" = cov_corrected_means_CI[[2]],
                          "cov_means_pvalue" = cov_corrected_mean_pvalue,
                          
                          #"true_cor" = round(cor_est.ne,3), # Correlation 
                          #"cor" = round(cor_est,3),
                          #"cor_lwrCI" = round(cor_CI[[1]],3),
                          #"cor_uprCI" = round(cor_CI[[2]],3),
                          #"cor_pvalue" = cor_pvalue,
                          
                          #"true_cor_means" = cor_est_means.ne, # Correlation -- means 
                          #"cor_means" = cor_est_means,
                          #"cor_means_lwrCI" = cor_means_CI[[1]],
                          #"cor_means_uprCI" = cor_means_CI[[2]],
                          #"cor_means_pvalue" = cor_mean_pvalue,
                          
                          "GxE_Anova" = round(aov.df1[3,6],3),
                          "true_GxE_emm" = round(GxE_emm.ne,3),# GxE Emmeans from loop
                          "GxE_emm" = round(GxE_emm,3), 
                          "GxE_emm_lwrCI" = round(GxE_emm_CI[[1]],3),
                          "GxE_emm_uprCI" = round(GxE_emm_CI[[2]],3),
                          "GxE_emm_pvalue" = round(GxE_emm_pvalue,3),
                          
                          #"true_GxE_old" = round(GxE_emm_original.ne,2), # Emmeans GxE using original method (unused)
                          #"GxE_emm_old" = round(GxE_emm_original,2),
                          #"GxE_emm_old_lwrCI" = round(GxE_orig_CI[[1]],2),
                          #"GxE_emm_old_uprCI" = round(GxE_orig_CI[[2]],2),
                          #"GxE_emm_old_pvalue" = round(GxE_emm_orig_pvalue,2),
                          
                          "true_GxE_omega" = round(omega2.ne,3), # Omega^2 
                          "GxE_omega" = round(omega2,3), 
                          "GxE_omega_lwrCI" = round(GxE_omega_CI[[1]],3),
                          "GxE_omega_uprCI" = round(GxE_omega_CI[[2]],3),
                          "GxE_omega_pvalue" = round(GxE_omega_pvalue,3), 
                          
                          #"true_GxE_eta" = round(eta2.ne,2), # Eta^2 
                          #"GxE_eta" = round(eta2,2),
                          #"GxE_eta_lwrCI" = round(GxE_eta_CI[[1]],2),
                          #"GxE_eta_uprCI" = round(GxE_eta_CI[[2]],2),
                          #"GxE_eta_pvalue" = round(GxE_eta_pvalue,2),
                          
                          #"true_GxE_SSq" = round(GxE_SSq.ne,2), # Sums of Squares
                          #"GxE_SSq" = round(GxE_SSq,2),
                          #"GxE_SSq_lwrCI" = round(GxE_SSq_CI[[1]],2),
                          #"GxE_SSq_uprCI" = round(GxE_SSq_CI[[2]],2),
                          
                          "true_GxE_means" = round(GxE_means.ne,3), # GxE on means
                          "GxE_means" = round(GxE_means,3),
                          "GxE_means_lwrCI" = round(GxE_means_CI[[1]],3),
                          "GxE_means_uprCI" = round(GxE_means_CI[[2]],3),
                          "GxE_means_pvalue" = round(GxE_mean_pvalue,3)) 

# Write Files
write.csv(output_data,paste0("/scratch/albecker/Power_analysis/power_output/Results_",row,".csv"))
write.csv(phen_out,paste0("/scratch/albecker/Power_analysis/phenotype_output/Phenotype_data",row,".csv"))
write.csv(perm_df,paste0("/scratch/albecker/Power_analysis/permutation_output/Permutation_data",row,".csv"))
write.csv(boot_df,paste0("/scratch/albecker/Power_analysis/bootstrap_output/Bootstrap_data",row,".csv"))
write.csv(Cov_Matrix_Output,paste0("/scratch/albecker/Power_analysis/GEmeans_output/covmatrix_",row,".csv"))
write.csv(model_info,paste0("/scratch/albecker/Power_analysis/Anova_output/model_SSQ_data",row,".csv"))
write.csv(aov_coefs1,paste0("/scratch/albecker/Power_analysis/Anova_output/anova_coef_data",row,".csv"))
write.csv(Variance.partition, paste0("/scratch/albecker/Power_analysis/power_output/variance_",row,".csv"))


