
# Load packages 
library("emmeans")
emm_options(msg.interaction = FALSE)
library("lme4")
library("rlist")
library("dplyr")
library("ggplot2")

# Load functions
#source("Cov_GxE_functions.R")

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
n_boot <- 100

# Output dataframes
output <- data.frame()
phen_out <- data.frame()
GEmeans_out <- data.frame()
boot_df <- data.frame()
perm_df <- data.frame()

# Simulate data
dfs <- df.foundations(row, replicate, delta_env, delta_gen, sample_size, n_env, std_dev, n_pop, interaction)

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
ggplot(model_df.ne, aes(x = exp_env_factor, y = phen_corrected_ne, group = gen_factor, fill = nat_env_factor,colour = nat_env_factor)) + geom_point() + geom_smooth() + theme_classic()

# Check: Mean Phenotype with no error
ggplot(mean_df.ne, aes(x = exp_env_factor, y = avg_phen_corrected_ne, group = gen_factor, fill = nat_env_factor,colour = nat_env_factor)) + geom_point() + geom_smooth() + theme_classic()

###########################
##   RAW DATA ANALYSES   ##
###########################

# Anova model fit & GxE estimates
m1 <- mod.GxE(model_df) # Insert raw phenotype data frame

# GxE 
cov_matrix <- m1[[1]]
GxE_emm_original <- m1[[2]]
GxE_emm <- m1[[3]]
GxE_loop_output <- m1[[4]] # GxE output 
omega2 <- m1[[5]]
eta2 <- m1[[6]]
GxE_SSq <- m1[[7]] 
aov.df1 <- m1[[8]] # Model output

# Covariance
cov_est = cov(cov_matrix$G_means,cov_matrix$E_means)
cor_est = cor(cov_matrix$G_means,cov_matrix$E_means)
correction_raw = max(sd(cov_matrix$E_means),sd(cov_matrix$G_means))
cov_corrected = round(cov(cov_matrix$G_means, cov_matrix$E_means)/(correction_raw^2),2)

# Check: GxE Loop output
hist(GxE_loop_output)

# Tracking: Anova Output
aov.df1$data_class <- rep("Raw_anova", nrow(aov.df1))


############################
##   MEAN DATA ANALYSES   ##
############################

# GxE estimates
m2 <- mean.GxE(mean_df) # Insert means data frame

# GxE 
Cov_mean_matrix <- m2[[1]]
GxE_means <- m2[[2]]
GxE_means_loop_output <- m2[[3]]

# Covariance
cov_est_means = cov(Cov_mean_matrix$G_means,Cov_mean_matrix$E_means)
cor_est_means = cor(Cov_mean_matrix$G_means,Cov_mean_matrix$E_means)
means_correction = max(sd(Cov_mean_matrix$E_means),sd(Cov_mean_matrix$G_means))
cov_means_corrected = round(cov(Cov_mean_matrix$G_means, Cov_mean_matrix$E_means)/(means_correction^2),2)

# Check: GxE means Loop output
hist(GxE_means_loop_output)

####################################
##   NO ERROR RAW DATA ANALYSES   ##
####################################

# Anova model fit & GxE estimates
m3 <- mod.GxE(model_df.ne) # Insert raw phenotype data frame with no error

# GxE 
cov_matrix.ne <- m3[[1]]
GxE_emm_original.ne <- m3[[2]]
GxE_emm.ne <- m3[[3]]
GxE_loop_output.ne <- m3[[4]] # GxE output 
omega2.ne <- m3[[5]]
eta2.ne <- m3[[6]]
GxE_SSq.ne <- m3[[7]] 
aov.df1.ne <- m3[[8]] # Model output

# Covariance
cov_est.ne = cov(cov_matrix.ne$G_means,cov_matrix.ne$E_means)
cor_est.ne = cor(cov_matrix.ne$G_means,cov_matrix.ne$E_means)
correction_raw.ne = max(sd(cov_matrix.ne$E_means),sd(cov_matrix.ne$G_means))
cov_corrected.ne = round(cov(cov_matrix.ne$E_means, cov_matrix.ne$G_means)/(correction_raw.ne^2),2)

# Check: GxE Loop output
hist(GxE_loop_output.ne)

# Tracking: Anova Output
aov.df1.ne$data_class <- rep("NoError_anova", nrow(aov.df1.ne))



