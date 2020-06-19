
# Load packages
library("emmeans")
emm_options(msg.interaction = FALSE)
library("lme4")
library("rlist")
library("dplyr")
library("ggplot2")

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
  
# Dataframe foundations
n_environments <- n_env 
gen <- rep(1:n_pop, each = sample_size*n_environments)
env <- rep(1:n_environments, times = n_pop, each = sample_size) 
  
# Random Noise
noise <- rnorm(sample_size * n_pop * n_environments, 0, sd = std_dev) 
  
# Interaction Term
int <- rep(rnorm(n_pop * n_environments, 0, sd = interaction), each = sample_size) # interaction term - one for each GE level
  
# Create model dataframe 
model_df <- data.frame(gen, env, noise, int)
model_df$gen_factor <- factor(paste("G", model_df$gen, sep = "_"))
model_df$exp_env_factor <- factor(paste("E", model_df$env, sep = "_"))
  
# Native environments
model_df$nat_env_factor = rep(unique(model_df$exp_env_factor), each = (sample_size*n_pop))
  
# Phenotype data using regression equation
phen = delta_env * model_df$env + delta_gen * model_df$gen  + model_df$noise + model_df$int 
model_df$phen = phen
  
# Generate means dataframe
mean_df = model_df %>%
  group_by(gen_factor, exp_env_factor) %>%
  summarize("avg_phen" = mean(phen),
            "deviation" = sd(phen))
mean_df$se <- mean_df$deviation/(sqrt(sample_size))
mean_df$nat_env_factor <- model_df$nat_env_factor[match(mean_df$gen_factor,model_df$gen_factor)]
  
# Standardize data
model_df$phen_corrected <- (phen - mean(phen))/sd(phen)
mean_df$avg_phen_corrected <- (mean_df$avg_phen - mean(mean_df$avg_phen))/sd(mean_df$avg_phen) 
  
# Create Phenotype output
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
  
# Raw Phenotype plot
#ggplot(phen_out, aes(x = exp_env_factor, y = phen, group = gen_factor, fill = nat_env_factor,colour = nat_env_factor)) + geom_point() + geom_smooth() + theme_classic()

# Mean Phenotype plot
#ggplot(mean_df, aes(x = exp_env_factor, y = avg_phen_corrected, group = gen_factor, fill = nat_env_factor,colour = nat_env_factor)) + geom_point() + geom_smooth() + theme_classic()

# Anova
test_temp <- lm(phen_corrected ~ exp_env_factor * gen_factor, data = model_df)

# Save Model Output
mod_df <- as.data.frame(summary(aov(test_temp))[[1]])
mod_df <- rownames_to_column(mod_df) 
colnames(mod_df)[1] <- "Fixed_effect"
mod_df$data_class <- rep("Raw_anova", nrow(mod_df))
GxE_anova_pvalue <- mod_df[3,6]

# Estimated Marginal Means
emm_E = as.data.frame(emmeans(test_temp,"exp_env_factor"))
emm_G = as.data.frame(emmeans(test_temp, "gen_factor"))
emm_GxE = as.data.frame(emmeans(test_temp, ~ exp_env_factor*gen_factor))

# Gmeans
G_matrix = data.frame()
for(h in 1:length(unique(emm_GxE$gen_factor))){
gtemp <- filter(emm_GxE, gen_factor == unique(emm_GxE$gen_factor)[h])
gmean <- mean(gtemp[,3])
tempdat = data.frame("G_means" = gmean,
                     "gen_factor" = unique(emm_GxE$gen_factor)[h])
G_matrix = rbind(G_matrix,tempdat)
}

# Emeans
E_matrix = data.frame()
for(j in 1:length(unique(emm_GxE$exp_env_factor))){
etemp <- filter(emm_GxE, exp_env_factor == unique(emm_GxE$exp_env_factor)[j])
emean <- mean(etemp[,3])
tempdat. = data.frame("E_means" = emean,
                      "exp_env_factor" = unique(emm_GxE$exp_env_factor)[j])
E_matrix = rbind(E_matrix,tempdat.)
}

# Match Genotypes to Native Environment
Cov_matrix = G_matrix
Cov_matrix$exp_env_factor <- model_df$nat_env_factor[match(G_matrix$gen_factor,model_df$gen_factor)]
Cov_matrix$E_means <- E_matrix$E_means[match(Cov_matrix$exp_env_factor,E_matrix$exp_env_factor)]

# Means of Means
E_means <- tapply(mean_df$avg_phen_corrected, mean_df$exp_env_factor, mean)
G_means <- tapply(mean_df$avg_phen_corrected, mean_df$gen_factor, mean)
Gmean_mat <- data.frame("G_means" = G_means, "gen_factor" = unique(mean_df$gen_factor))
Emean_mat <- data.frame("E_means" = E_means, "exp_env_factor" = unique(mean_df$exp_env_factor))

# Match means to native
Cov_mean_matrix = Gmean_mat
Cov_mean_matrix$exp_env_factor <- mean_df$nat_env_factor[match(Cov_mean_matrix$gen_factor,mean_df$gen_factor)]
Cov_mean_matrix$E_means <- Emean_mat$E_means[match(Cov_mean_matrix$exp_env_factor,Emean_mat$exp_env_factor)]

# Covariances
cov_est = cov(Cov_matrix$G_means,Cov_matrix$E_means)
cor_est = cor(Cov_matrix$G_means,Cov_matrix$E_means)
correction_ratio = max(sd(Cov_matrix$E_means),sd(Cov_matrix$G_means))
cov_corrected = round(cov(Cov_matrix$E_means, Cov_matrix$G_means)/(correction_ratio^2),2)

cov_est_means = cov(Cov_mean_matrix$G_means,Cov_mean_matrix$E_means)
cor_est_means = cor(Cov_mean_matrix$G_means,Cov_mean_matrix$E_means)
means_correction = max(sd(Cov_mean_matrix$E_means),sd(Cov_mean_matrix$G_means))
cov_means_corrected = round(cov(Cov_mean_matrix$E_means, Cov_mean_matrix$G_means)/(means_correction^2),2)

# Magnitude of GxE -- EMMs
GxE_emm_original<- abs(emm_GxE$emmean[emm_GxE$gen_factor == "G_1" & emm_GxE$exp_env_factor == "E_1"] - # GxE (Phenotype of ith genotype in jth environment)
                       emm_G$emmean[emm_G$gen_factor == "G_1"] - # phenotype of ith Genotype
                       emm_E$emmean[emm_E$exp_env_factor == "E_1"] + # phenotype of jth Environment
                       mean(emm_GxE$emmean)) # Overall mean phenotype

# Magnitude of GxE -- Loop
allGE <- c()
for (i in 1:nlevels(emm_GxE$gen_factor)){
  for (j in 1:nlevels(emm_GxE$exp_env_factor)){
    G_levels <- levels(emm_GxE$gen_factor)
    E_levels <- levels(emm_GxE$exp_env_factor)
    loopGxE <- abs(emm_GxE$emmean[emm_GxE$gen_factor == G_levels[i] & emm_GxE$exp_env_factor == E_levels[j]] - # GxE (Phenotype of ith genotype in jth environment)
                  emm_G$emmean[emm_G$gen_factor == G_levels[i]] - # phenotype of ith Genotype
                  emm_E$emmean[emm_E$exp_env_factor == E_levels[j]] + # phenotype of jth Environment
                  mean(emm_GxE$emmean)) # Overall mean
  allGE <- c(allGE, loopGxE)
}
}
#hist(allGE)
GxE_emm_loop = mean(allGE)

# Magnitude of GxE -- Omega^2
w2_GxE = (summary(aov(test_temp))[[1]][3,2] - #(SS_effect -
         (summary(aov(test_temp))[[1]][3,1]*summary(aov(test_temp))[[1]][4,3])) / #(Df_effect * MS_error))/
         (sum(summary(aov(test_temp))[[1]][,2]) + # (SS_total+
         (summary(aov(test_temp))[[1]][4,3])) # MS_error)

# Magnitude of GxE -- Eta2 
eta_GxE = summary(aov(test_temp))[[1]][3,2]/sum(summary(aov(test_temp))[[1]][,2])

# Residual Variation -- proportion of sums of squares without error
ResVar = summary(aov(test_temp))[[1]][3,2]/sum(summary(aov(test_temp))[[1]][c(1:3),2])

# Magnitude of GxE -- Loop -- Means
allGEmeans <- c()
for (i in 1:nlevels(mean_df$gen_factor)){
for (j in 1:nlevels(mean_df$exp_env_factor)){
  G_levels <- levels(mean_df$gen_factor)
  E_levels <- levels(mean_df$exp_env_factor)
  GxE_emm_mean1 <- abs(mean_df$avg_phen_corrected[mean_df$gen_factor == G_levels[i] & mean_df$exp_env_factor == E_levels[j]] - 
                      mean(mean_df$avg_phen_corrected[mean_df$gen_factor == G_levels[i]]) - # phenotype of ith Genotype
                      mean(mean_df$avg_phen_corrected[mean_df$exp_env_factor == E_levels[j]]) + # phenotype of jth Environment
                      mean(mean_df$avg_phen_corrected)) # Overall mean
  allGEmeans <- c(allGEmeans, GxE_emm_mean1)
}
}
#hist(allGEmeans)
GxE_means = mean(allGEmeans)

#################################
## No Error Covariance and GxE ##
#################################

# Phenotype data with no standard deviation
no_err_phen = delta_env * model_df$env + delta_gen * model_df$gen + model_df$int
model_df$no_err_phen = no_err_phen

# Generate means dataframe
mean_df_ne = model_df %>%
group_by(gen_factor, exp_env_factor) %>%
summarize(avg_phen_ne = mean(no_err_phen),
          deviation_ne = sd(no_err_phen))
mean_df_ne$se_ne = mean_df_ne$deviation_ne/(sqrt(sample_size)) # Should be zero 

# Add native environments
mean_df_ne$nat_env_factor = model_df$nat_env_factor[match(mean_df_ne$gen_factor,model_df$gen_factor)]

# Standardize -- no error data
model_df$phen_corrected_ne = (no_err_phen - mean(no_err_phen))/sd(no_err_phen)
mean_df_ne$avg_phen_corrected_ne = (mean_df_ne$avg_phen_ne - mean(mean_df_ne$avg_phen_ne))/sd(mean_df_ne$avg_phen_ne) 

# Anova -- no error
test_temp_noerror <- lm(phen_corrected_ne ~ exp_env_factor * gen_factor, data = model_df)

# Save Model Output
mod_df_ne = as.data.frame(summary(aov(test_temp_noerror))[[1]])
mod_df_ne = tibble::rownames_to_column(mod_df_ne) 
colnames(mod_df_ne)[1] <- "Fixed_effect"
mod_df_ne$data_class = rep("Raw_Anova_noError", nrow(mod_df_ne))

# Estimated Marginal Means -- no error
emm_E_noerror = as.data.frame(emmeans(test_temp_noerror,"exp_env_factor"))
emm_G_noerror = as.data.frame(emmeans(test_temp_noerror, "gen_factor"))
emm_GxE_noerror = as.data.frame(emmeans(test_temp_noerror, ~ exp_env_factor*gen_factor))

# Gmeans -- no error
G_matrix_noerror = data.frame()
for(y in 1:length(unique(emm_GxE_noerror$gen_factor))){
gtemp_ne <- filter(emm_GxE_noerror, gen_factor == unique(emm_GxE_noerror$gen_factor)[y])
gmean_ne <- mean(gtemp_ne[,3])
tempdat_ne = data.frame("G_means" = gmean_ne,
                        "gen_factor" = unique(emm_GxE_noerror$gen_factor)[y])
G_matrix_noerror = rbind(G_matrix_noerror,tempdat_ne)
}

# Emeans -- no error
E_matrix_noerror = data.frame()
for(f in 1:length(unique(emm_GxE_noerror$exp_env_factor))){
etemp_ne <- filter(emm_GxE_noerror, exp_env_factor == unique(emm_GxE_noerror$exp_env_factor)[f])
emean_ne <- mean(etemp_ne[,3])
tempdat._ne = data.frame("E_means" = emean_ne,
                        "exp_env_factor" = unique(emm_GxE_noerror$exp_env_factor)[f])
E_matrix_noerror = rbind(E_matrix_noerror,tempdat._ne)
}

# Match Genotypes to Native Environment
Cov_matrix_ne = G_matrix_noerror
Cov_matrix_ne$exp_env_factor <- model_df$nat_env_factor[match(G_matrix_noerror$gen_factor,model_df$gen_factor)]
Cov_matrix_ne$E_means <- E_matrix_noerror$E_means[match(Cov_matrix_ne$exp_env_factor,E_matrix_noerror$exp_env_factor)]

# Means on Means -- stacks on stacks of means
E_means_ne <- tapply(mean_df_ne$avg_phen_corrected_ne, mean_df_ne$exp_env_factor, mean)
G_means_ne <- tapply(mean_df_ne$avg_phen_corrected_ne, mean_df_ne$gen_factor, mean)
Gmean_mat_ne <- data.frame("G_means" = G_means_ne, "gen_factor" = unique(mean_df$gen_factor))
Emean_mat_ne <- data.frame("E_means" = E_means_ne, "exp_env_factor" = unique(mean_df$exp_env_factor))

#Match means to native
Cov_mean_matrix_ne = Gmean_mat_ne
Cov_mean_matrix_ne$exp_env_factor <- mean_df$nat_env_factor[match(Cov_mean_matrix_ne$gen_factor,mean_df$gen_factor)]
Cov_mean_matrix_ne$E_means <- Emean_mat_ne$E_means[match(Cov_mean_matrix_ne$exp_env_factor,Emean_mat_ne$exp_env_factor)]

# True Covariances 
true_cov = cov(Cov_matrix_ne$G_means,Cov_matrix_ne$E_means)
true_cor = cor(Cov_matrix_ne$G_means,Cov_matrix_ne$E_means)
correction_ratio_noerror = max(sd(Cov_matrix_ne$E_means),sd(Cov_matrix_ne$G_means))
true_cov_corrected = round(cov(Cov_matrix_ne$G_means, Cov_matrix_ne$E_means)/(correction_ratio_noerror^2),2)

true_cov_means = cov(Cov_mean_matrix_ne$G_means, Cov_mean_matrix_ne$E_means)
true_cor_means = cor(Cov_mean_matrix_ne$G_means, Cov_mean_matrix_ne$E_means)
means_correction_ne = max(sd(Cov_mean_matrix_ne$E_means),sd(Cov_mean_matrix_ne$G_means))
true_cov_means_corrected = round(cov(Cov_mean_matrix_ne$G_means, Cov_mean_matrix_ne$E_means)/(means_correction_ne^2),2)

# Magnitude of GxE -- EMM -- no error
true_GxE_original <- abs(mean(emm_GxE_noerror$emmean) - # Overall mean
                         mean(emm_GxE_noerror$emmean[emm_GxE_noerror$gen_factor=="G_1"])- # G
                         mean(emm_GxE_noerror$emmean[emm_GxE_noerror$exp_env_factor=="E_1"])+ # E
                         emm_GxE_noerror$emmean[emm_GxE_noerror$gen_factor == "G_1" & emm_GxE_noerror$exp_env_factor == "E_1"]) # GxE

# Magnitude of True GxE -- Loop
allGE_ne <- c()
for (i in 1:nlevels(model_df$gen_factor)){
for (j in 1:nlevels(model_df$exp_env_factor)){
  G_levels <- levels(model_df$gen_factor)
  E_levels <- levels(model_df$exp_env_factor)
  loopGxE_ne  <- abs(mean(emm_GxE_noerror$emmean) - # Overall mean
                     emm_G_noerror$emmean[emm_G_noerror$gen_factor == G_levels[i]] - # average phenotype of ith Genotype
                     emm_E_noerror$emmean[emm_E_noerror$exp_env_factor == E_levels[j]] + # average phenotype of jth Environment
                     emm_GxE_noerror$emmean[emm_GxE_noerror$gen_factor == G_levels[i] & emm_GxE_noerror$exp_env_factor == E_levels[j]]) # Phenotype of interaction
  allGE_ne <- c(allGE_ne, loopGxE_ne)
}
}
#hist(allGE_ne)
GxE_emm_ne_loop = mean(allGE_ne)

# Magnitude of GxE -- Omega^2 -- no error
true_w2_GxE = (summary(aov(test_temp_noerror))[[1]][3,2] - #SS_effect -
            (summary(aov(test_temp_noerror))[[1]][3,1]*summary(aov(test_temp_noerror))[[1]][4,3])) / #Df_effect * MS_error
            (sum(summary(aov(test_temp_noerror))[[1]][,2]) + # SS_total+
            (summary(aov(test_temp_noerror))[[1]][4,3])) # MS_error

# Magnitude of GxE -- Eta2 -- No Error
eta_GxE_ne = summary(aov(test_temp_noerror))[[1]][3,2]/sum(summary(aov(test_temp_noerror))[[1]][,2])

# Residual Variation
ResVar_ne = summary(aov(test_temp_noerror))[[1]][3,2]/(summary(aov(test_temp_noerror))[[1]][1,2]+summary(aov(test_temp_noerror))[[1]][2,2]+summary(aov(test_temp_noerror))[[1]][3,2])

# Magnitude of GxE -- Looped Means -- No Error
allGEmeans_ne <- c()
for (i in 1:nlevels(Cov_mean_matrix_ne$gen_factor)){
for (j in 1:nlevels(Cov_mean_matrix_ne$exp_env_factor)){
  G_levels <- levels(Cov_mean_matrix_ne$gen_factor)
  E_levels <- levels(Cov_mean_matrix_ne$exp_env_factor)
  GxE_emm_mean_ne <- abs(mean(mean_df_ne$avg_phen_corrected_ne) - # Overall mean
                        Cov_mean_matrix_ne$G_means[Cov_mean_matrix_ne$gen_factor == G_levels[i]] - # phenotype of ith Genotype
                        unique(Cov_mean_matrix_ne$E_means[Cov_mean_matrix_ne$exp_env_factor == E_levels[j]]) + # phenotype of jth Environment
                        mean_df_ne$avg_phen_corrected_ne[mean_df_ne$gen_factor == G_levels[i] & mean_df_ne$exp_env_factor == E_levels[j]]) # phen of ith and jth
  
  allGEmeans_ne <- c(allGEmeans_ne, GxE_emm_mean_ne)
}
}
#hist(allGE)
GxE_means_ne = mean(allGEmeans_ne)

## Generate Means output
means_out <- data.frame("row" = rep(unique(row),length(G_means)), 
                        "replicate" = rep(unique(replicate),length(G_means)), 
                        "delta_env" = rep(unique(delta_env),length(G_means)), 
                        "delta_gen" = rep(unique(delta_gen),length(G_means)), 
                        "sample_size" = rep(unique(sample_size),length(G_means)), 
                        "n_pop" = rep(unique(n_pop),length(G_means)), 
                        "n_env" = rep(unique(n_env),length(G_means)),
                        "std_dev" = rep(unique(std_dev),length(G_means)), 
                        "interaction" = rep(unique(interaction),length(G_means)),
                        "gen_factor" = Cov_matrix$gen_factor,
                        "Anova.Gmeans_EMM" = Cov_matrix$G_means,
                        "Anova.Gmeans_EMM_ne" = Cov_matrix_ne$G_means,
                        "Means.Gmeans" = Cov_mean_matrix$G_means,
                        "Means.Gmeans_ne" = Cov_mean_matrix_ne$G_means,
                        "exp_env_factor" = Cov_matrix$exp_env_factor, 
                        "Anova.Emeans_EMM" = Cov_matrix$E_means,
                        "Anova.Emeans_EMM_ne" = Cov_matrix_ne$E_means,
                        "Means.Emeans" = Cov_mean_matrix$E_means,
                        "Means.Emeans_ne" = Cov_mean_matrix_ne$E_means)
                        
GEmeans_out <- rbind(GEmeans_out,means_out)

###############
## Bootstrap ##
###############

# Resampling Loop
for(a in 1:n_boot){

new_phen <- NULL
shuffle_dat <- data.frame()

# Resample data within each genotype and environment
for (l in 1:nlevels(model_df$gen_factor)){
  for (j in 1:nlevels(model_df$exp_env_factor)){
    
    cond <- model_df %>%
              filter(gen_factor == unique(model_df$gen_factor)[l]) %>%
              filter(exp_env_factor == unique(model_df$exp_env_factor)[j])
    
    # Shuffle data 
    new_phen <- sample(cond$phen_corrected, size=nrow(cond), replace=TRUE)
    
    # Output    
    shuffle_dat_temp <- data.frame("gen_factor" = cond$gen_factor,
                                   "exp_env_factor" = cond$exp_env_factor,
                                   "phen_corrected" = new_phen)
    shuffle_dat <- rbind(shuffle_dat, shuffle_dat_temp)
  }
}

# Resample means dataframe
new_means <- data.frame()

for (u in 1:nlevels(mean_df$gen_factor)){
  for (r in 1:nlevels(mean_df$exp_env_factor)){
    
    # Retain levels
    cond <- mean_df %>%
      filter(gen_factor == unique(mean_df$gen_factor)[u]) %>%
      filter(exp_env_factor == unique(mean_df$exp_env_factor)[r])
    
    # Create new means data
    new_phen <- rnorm(nrow(cond), mean = cond$avg_phen, sd = cond$se)
    
    # Output
    new_mean_temp <- data.frame("gen_factor" = cond$gen_factor,
                                "exp_env_factor" = cond$exp_env_factor,
                                "mean_phen" = new_phen)
    new_means <- rbind(new_means, new_mean_temp)
  }
}

# Standardize resampled means
new_means$new_mean_corrected = (new_means$mean_phen - mean(new_means$mean_phen))/sd(new_means$mean_phen) 

# Anova -- bootstrap
test_boot <- lm(phen_corrected ~ exp_env_factor * gen_factor, data = shuffle_dat)

# Estimated Marginal Means -- bootstrap
emm_E_boot = as.data.frame(emmeans(test_boot,"exp_env_factor"))
emm_G_boot = as.data.frame(emmeans(test_boot, "gen_factor"))
emm_GxE_boot = as.data.frame(emmeans(test_boot, ~ exp_env_factor*gen_factor))

# Gmeans -- Bootstrap
G_matrix_boot = data.frame()
for(p in 1:length(unique(emm_GxE_boot$gen_factor))){
  g_boot <- filter(emm_GxE_boot, gen_factor == unique(emm_GxE_boot$gen_factor)[p])
  g_mean_boot <- mean(g_boot[,3])
  tempdat = data.frame("G_means" = g_mean_boot,
                       "gen_factor" = unique(g_boot$gen_factor))
  G_matrix_boot = rbind(G_matrix_boot,tempdat)
}

# Emeans -- Bootstrap
E_matrix_boot = data.frame()
for(q in 1:length(unique(emm_GxE_boot$exp_env_factor))){
  e_boot <- filter(emm_GxE_boot, exp_env_factor == unique(emm_GxE_boot$exp_env_factor)[q])
  e_mean_boot <- mean(e_boot[,3])
  tempdat. = data.frame("E_means" = e_mean_boot,
                        "exp_env_factor" = unique(e_boot$exp_env_factor))
  E_matrix_boot = rbind(E_matrix_boot,tempdat.)
}

# Match Genotypes to Native Environment
Cov_matrix_boot = G_matrix_boot
Cov_matrix_boot$exp_env_factor <- model_df$nat_env_factor[match(G_matrix_boot$gen_factor,model_df$gen_factor)]
Cov_matrix_boot$E_means <- E_matrix_boot$E_means[match(Cov_matrix_boot$exp_env_factor,E_matrix_boot$exp_env_factor)]

# Means on Means -- bootstrap
E_means_shuffle <- tapply(new_means$new_mean_corrected, new_means$exp_env_factor, mean)
G_means_shuffle <- tapply(new_means$new_mean_corrected, new_means$gen_factor, mean)
Gmean_mat_boot <- data.frame("G_means" = G_means_shuffle, "gen_factor" = unique(new_means$gen_factor))
Emean_mat_boot <- data.frame("E_means" = E_means_shuffle, "exp_env_factor" = unique(new_means$exp_env_factor))

# Match means to Native
Cov_mean_matrix_boot = Gmean_mat_boot
Cov_mean_matrix_boot$exp_env_factor <- mean_df$nat_env_factor[match(Cov_mean_matrix_boot$gen_factor,mean_df$gen_factor)]
Cov_mean_matrix_boot$E_means <- Emean_mat_boot$E_means[match(Cov_mean_matrix_boot$exp_env_factor,Emean_mat_boot$exp_env_factor)]

# Covariances -- Bootstrap 
cov_est_boot = cov(Cov_matrix_boot$G_means,Cov_matrix_boot$E_means)
cor_est_boot = cor(Cov_matrix_boot$G_means,Cov_matrix_boot$E_means)
correction_ratio_boot = max(sd(Cov_matrix_boot$E_means),sd(Cov_matrix_boot$G_means))
cov_corrected_boot = round(cov(Cov_matrix_boot$G_means, Cov_matrix_boot$E_means)/(correction_ratio_boot^2),2)

cov_means_boot = cov(Cov_mean_matrix_boot$G_means,Cov_mean_matrix_boot$E_means) 
cor_means_boot = cor(Cov_mean_matrix_boot$G_means,Cov_mean_matrix_boot$E_means)
means_correction_shuffle = max(sd(Cov_mean_matrix_boot$E_means),sd(Cov_mean_matrix_boot$G_means)) 
cov_means_correct_boot = round(cov(Cov_mean_matrix_boot$G_means,Cov_mean_matrix_boot$E_means)/(means_correction_shuffle^2),2)

# Magnitude of GxE -- EMMs -- Bootstrap
GxE_emm_boot <- abs(mean(emm_GxE_boot$emmean) - # Overall mean
                    emm_G_boot$emmean[emm_G_boot$gen_factor=="G_1"]- # G
                    emm_E_boot$emmean[emm_E_boot$exp_env_factor=="E_1"]+ # E
                    emm_GxE_boot$emmean[emm_GxE_boot$gen_factor=="G_1"&emm_GxE_boot$exp_env_factor=="E_1"]) # GxE

# Magnitude of GxE -- Loop
allGE_boot <- c()
for (i in 1:nlevels(shuffle_dat$gen_factor)){
  for (j in 1:nlevels(shuffle_dat$exp_env_factor)){
    G_levels <- levels(shuffle_dat$gen_factor)
    E_levels <- levels(shuffle_dat$exp_env_factor)
    GxE_emm_boot1 <- abs(mean(emm_GxE_boot$emmean) - 
                        emm_G_boot$emmean[emm_G_boot$gen_factor == G_levels[i]] - # phenotype of ith Genotype
                        emm_E_boot$emmean[emm_E_boot$exp_env_factor == E_levels[j]] + # phenotype of jth Environment
                        emm_GxE_boot$emmean[emm_GxE_boot$gen_factor == G_levels[i] & emm_GxE_boot$exp_env_factor == E_levels[j]])
    allGE_boot <- c(allGE_boot, GxE_emm_boot1)
  }
}
#hist(allGE_boot)
loopGxE_boot = mean(allGE_boot)

# Magnitude of GxE -- Omega2 -- Bootstrap
boot_w2_GxE = (summary(aov(test_boot))[[1]][3,2] - #SS_effect -
              (summary(aov(test_boot))[[1]][3,1]*summary(aov(test_boot))[[1]][4,3])) / #Df_effect * MS_error
              (sum(summary(aov(test_boot))[[1]][,2]) + # SS_total+
              (summary(aov(test_boot))[[1]][4,3])) # MS_error

# Magnitude of GxE -- Eta2 -- Bootstrap
boot_eta_GxE = summary(aov(test_boot))[[1]][3,2]/sum(summary(aov(test_boot))[[1]][,2]) 

# Residual Variation
ResVar_boot = summary(aov(test_boot))[[1]][3,2]/(summary(aov(test_boot))[[1]][1,2]+summary(aov(test_boot))[[1]][2,2]+summary(aov(test_boot))[[1]][3,2])
     
# Magnitude of GxE -- Means -- Bootstrap
allGE_means_boot <- c()
for (i in 1:nlevels(new_means$gen_factor)){
  for (j in 1:nlevels(new_means$exp_env_factor)){
    G_levels <- levels(new_means$gen_factor)
    E_levels <- levels(new_means$exp_env_factor)
    GxE_emm_boot_means <- abs(mean(new_means$new_mean_corrected) - # Overall mean
                              mean(new_means$new_mean_corrected[new_means$gen_factor == G_levels[i]]) - # phenotype of ith Genotype
                              mean(new_means$new_mean_corrected[new_means$exp_env_factor == E_levels[j]]) + # phenotype of jth Environment
                              new_means$new_mean_corrected[new_means$gen_factor == G_levels[i] & new_means$exp_env_factor == E_levels[j]]) # phen of ith and jth
    allGE_means_boot <- c(allGE_means_boot, GxE_emm_boot_means)
  }
}
#hist(allGE_means_boot)
GxE_means_boot = mean(allGE_means_boot)

# Bootstrap dataframe
boot_dat <- data.frame("covariance" = cov_est_boot,
                       "cor_est_boot" = cor_est_boot,
                       "cov_corrected_boot" = cov_corrected_boot,
                       "cov_means_boot" = cov_means_boot,
                       "cor_means_boot" = cor_means_boot,
                       "cov_means_correct_boot" = cov_means_correct_boot,
                       "GxE_mag" = GxE_emm_boot,
                       "GxE_mag_loop" = loopGxE_boot,
                       "GxE_mag_omega" = boot_w2_GxE,
                       "GxE_eta_boot" = boot_eta_GxE,
                       "GxE_resvar_boot" = ResVar_boot,
                       "GxE_means_boot" = GxE_means_boot)
boot_df <- rbind(boot_df,boot_dat)
}

# Covariance Confidence Intervals
cov_CI = quantile(boot_df$covariance, probs=c(0.025, 0.975), type=1) 
cor_CI = quantile(boot_df$cor_est_boot, probs=c(0.025, 0.975), type=1) 
cov_corrected_CI = quantile(boot_df$cov_corrected_boot, probs=c(0.025, 0.975), type=1) 
cov_means_CI = quantile(boot_df$cov_means_boot, probs=c(0.025, 0.975), type=1) 
cor_means_CI = quantile(boot_df$cor_means_boot, probs=c(0.025, 0.975), type=1) 
cov_means_correct_CI = quantile(boot_df$cov_means_correct_boot, probs=c(0.025, 0.975), type=1) 

# GxE Confidence Intervals
GxE_CI = quantile(boot_df$GxE_mag, probs=c(0.025, 0.975), type=1) 
GxE_loop_CI = quantile(boot_df$GxE_mag_loop, probs = c(0.025, 0.975), type=1)
GxE_omega_CI = quantile(boot_df$GxE_mag_omega, probs=c(0.025, 0.975), type=1)
GxE_eta_CI = quantile(boot_df$GxE_eta_boot, probs=c(0.025,0.975), type = 1)
GxE_resvar_CI = quantile(boot_df$GxE_resvar_boot, probs = c(0.025,0.975), type = 1)
GxE_means_CI = quantile(boot_df$GxE_means_boot, probs=c(0.025, 0.975), type=1)

#################
## Permutation ##
#################

for(b in 1:n_boot){  

# Shuffle raw data
null_temp <- sample(model_df$phen_corrected, size=nrow(model_df), replace=FALSE)
null_temp_ne <- sample(model_df$phen_corrected_ne, size=nrow(model_df), replace=FALSE)

perm_dat <- data.frame("gen_factor" = model_df$gen_factor,
                       "exp_env_factor" = model_df$exp_env_factor,
                       "phen" = null_temp,
                       "phen_ne" = null_temp_ne)

# Shuffle means data
#null_means. <- rnorm(nrow(mean_df), mean = mean_df$avg_phen, sd = mean_df$se)
#null_means <- sample(null_means., size=length(null_means.), replace=FALSE)
null_means <- sample(mean_df$avg_phen_corrected, size=nrow(mean_df), replace=FALSE)      
null_means_ne <- sample(mean_df_ne$avg_phen_corrected_ne, size=nrow(mean_df_ne), replace=FALSE)

perm_means <- data.frame("gen_factor" = mean_df$gen_factor,
                         "exp_env_factor" = mean_df$exp_env_factor,
                         "phen_data" = null_means,
                         "phen_data_ne" = null_means_ne)

# Re-Standardize 
perm_dat$phen_corrected = perm_dat$phen#(perm_dat$phen - mean(perm_dat$phen))/sd(perm_dat$phen)
perm_dat$phen_corrected_ne = perm_dat$phen_ne#(perm_dat$phen_ne - mean(perm_dat$phen_ne))/sd(perm_dat$phen_ne)
perm_means$avg_phen_corrected = perm_means$phen_data#(perm_means$phen_data - mean(perm_means$phen_data))/sd(perm_means$phen_data)
perm_means$avg_phen_corrected_ne = perm_means$phen_data_ne#(perm_means$phen_data_ne - mean(perm_means$phen_data_ne))/sd(perm_means$phen_data_ne)

# Anovas
test_perm <- lm(phen_corrected ~ exp_env_factor * gen_factor, data = perm_dat)
test_perm_ne <- lm(phen_corrected_ne ~ exp_env_factor * gen_factor, data = perm_dat)

# Estimated Marginal Means -- Permutation
emm_E_perm = as.data.frame(emmeans(test_perm,"exp_env_factor"))
emm_G_perm = as.data.frame(emmeans(test_perm, "gen_factor"))
emm_GxE_perm = as.data.frame(emmeans(test_perm, ~ exp_env_factor*gen_factor))

# Estimated Marginal Means -- Permutation -- No Error
emm_E_perm_ne = as.data.frame(emmeans(test_perm_ne,"exp_env_factor"))
emm_G_perm_ne = as.data.frame(emmeans(test_perm_ne, "gen_factor"))
emm_GxE_perm_ne = as.data.frame(emmeans(test_perm_ne, ~ exp_env_factor*gen_factor))

# Gmeans - Permutation - EMMs
G_matrix_perm = data.frame()
for(r in 1:length(unique(emm_GxE_perm$gen_factor))){
  g_perm <- filter(emm_GxE_perm, gen_factor == unique(emm_GxE_perm$gen_factor)[r])
  g_mean_perm <- sum(g_perm[,3])/length(unique(emm_GxE_perm$gen_factor))
  tempdat = data.frame("G_means" = g_mean_perm,
                       "gen_factor" = unique(g_perm$gen_factor))
  G_matrix_perm = rbind(G_matrix_perm,tempdat)
}

# Emeans - Permutation - EMMs
E_matrix_perm = data.frame()
for(s in 1:length(unique(emm_GxE_perm$exp_env_factor))){
  e_perm <- filter(emm_GxE_perm, exp_env_factor == unique(emm_GxE_perm$exp_env_factor)[s])
  e_mean_perm <- sum(e_perm[,3])/length(unique(emm_GxE_perm$exp_env_factor))
  tempdat. = data.frame("E_means" = e_mean_perm,
                        "exp_env_factor" = unique(e_perm$exp_env_factor))
  E_matrix_perm = rbind(E_matrix_perm,tempdat.)
}

# Match Genotypes to Native Environment
Cov_matrix_perm = G_matrix_perm
Cov_matrix_perm$exp_env_factor <- model_df$nat_env_factor[match(G_matrix_perm$gen_factor,model_df$gen_factor)]
Cov_matrix_perm$E_means <- E_matrix_perm$E_means[match(Cov_matrix_boot$exp_env_factor,E_matrix_perm$exp_env_factor)]

# Gmeans - Permutation - No error
G_matrix_perm_ne = data.frame()
for(e in 1:length(unique(emm_GxE_perm_ne$gen_factor))){
  g_perm_ne <- filter(emm_GxE_perm_ne, gen_factor == unique(emm_GxE_perm_ne$gen_factor)[e])
  g_mean_perm_ne <- sum(g_perm_ne[,3])/length(unique(emm_GxE_perm_ne$gen_factor))
  tempdat_ne = data.frame("G_means" = g_mean_perm_ne,
                          "gen_factor" = unique(g_perm_ne$gen_factor))
  G_matrix_perm_ne = rbind(G_matrix_perm_ne,tempdat_ne)
}

# Emeans - Permutation - No Error
E_matrix_perm_ne = data.frame()
for(o in 1:length(unique(emm_GxE_perm_ne$exp_env_factor))){
  e_perm_ne <- filter(emm_GxE_perm_ne, exp_env_factor == unique(emm_GxE_perm_ne$exp_env_factor)[o])
  e_mean_perm_ne <- sum(e_perm_ne[,3])/length(unique(emm_GxE_perm_ne$exp_env_factor))
  tempdat._ne = data.frame("E_means" = e_mean_perm_ne,
                           "exp_env_factor" = unique(e_perm_ne$exp_env_factor))
  E_matrix_perm_ne = rbind(E_matrix_perm_ne,tempdat._ne)
}

# Match Genotypes to Native Environment - No Error
Cov_matrix_perm_ne = G_matrix_perm_ne
Cov_matrix_perm_ne$exp_env_factor <- model_df$nat_env_factor[match(G_matrix_perm_ne$gen_factor,model_df$gen_factor)]
Cov_matrix_perm_ne$E_means <- E_matrix_perm_ne$E_means[match(Cov_matrix_perm_ne$exp_env_factor,E_matrix_perm_ne$exp_env_factor)]

# Means on Means -- Permutation
E_means_perm <- tapply(perm_means$avg_phen_corrected, perm_means$exp_env_factor, mean)
G_means_perm <- tapply(perm_means$avg_phen_corrected, perm_means$gen_factor, mean)
Gmean_mat_perm <- data.frame("G_means" = G_means_perm, "gen_factor" = unique(perm_means$gen_factor))
Emean_mat_perm <- data.frame("E_means" = E_means_perm, "exp_env_factor" = unique(perm_means$exp_env_factor))

# Match means to Native
Cov_mean_matrix_perm = Gmean_mat_perm
Cov_mean_matrix_perm$exp_env_factor <- mean_df$nat_env_factor[match(Cov_mean_matrix_perm$gen_factor,mean_df$gen_factor)]
Cov_mean_matrix_perm$E_means <- Emean_mat_perm$E_means[match(Cov_mean_matrix_perm$exp_env_factor,Emean_mat_perm$exp_env_factor)]

# Covariances -- Permutation
cov_est_perm = cov(Cov_matrix_perm$G_means,Cov_matrix_perm$E_means)
true_cov_est_perm = cov(Cov_matrix_perm_ne$G_means,Cov_matrix_perm_ne$E_means)

cor_est_perm = cor(Cov_matrix_perm$G_means,Cov_matrix_perm$E_means)
true_cor_est_perm = cor(Cov_matrix_perm_ne$G_means,Cov_matrix_perm_ne$E_means)

correction_ratio_perm = max(sd(Cov_matrix_perm$E_means),sd(Cov_matrix_perm$G_means))
cov_corrected_perm = round(cov(Cov_matrix_perm$G_means,Cov_matrix_perm$E_means)/(correction_ratio_perm^2),2)

true_correction_ratio_perm = max(sd(Cov_matrix_perm_ne$E_means),sd(Cov_matrix_perm_ne$G_means))
true_cov_corrected_perm = round(cov(Cov_matrix_perm_ne$E_means, Cov_matrix_perm_ne$G_means)/(true_correction_ratio_perm^2),2)

cov_means_perm = cov(Cov_mean_matrix_perm$G_means,Cov_mean_matrix_perm$E_means)
cor_means_perm = cor(Cov_mean_matrix_perm$G_means,Cov_mean_matrix_perm$E_means)

means_correction_perm = max(sd(Cov_mean_matrix_perm$G_means),sd(Cov_mean_matrix_perm$E_means)) 
cov_means_correct_perm = round(cov(Cov_mean_matrix_perm$G_means,Cov_mean_matrix_perm$E_means)/(means_correction_perm^2),2)

# Magnitude of GxE -- EMMs -- Permutation
GxE_emm_perm <- abs(mean(emm_GxE_perm$emmean) - # Overall mean
                    emm_G_perm$emmean[emm_G_perm$gen_factor=="G_1"]- # G
                    emm_E_perm$emmean[emm_E_perm$exp_env_factor=="E_1"]+ # E
                    emm_GxE_perm$emmean[emm_GxE_perm$gen_factor=="G_1" & emm_GxE_perm$exp_env_factor=="E_1"]) # GxE

# Magnitude of GxE -- Loop
allGE_perm <- c()
for (i in 1:nlevels(perm_dat$gen_factor)){
  for (j in 1:nlevels(perm_dat$exp_env_factor)){
    G_levels <- levels(perm_dat$gen_factor)
    E_levels <- levels(perm_dat$exp_env_factor)
    GxE_emm_perm1 <- abs(mean(emm_GxE_perm$emmean) - # Overall mean
                        emm_G_perm$emmean[emm_G_perm$gen_factor == G_levels[i]]- # G
                        emm_E_perm$emmean[emm_E_perm$exp_env_factor == E_levels[j]]+ # E
                        emm_GxE_perm$emmean[emm_GxE_perm$gen_facto == G_levels[i] & emm_GxE_perm$exp_env_factor == E_levels[j]]) # GxE
    allGE_perm <- c(allGE_perm, GxE_emm_perm1)
  }
}
#hist(allGE_perm)
loopGxE_perm = mean(allGE_perm)

# Magnitude of GxE -- EMMs -- Permutation -- No Error
GxE_emm_perm_ne <- abs(mean(emm_GxE_perm_ne$emmean) - # Overall mean
                      emm_G_perm_ne$emmean[emm_G_perm_ne$gen_factor=="G_1"]- # G
                      emm_E_perm_ne$emmean[emm_E_perm_ne$exp_env_factor=="E_1"]+ # E
                      emm_GxE_perm_ne$emmean[emm_GxE_perm_ne$gen_factor=="G_1" & emm_GxE_perm_ne$exp_env_factor=="E_1"]) # GxE

# Magnitude of GxE -- Loop -- EMMs -- Permutation -- No Error
allGE_true_perm <- c()
for (i in 1:nlevels(perm_dat$gen_factor)){
  for (j in 1:nlevels(perm_dat$exp_env_factor)){
    G_levels <- levels(perm_dat$gen_factor)
    E_levels <- levels(perm_dat$exp_env_factor)
    true_GxE_emm_perm <- abs(mean(emm_GxE_perm_ne$emmean) - 
                             emm_G_perm_ne$emmean[emm_G_perm_ne$gen_factor == G_levels[i]] - # phenotype of ith Genotype
                             emm_E_perm_ne$emmean[emm_E_perm_ne$exp_env_factor == E_levels[j]] + # phenotype of jth Environment
                             emm_GxE_perm_ne$emmean[emm_GxE_perm_ne$gen_factor == G_levels[i] & emm_GxE_perm_ne$exp_env_factor == E_levels[j]])
    allGE_true_perm <- c(allGE_true_perm, true_GxE_emm_perm)
  }
}
#hist(allGE_true_perm)
loopGxE_perm_true = mean(allGE_true_perm)

# Magnitude of GxE -- Omega2 -- Permutation
perm_w2_GxE = (summary(aov(test_perm))[[1]][3,2] - #SS_effect -
              (summary(aov(test_perm))[[1]][3,1]*summary(aov(test_perm))[[1]][4,3])) / #Df_effect * MS_error
              (sum(summary(aov(test_perm))[[1]][,2]) + # SS_total+
              (summary(aov(test_perm))[[1]][4,3])) # MS_error

# Magnitude of GxE -- Omega2 -- Permutation -- No Error
perm_w2_GxE_ne = (summary(aov(test_perm_ne))[[1]][3,2] - #SS_effect -
                 (summary(aov(test_perm_ne))[[1]][3,1]*summary(aov(test_perm_ne))[[1]][4,3])) / #Df_effect * MS_error
                 (sum(summary(aov(test_perm_ne))[[1]][,2]) + # SS_total+
                 (summary(aov(test_perm_ne))[[1]][4,3])) # MS_error

# Magnitude of GxE -- Eta2 -- Permutation 
perm_eta_GxE = summary(aov(test_perm))[[1]][3,2]/sum(summary(aov(test_perm))[[1]][,2])

# Magnitude of GxE -- Eta2 -- Permutation -- No Error
perm_eta_GxE_ne = summary(aov(test_perm_ne))[[1]][3,2]/sum(summary(aov(test_perm_ne))[[1]][,2])

# Residual Variation
ResVar_perm = summary(aov(test_perm))[[1]][3,2]/sum(summary(aov(test_perm))[[1]][c(1:3),2])

# Magnitude of GxE -- Means -- Permutation
allGE_means_perm <- c()
for (i in 1:nlevels(perm_means$gen_factor)){
  for (j in 1:nlevels(perm_means$exp_env_factor)){
    G_levels <- levels(perm_means$gen_factor)
    E_levels <- levels(perm_means$exp_env_factor)
    GxE_emm_perm_means <- abs(mean(perm_means$avg_phen_corrected) - # phen of ith and jth
                              mean(perm_means$avg_phen_corrected[perm_means$exp_env_factor == E_levels[j]]) - # phenotype of jth Environment
                              mean(perm_means$avg_phen_corrected[perm_means$gen_factor == G_levels[i]]) + # phenotype of ith Genotype
                              perm_means$avg_phen_corrected[perm_means$gen_factor == G_levels[i] & perm_means$exp_env_factor == E_levels[j]])
    allGE_means_perm <- c(allGE_means_perm, GxE_emm_perm_means)
  }
}
#hist(allGE_means_perm)
GxE_means_perm = mean(allGE_means_perm)

# Magnitude of GxE -- Means -- Permutation -- No error
allGE_means_perm_ne <- c()
for (i in 1:nlevels(perm_means$gen_factor)){
  for (j in 1:nlevels(perm_means$exp_env_factor)){
    G_levels <- levels(perm_means$gen_factor)
    E_levels <- levels(perm_means$exp_env_factor)
    GxE_emm_means_perm_ne <- abs(mean(perm_means$avg_phen_corrected_ne) - # Overall mean
                                 mean(perm_means$avg_phen_corrected_ne[perm_means$gen_factor == G_levels[i]]) - # phenotype of ith Genotype
                                 mean(perm_means$avg_phen_corrected_ne[perm_means$exp_env_factor == E_levels[j]]) + # phenotype of jth Environment
                                 perm_means$avg_phen_corrected_ne[perm_means$gen_factor == G_levels[i] & perm_means$exp_env_factor == E_levels[j]]) # phen of ith and jth
    allGE_means_perm_ne <- c(allGE_means_perm_ne, GxE_emm_means_perm_ne)
  }
}
#hist(allGE_means_perm_ne)
GxE_means_perm_ne = mean(allGE_means_perm_ne)

# Permutation Output
perm_dat. <- data.frame("covariance_perm" = cov_est_perm,
                        "true_covariance_perm" = true_cov_est_perm,
                        "cor_est_perm" = cor_est_perm,
                        "true_cor_est_perm" = true_cor_est_perm,
                        "cov_corrected_perm" = cov_corrected_perm,
                        "true_cov_corrected_perm" = true_cov_corrected_perm,
                        "cov_means_perm" = cov_means_perm,
                        "cor_means_perm" = cor_means_perm,
                        "cov_means_correct_perm"=cov_means_correct_perm,
                        "GxE_mag_perm" = GxE_emm_perm,
                        "true_GxE_mag_perm" = GxE_emm_perm_ne,
                        "GxE_emm_loop_perm" = loopGxE_perm,
                        "true_GxE_emm_loop_perm" = loopGxE_perm_true,
                        "GxE_omega_perm" = perm_w2_GxE,
                        "true_GxE_omega_perm" = perm_w2_GxE_ne,
                        "GxE_eta_perm" = perm_eta_GxE,
                        "true_GxE_eta_perm" = perm_eta_GxE_ne,
                        #"Resid_var_perm" = ResVar_perm,
                        "GxE_means_perm" = GxE_means_perm,
                        "true_GxE_means_perm" = GxE_means_perm_ne)
perm_df <- rbind(perm_df,perm_dat.)
}

## Sanity Check: Null distribution histogram
#ggplot(perm_df, aes(x = Resid_var_perm))+geom_histogram() # Insert variable of interest as x

# Covariance p-value
ptemp1 = (rank(c(cov_est,perm_df$covariance_perm))[1])/(n_boot+1) 
cov_pvalue = NULL
if(ptemp1 < 0.5){cov_pvalue = ptemp1}else{cov_pvalue = (1-ptemp1)} # 2-tailed

# True Covariance p-value
ptemp1_ne = (rank(c(true_cov,perm_df$true_covariance_perm))[1])/(n_boot+1) 
true_cov_pvalue = NULL
if(ptemp1_ne < 0.5){true_cov_pvalue = ptemp1_ne}else{true_cov_pvalue = (1-ptemp1_ne)} # 2-tailed

# Correlation p-value
ptemp2 = (rank(c(cor_est,perm_df$cor_est_perm))[1])/(n_boot+1) 
cor_pvalue = NULL
if(ptemp2 < 0.5){cor_pvalue = ptemp2}else{cor_pvalue = (1-ptemp2)} # 2-tailed

# True Correlation p-value
ptemp2_ne = (rank(c(true_cor,perm_df$cor_est_perm))[1])/(n_boot+1) 
true_cor_pvalue = NULL
if(ptemp2_ne < 0.5){true_cor_pvalue = ptemp2_ne}else{true_cor_pvalue = (1-ptemp2_ne)} # 2-tailed

# Corrected Covariance p-value
ptemp3 = (rank(c(cov_corrected,perm_df$cov_corrected_perm))[1])/(n_boot+1) 
cov_corrected_pvalue = NULL
if(ptemp3 < 0.5){cov_corrected_pvalue = ptemp3}else{cov_corrected_pvalue = (1-ptemp3)} # 2-tailed

# True Corrected Covariance p-value
ptemp3_ne = (rank(c(true_cov_corrected,perm_df$cov_corrected_perm))[1])/(n_boot+1) 
true_cov_corrected_pvalue = NULL
if(ptemp3_ne < 0.5){true_cov_corrected_pvalue = ptemp3_ne}else{true_cov_corrected_pvalue = (1-ptemp3_ne)} # 2-tailed

# Covariance - means - p-value
ptemp4 = (rank(c(cov_est_means,perm_df$cov_means_perm))[1])/(n_boot+1) 
cov_est_means_pvalue = NULL
if(ptemp4 < 0.5){cov_est_means_pvalue = ptemp4}else{cov_est_means_pvalue = (1-ptemp4)} # 2-tailed

# Correlation - means - p-value
ptemp5 = (rank(c(cor_est_means,perm_df$cor_means_perm))[1])/(n_boot+1) 
cor_means_pvalue = NULL
if(ptemp5 < 0.5){cor_means_pvalue = ptemp5}else{cor_means_pvalue = (1-ptemp5)} # 2-tailed

# Corrected covariance - means - p-value
ptemp6 = (rank(c(cov_means_corrected,perm_df$cov_means_correct_perm))[1])/(n_boot+1) 
cov_means_correct_pvalue = NULL
if(ptemp6 < 0.5){cov_means_correct_pvalue = ptemp6}else{cov_means_correct_pvalue = (1-ptemp6)} # 2-tailed

# GxE - EMM - P-value
ptemp7 = (rank(c(GxE_emm_original,perm_df$GxE_mag))[1])/(n_boot+1) 
GxE_pvalue = 1-ptemp7 # Right-tailed

# GxE - EMM - No Error - P-value
ptemp7_ne = (rank(c(true_GxE_original,perm_df$true_GxE_mag_perm))[1])/(n_boot+1) 
true_GxE_pvalue = 1-ptemp7_ne # Right-tailed

# GxE - Omega2 - p-value
ptemp8 = (rank(c(w2_GxE,perm_df$GxE_omega_perm))[1])/(n_boot+1) 
GxE_omega_pvalue = (1-ptemp8) # Right-tailed

# GxE - Omega2 - No error - p-value
ptemp8_ne = (rank(c(true_w2_GxE,perm_df$true_GxE_omega_perm))[1])/(n_boot+1) 
true_GxE_omega_pvalue = (1-ptemp8_ne) # Right-tailed

# GxE - Means - p-value
ptemp9 = (rank(c(GxE_means,perm_df$GxE_means_perm))[1])/(n_boot+1) 
GxE_means_pvalue = (1-ptemp9) # Right-tailed

# GxE - Means - No error - p-value
ptemp9_ne= (rank(c(GxE_means_ne,perm_df$true_GxE_means_perm))[1])/(n_boot+1) 
true_GxE_means_pvalue = (1-ptemp9_ne)

# GxE - Looped Emms
ptemp10 = (rank(c(GxE_emm_loop,perm_df$GxE_emm_loop_perm))[1])/(n_boot+1) 
GxE_loop_pvalue = (1-ptemp10)

# GxE - Looped Emms - no error
ptemp10_ne = (rank(c(GxE_emm_ne_loop,perm_df$true_GxE_emm_loop_perm))[1])/(n_boot+1) 
true_GxE_loop_pvalue = (1-ptemp10_ne)

# GxE - Eta Squared
ptemp11 = (rank(c(eta_GxE,perm_df$GxE_eta_perm))[1])/(n_boot+1) 
GxE_eta_pvalue = (1-ptemp11)

# GxE - Eta Squared - no error
ptemp11_ne = (rank(c(eta_GxE_ne,perm_df$true_GxE_eta_perm))[1])/(n_boot+1) 
GxE_eta_pvalue_ne = (1-ptemp11_ne)

# GxE - Residual Variation
#ptemp12 = (rank(c(ResVar, perm_df$Resid_var_perm))[1])/(n_boot + 1)
#Resid_var_pvalue = (1-ptemp12)

# Model Output
model_info = rbind(mod_df,mod_df_ne)

end.time <- Sys.time()
time.taken <- end.time - start.time

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
                   "true_cov" = true_cov, #Covariance Original
                   "true_cov_pvalue" = true_cov_pvalue,
                   "cov_estimate" = cov_est,
                   "cov_lwrCI" = cov_CI[[1]],
                   "cov_uprCI" = cov_CI[[2]],
                   "cov_pvalue" = cov_pvalue,  
                   
                   "true_cor" = true_cor, #Correlation Original
                   "true_cor_pvalue" = true_cor_pvalue,
                   "cor_estimate" = cor_est,
                   "cor_lwrCI" = cor_CI[[1]],
                   "cor_uprCI" = cor_CI[[2]],
                   "cor_pvalue" = cor_pvalue,
                   
                   "true_cov_corrected" = true_cov_corrected, #Corrected Covariance Estimates
                   "true_cov_corrected_pvalue" = true_cov_corrected_pvalue,
                   "cov_corrected" = cov_corrected, 
                   "cov_corrected_lwrCI" = cov_corrected_CI[[1]],
                   "cov_corrected_uprCI" = cov_corrected_CI[[2]],
                   "cov_corrected_pvalue" = cov_corrected_pvalue,
                   
                   "true_cov_means" = true_cov_means, # Covariance estimate on means (no pvalue on true means)
                   "cov_means" = cov_est_means, 
                   "cov_means_lwrCI" = cov_means_CI[[1]],
                   "cov_means_uprCI" = cov_means_CI[[2]],
                   "cov_means_pvalue" = cov_est_means_pvalue,
                   
                   "true_cor_means" = true_cor_means, # Correlation estimate on means 
                   "cor_means" = cor_est_means,
                   "cor_means_lwrCI" = cor_means_CI[[1]],
                   "cor_means_uprCI" = cor_means_CI[[2]],
                   "cor_means_pvalue" = cor_means_pvalue,
                   
                   "true_cov_means_correct" = true_cov_means_corrected, # Corrected Covariance on means
                   "cov_means_correct" = cov_means_corrected,
                   "cov_means_correct_lwrCI" = cov_means_correct_CI[[1]],
                   "cov_means_correct_uprCI" = cov_means_correct_CI[[2]],
                   "cov_means_correct_pvalue" = cov_means_correct_pvalue)

GxE <- data.frame("row" = row,
                  "true_GxE" = round(true_GxE_original,2), # Emmeans GxE 
                  "true_GxE_pvalue" = round(true_GxE_pvalue,2),
                  "GxE_emm_estimate" = round(GxE_emm_original,2),
                  "GxE_emm_lwrCI" = round(GxE_CI[[1]],2),
                  "GxE_emm_uprCI" = round(GxE_CI[[2]],2),
                  "GxE_emm_pvalue" = round(GxE_pvalue,2),
                     
                  "GxE_emm_loop" = round(GxE_emm_loop,2), # GxE Emmeans from loop
                  "true_GxE_emm_loop" = round(GxE_emm_ne_loop,2),
                  "true_GxE_emm_loop_pvalue" = round(true_GxE_loop_pvalue,2),
                  "GxE_emm_loop_lwrCI" = round(GxE_loop_CI[[1]],2),
                  "GxE_emm_loop_uprCI" = round(GxE_loop_CI[[2]],2),
                  "GxE_emm_loop_pvalue" = round(GxE_loop_pvalue,2),
                 
                  "true_omega_GxE" = round(true_w2_GxE,2), # Omega^2 GxE 
                  "true_omega_GxE_pvalue" = round(true_GxE_omega_pvalue,2),
                  "GxE_omega_estimate" = round(w2_GxE,2), 
                  "GxE_omega_lwrCI" = round(GxE_omega_CI[[1]],2),
                  "GxE_omega_uprCI" = round(GxE_omega_CI[[2]],2),
                  "GxE_omega_pvalue" = round(GxE_omega_pvalue,2), 
                 
                  "true_eta_GxE" = round(eta_GxE_ne,2),
                  "true_eta_GxE_pvalue" = round(GxE_eta_pvalue_ne,2),
                  "GxE_eta" = round(eta_GxE,2),
                  "GxE_eta_lwrCI" = round(GxE_eta_CI[[1]],2),
                  "GxE_eta_uprCI" = round(GxE_eta_CI[[2]],2),
                  "GxE_eta_pvalue" = round(GxE_eta_pvalue,2),
              
                  "True_resid_variation" = round(ResVar_ne,2),   
                  "Resid_variation" = round(ResVar,2),
                  "GxE_residVar_lwrCI" = round(GxE_resvar_CI[[1]],2),
                  "GxE_residVar_uprCI" = round(GxE_resvar_CI[[2]],2),
                  #"Resid_variation_pvalue" = Resid_var_pvalue,
                 
                  "true_GxE_means" = round(GxE_means_ne,2), # GxE on means
                  "true_GxE_means_pvalue"= round(true_GxE_means_pvalue,2),
                  "GxE_means" = round(GxE_means,2),
                  "GxE_means_lwrCI" = round(GxE_means_CI[[1]],2),
                  "GxE_means_uprCI" = round(GxE_means_CI[[2]],2),
                  "GxE_means_pvalue" = round(GxE_means_pvalue,2)
                  ) 

# Output
write.csv(GxE,paste0("/scratch/albecker/Power_analysis/power_output/GxE_",row,"_output.csv"))
write.csv(Covariance,paste0("/scratch/albecker/Power_analysis/power_output/Covariance_",row,"_output.csv"))
write.csv(Parameters,paste0("/scratch/albecker/Power_analysis/power_output/Parameters_",row,"_output.csv"))

write.csv(phen_out,paste0("/scratch/albecker/Power_analysis/phenotype_output/Phenotype_data",row,"_output.csv"))
write.csv(perm_df,paste0("/scratch/albecker/Power_analysis/permutation_output/Permutation_data",row,"_output.csv"))
write.csv(boot_df,paste0("/scratch/albecker/Power_analysis/bootstrap_output/Bootstrap_data",row,"_output.csv"))
write.csv(GEmeans_out,paste0("/scratch/albecker/Power_analysis/GEmeans_output/GEmeans_data",row,"_output.csv"))
write.csv(model_info,paste0("/scratch/albecker/Power_analysis/Anova_output/model_info_data",row,"_output.csv"))


#####################
### DATA PLOTTING ###
#####################

# Compile Data
require(readr)
#Param_temp <- list.files(path = "~/Desktop/power_output/",pattern= "*Parameter") 
#Cov_temp <-list.files(path = "~/Desktop/power_output/",pattern= "*Covariance") 
#GxE_temp <- list.files(path = "~/Desktop/power_output/",pattern= "*GxE") 
#setwd("~/Desktop/power_output/")
#GE = plyr::ldply(GxE_temp, read_csv)
#Cov = plyr::ldply(Cov_temp, read_csv)
#Par = plyr::ldply(Param_temp, read_csv)

#dat_csv. = full_join(GE, Cov, by= "row")
#dat_csv = full_join(dat_csv., Par, by = "row")

#write.csv(dat_csv,"~/Desktop/dat_csv.csv")
dat_csv <- read.csv("~/Desktop/dat_csv.csv")

#temp2 <- list.files(path = "~/Desktop/Simulation_output/phenotype_output/",pattern= "*Phenotype_data")
#setwd("~/Desktop/Simulation_output/phenotype_output/")
#phendf = plyr::ldply(temp2, read_csv)
#write.csv(phendf,"~/Desktop/phendf.csv")
#phendf <- read.csv("~/Desktop/phendf.csv")

#temp3 <- list.files(path = "~/Desktop/Power_analysis_output/GEmeans_output/",pattern= "*GEmeans_data")
#setwd("~/Desktop/Power_analysis_output/GEmeans_output/")
#gemeans = plyr::ldply(temp3, read_csv)

dat_csv$col = NULL
# Color assignment code
for(i in 1:nrow(dat_csv)){
  if(dat_csv$covariance_pvalue[i] <= 0.05 & dat_csv$GxE_emm_pvalue[i] <= 0.05){dat_csv$col[i] = "red" # Both significant
  }else if(dat_csv$covariance_pvalue[i] <= 0.05 & dat_csv$GxE_emm_pvalue[i] > 0.05){dat_csv$col[i] = "darkgreen" # Cov significant
  }else if(dat_csv$covariance_pvalue[i] > 0.05 & dat_csv$GxE_emm_pvalue[i] <= 0.05){dat_csv$col[i] = "dodgerblue4" # GxE significant
  }else{dat_csv$col[i] = "grey"} # None significant
}

# Summarize by proportion of significant
dat_csv$cull = NULL
for(i in 1:nrow(dat_csv)){
  if(dat_csv$cov_corrected[i] >= 0.4 & dat_csv$cov_corrected[i]<=0.6){dat_csv$cull[i] = TRUE
  }else if(dat_csv$GxE_emm_estimate[i] >= 0.4 & dat_csv$GxE_emm_estimate[i] <= 0.6){dat_csv$cull[i] = TRUE
    }else{dat_csv$cull[i] = FALSE}
}
col_df = dat_csv %>%
  filter(cull == TRUE) %>%
  group_by(n_pop,sample_size,col) %>%
  summarize(frequency = n())

new_df = data.frame()
for(i in unique(col_df$sample_size)){
  for(j in unique(col_df$n_pop)){
  subset = col_df[col_df$sample_size==i,]
  subsetsubset = subset[subset$n_pop == j,]
  total = sum(subsetsubset$frequency)
  colproportion = subsetsubset$frequency/total
  new_df. = data.frame(subsetsubset,"proportion"=colproportion)
  new_df = rbind(new_df,new_df.)
  }
}
red_dat = filter(new_df,col=="red")
ggplot(red_dat,aes(x = factor(n_pop), y = proportion,group = factor(sample_size), colour = factor(sample_size)))+
  geom_point(size = 4)+geom_line()+#theme_classic()+
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(axis.text.x = element_text(size=16,colour = "black"),
        axis.title.x = element_text(size=18,face="bold")) +
  theme(axis.text.y = element_text(size=16,colour = "black"),
        axis.title.y = element_text(size=18,face="bold")) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2))+
  xlab("Number of Populations")+ylab("Proportion of significant CovGE and GxE (alpha = 0.05)")+
  labs(colour = "Sample Size")

blue_dat = filter(new_df,col=="dodgerblue4")
ggplot(blue_dat,aes(x = factor(n_pop), y = proportion,group = factor(sample_size), colour = factor(sample_size)))+
  geom_point(size = 4)+geom_line()+#theme_classic()+
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(axis.text.x = element_text(size=16,colour = "black"),
        axis.title.x = element_text(size=18,face="bold")) +
  theme(axis.text.y = element_text(size=16,colour = "black"),
        axis.title.y = element_text(size=18,face="bold")) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2))+
  xlab("Number of Populations")+ylab("Proportion of significant GxE (alpha = 0.05)")+
  labs(colour = "Sample Size")

green_dat = filter(new_df,col=="darkgreen")
ggplot(green_dat,aes(x = factor(n_pop), y = proportion,group = factor(sample_size), colour = factor(sample_size)))+
  geom_point(size = 4)+geom_line()+#theme_classic()+
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(axis.text.x = element_text(size=16,colour = "black"),
        axis.title.x = element_text(size=18,face="bold")) +
  theme(axis.text.y = element_text(size=16,colour = "black"),
        axis.title.y = element_text(size=18,face="bold")) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2))+
  xlab("Number of Populations")+ylab("Proportion of significant CovGE (alpha = 0.05)")+
  labs(colour = "Sample Size")

## Sanity Check: Makesure CI for estimates overlaps with TruCov and TruGxE
dat_csv$probgxe = NULL
dat_csv$probcov = NULL

for(i in 1:nrow(dat_csv)){
  if((dat_csv$true_GxE[i] < dat_csv$GxE_emm_lwrCI[i]) | (dat_csv$true_GxE[i] >dat_csv$GxE_emm_uprCI[i])){dat_csv$probgxe[i] = TRUE
  }else{dat_csv$probgxe[i] = FALSE}
  if((dat_csv$true_cov_corrected[i] < dat_csv$cov_corrected_lwrCI[i]) | (dat_csv$true_cov_corrected[i] >dat_csv$cov_corrected_uprCI[i])){dat_csv$probcov[i] = TRUE
  }else{dat_csv$probcov[i] = FALSE}
}

GxEanoms = dat_csv %>%
  filter(probgxe == TRUE)
dim(GxEanoms)
length(which(GxEanoms$GxE_emm_pvalue<=0.05))

cov_anom = dat_csv %>%
  filter(probcov ==TRUE)
length(which(cov_anom$cov_corrected_pvalue<=0.05))

overlapper_cov = ggplot(cov_anom,aes(x = row))+ theme_classic() + 
  geom_errorbar(aes(ymin = cov_corrected_lwrCI,ymax = cov_corrected_uprCI),color = "black")+
  geom_point(aes(y = true_cov_corrected),size = 2, color = "red",alpha = 0.5)+
  xlab("Unique Parameter Set (row)") + ylab("Covariance")
  overlapper_cov

overlapper_GxE = ggplot(GxEanoms,aes(x = row))+ theme_classic() + 
    geom_errorbar(aes(ymin = GxE_emm_lwrCI,ymax = GxE_emm_uprCI),color = "black")+
    geom_point(aes(y = true_GxE),size = 2, color = "red",alpha = 0.5)+
    xlab("Unique Parameter Set (row)") + ylab("GxE - Estimated Marginal Mean")
overlapper_GxE

## Sanity check: Plot means against raw estimates... should fall along 1:1 line
dat_csv$meancoverror = abs(dat_csv$cov_means_correct_uprCI - dat_csv$cov_means_correct_lwrCI)
dat_csv$coverror = abs(dat_csv$covariance_uprCI - dat_csv$covariance_lwrCI)
dat_csv$meangxeerror = dat_csv$GxE_means_uprCI - dat_csv$GxE_means_lwrCI
dat_csv$gxeerror = dat_csv$GxE_emm_uprCI - dat_csv$GxE_emm_lwrCI

(covmeancheck = ggplot(dat_csv,aes(x = true_cov, y = true_cov_means_correct))+
    geom_point()+theme_classic()+ylab("Covariance from Means")+xlab("Covariance from Raw")+
    geom_abline(slope = 1, intercept = 0,colour = "red"))
(covcheck = ggplot(dat_csv,aes(x = true_cov, y = covariance))+
    geom_point()+theme_classic()+ylab("Covariacne")+xlab("True Covariance")+
    geom_abline(slope = 1, intercept = 0,colour = "red"))
suspect = dat_csv %>%
  filter(true_cov == 1)

(gxemeancheck = ggplot(dat_csv,aes(x = true_GxE_emm, y = true_GxE_means))+
    geom_point()+theme_classic()+ylab("GxE from Means")+xlab("GxE from Raw")+
    geom_abline(slope = 1, intercept = 0,colour = "red"))

(coverrorcheck = ggplot(dat_csv, aes(x = coverror,y = meancoverror)) + 
    geom_point(alpha = 0.5) + theme_classic() + ylab("Length of 95% CI for CovGE means") + xlab("Length of 95% CI for CovGE raw")+
    geom_abline(slope = 1, intercept = 0,colour = "red")+scale_colour_identity())
(gxeerrorcheck = ggplot(dat_csv, aes(x = gxeerror,y = meangxeerror)) + 
    geom_point(alpha = 0.3) + theme_classic() + ylab("Length of 95% CI for GxE means") + xlab("Length of 95% CI for GxE raw")+
    geom_abline(slope = 1, intercept = 0,colour = "red")+scale_colour_identity())

dat_csv$testcol = NULL

# Color assignment code FOR TEST
for(i in 1:nrow(dat_csv)){
  if(dat_csv$meangxeerror[i] <= 0.05 & dat_csv$gxeerror[i] <= 0.05){dat_csv$testcol[i] = "red" # Both significant
  }else if(dat_csv$meangxeerror[i] <= 0.05 & dat_csv$gxeerror[i] > 0.05){dat_csv$testcol[i] = "darkgreen" # Just means significant
  }else if(dat_csv$meangxeerror[i] > 0.05 & dat_csv$gxeerror[i] <= 0.05){dat_csv$testcol[i] = "dodgerblue4" # Just raw is significant
  }else{dat_csv$testcol[i] = "grey"} # None significant
}
    
# Omega^2 value against significance:
dat_csv$omega_colour = NULL
for(i in 1:nrow(dat_csv)){
  if(dat_csv$GxE_omega_pvalue[i] <= 0.05){dat_csv$omega_colour[i] = "purple" # Both significant
  }else{dat_csv$omega_colour[i] = "grey"} # None significant
}
ggplot(dat_csv, aes(x = GxE_emm_estimate, y = GxE_omega_estimate))+geom_point(colour = "black",alpha = 0.3)+theme_classic()#+scale_colour_identity()

tempdf$exp_env_factor <- factor(tempdf$exp_env_factor, levels =c("E_1","E_2","E_3","E_4","E_5","E_6","E_7","E_8","E_9","E_10","E_11","E_12","E_13","E_14","E_15"))

plotspot.xaxis <- c("E_1" = "1","E_2" = "2", "E_3" = "3","E_4" = "4","E_5" = "5",
                "E_6" = "6","E_7" = "7", "E_8" = "8","E_9" = "9","E_10" = "10",
                "E_11" = "11","E_12" = "12", "E_13" = "13","E_14" = "14","E_15" = "15")
nb.cols <- 15
mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(nb.cols)
ggplot(dat_csv, aes(x= interaction , y = GxE_emm_estimate))+ geom_point()
biglittle = dat_csv %>% #Filter to small EMM but large Omegas. 
  filter(GxE_emm_estimate < 0.1) %>%
  filter(GxE_omega_estimate > 0.8)

tempdf. = filter(phendf, row == 11480) 
tempdf = left_join(tempdf., biglittle[biglittle$row == 11480,],by= "row")
tempdf$envgen = paste(tempdf$gen_factor,tempdf$exp_env_factor,by = "_")
tempplot = ggplot(tempdf, aes(x = exp_env_factor,y = phen_corrected, group = envgen, colour = gen_factor))+
  geom_boxplot()+
  #geom_line(aes(group = gen_factor))+
  scale_colour_manual(values = mycolors)+
    theme_classic()+
    scale_x_discrete(labels = plotspot.xaxis)+
    ggtitle(paste0("Row: ", unique(tempdf$row)))+
  xlab("Environment")+ylab("Phenotype (standardized)")+
    annotate("text", x = "E_1", y = 2.6, 
             label = paste0("Omega2 = ", unique(tempdf$GxE_omega_estimate),
                            ", p = ", unique(tempdf$GxE_omega_pvalue),
                            ", Emm = ",unique(tempdf$GxE_emm_estimate), 
                            ", p = ", unique(tempdf$GxE_emm_pvalue)), size = 4, hjust = 0)
tempplot

# New facet label names for supp variable
sample.labs <- c("2" = "2 Samples", "3" = "3 Samples", "5"= "5 Samples","10"="10 Samples","15"="15 Samples")
pop.labs <- c("5"="5 Populations", "10"="10 Populations", "20"="20 Populations")


## Corrected Covariance by GxE
require(ggplot2)
ggplot(dat_csv, aes(x = covariance, y = GxE_emm, group = factor(n_pop), alpha = 0.1,colour = col)) + 
  geom_jitter() + theme_classic() + ylim(0,1) + xlim(-1,1)+
  xlab("Covariance Estimate") + ylab("GxE Estimate") +
  theme(legend.position = "none")+
  scale_colour_identity()+
  facet_grid(sample_size~n_pop)#,
             #labeller=labeller(sample_size = as_labeller(sample.labs),
                              # n_pop = as_labeller(pop.labs)))

## Corrected Covariance by Omege
require(ggplot2)
ggplot(dat_csv, aes(x = covariance, y = GxE_omega, group = factor(n_pop), alpha = 0.1,colour = col)) + 
  geom_point() + theme_classic() + #ylim(0,1.2)+ xlim(-1,1)+
  xlab("Covariance Estimate") + ylab("GxE Estimate (Omega^2)") +
  theme(legend.position = "none")+
  scale_colour_identity()+
  facet_grid(sample_size~n_pop)#,

## Relationship between GxE mag and Covariance (Is there a tradeoff? Yarp.)
sigGxE = filter(dat_csv, GxE_emm_pvalue <=0.05 | covariance_pvalue <= 0.05)
ggplot(sigGxE, aes(x = GxE_emm, y = covtick))+
  geom_smooth(method = "glm",method.args = list(family = "binomial"),se = T,colour = "black") + 
  xlab("Magnitude of GxE")+ylab("Proportion of significant CovGE values (p < 0.05)")+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(axis.text.x = element_text(size=16,colour = "black"),
        axis.title.x = element_text(size=18,face="bold")) +
  theme(axis.text.y = element_text(size=16,colour = "black"),
        axis.title.y = element_text(size=18,face="bold")) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2))

# Hex plot
ggplot(dat_csv, aes(x = covariance, y = GxE_emm)) + 
  geom_hex()+
    theme_classic() +# ylim(0,1.2)+ xlim(-1,1)+
  xlab("Covariance Estimate") + ylab("GxE Estimate") +
  ggtitle("HexPlot")+
  #theme(legend.position = "none")+
  #scale_colour_identity()+ 
  facet_grid(sample_size~n_pop)#,
#labeller=labeller(sample_size = as_labeller(sample.labs),
#n_pop = as_labeller(pop.labs)))

## Difference in covariance measures
dat_csv$cov_diff = dat_csv$true_cov-dat_csv$true_cov_corrected
ggplot(dat_csv, aes(x = true_GxE, y = cov_diff, colour = col)) + 
  geom_point() + theme_classic() + #ylim(0,3)+ xlim(-1,1)+
  xlab("True Covariance Corrected") + ylab("Difference in covariance (Uncorrected - Corrected)") +
  scale_colour_identity()+
  facet_grid(sample_size~n_pop)

## Difference in covariance measures
dat_csv$means_diff = dat_csv$true_cov_corrected-dat_csv$true_cov_means_correct
ggplot(dat_csv,aes(x = means_diff))+geom_histogram()+theme_classic()+xlab("Difference in Cov (Raw - Means)")+ylab("Count")

dat_csv$gxe_means_diff = dat_csv$true_GxE-dat_csv$true_GxE_means
ggplot(dat_csv,aes(x = gxe_means_diff))+geom_histogram()+theme_classic()+xlab("Difference in GxE (Raw - Means)")+ylab("Count")

# Plot phenotype according to environment 

setwd("~/Desktop/Simulation_output/phenotype_output/plots/")

tempData = dat_csv %>% # High GxE; CnGV
  filter(sample_size == 5) %>%
  filter(std_dev == 0.5) %>%
  filter(n_pop == 15) %>%
  #filter(delta_gen == -1) %>%
  #filter(delta_env == 1) %>%
  filter(interaction == 15) %>%
  filter(GxE_emm_pvalue <= 0.05) %>%
  filter(cov_corrected_pvalue <= 0.05) %>%
  filter(cov_corrected > 0)

if(length(unique(tempData$row))>1){
tempData=tempData[1,]}else{tempData = tempData}

sumdf = merge(phendf,tempData, by = "row")
plotdf = sumdf

plotspot.x = plotspot.y = NULL
mycolors <- NULL
  if(unique(plotdf$n_pop.x) == 2){
    plotdf$exp_env_factor <- factor(plotdf$exp_env_factor, levels =c("E_1","E_2"))
    plotspot.x <- c("E_1" = "1","E_2" = "2")
    plotspot.y <- c("G_1" = "Genotype 1","G_2" = "Genotype 2")
    nb.cols <- 2
    mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(nb.cols)
  }else if(unique(plotdf$n_pop.x) == 3){
    plotdf$exp_env_factor <- factor(plotdf$exp_env_factor, levels =c("E_1","E_2","E_3"))
    plotspot.x <- c("E_1" = "1","E_2" = "2", "E_3" = "3")
    plotspot.y <- c("G_1" = "Genotype 1","G_2" = "Genotype 2","G_3" = "Genotype 3")
    nb.cols <- 3
    mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(nb.cols)
  }else if(unique(plotdf$n_pop.x) == 5){
    plotdf$exp_env_factor <- factor(plotdf$exp_env_factor, levels =c("E_1","E_2","E_3","E_4","E_5"))
    plotspot.x <- c("E_1" = "1","E_2" = "2", "E_3" = "3","E_4" = "4","E_5" = "5")  
    plotspot.y <- c("G_1" = "Genotype 1","G_2" = "Genotype 2","G_3" = "Genotype 3","G_4" = "Genotype 4","G_5" = "Genotype 5")
    nb.cols <- 5
    mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(nb.cols)
  }else if(unique(plotdf$n_pop.x) == 10){
    plotdf$exp_env_factor <- factor(plotdf$exp_env_factor, levels =c("E_1","E_2","E_3","E_4","E_5","E_6","E_7","E_8","E_9","E_10"))
    plotdf$gen_factor <- factor(plotdf$gen_factor, levels =c("G_1","G_2","G_3","G_4","G_5","G_6","G_7","G_8","G_9","G_10"))
    plotspot.x <- c("E_1" = "1","E_2" = "2", "E_3" = "3","E_4" = "4","E_5" = "5",
                    "E_6" = "6","E_7" = "7", "E_8" = "8","E_9" = "9","E_10" = "10")
    plotspot.y <- c("G_1" = "Genotype 1","G_2" = "Genotype 2","G_3" = "Genotype 3","G_4" = "Genotype 4","G_5" = "Genotype 5",
                    "G_6" = "Genotype 6","G_7" = "Genotype 7","G_8" = "Genotype 8","G_9" = "Genotype 9","G_10" = "Genotype 10")
    nb.cols <- 10
    mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(nb.cols)
  }else{
    plotdf$exp_env_factor <- factor(plotdf$exp_env_factor, levels =c("E_1","E_2","E_3","E_4","E_5","E_6","E_7","E_8","E_9","E_10","E_11","E_12","E_13","E_14","E_15"))
    plotdf$gen_factor <- factor(plotdf$gen_factor, levels =c("G_1","G_2","G_3","G_4","G_5","G_6","G_7","G_8","G_9","G_10","G_11","G_12","G_13","G_14","G_15"))
    plotspot.x <- c("E_1" = "1","E_2" = "2", "E_3" = "3","E_4" = "4","E_5" = "5",
                    "E_6" = "6","E_7" = "7", "E_8" = "8","E_9" = "9","E_10" = "10",
                    "E_11" = "11","E_12" = "12", "E_13" = "13","E_14" = "14","E_15" = "15")
    plotspot.y <- c("G_1" = "Genotype 1","G_2" = "Genotype 2","G_3" = "Genotype 3","G_4" = "Genotype 4","G_5" = "Genotype 5",
                    "G_6" = "Genotype 6","G_7" = "Genotype 7","G_8" = "Genotype 8","G_9" = "Genotype 9","G_10" = "Genotype 10",
                    "G_11" = "Genotype 11","G_12" = "Genotype 12","G_13" = "Genotype 13","G_14" = "Genotype 14","G_15" = "Genotype 15")
    nb.cols <- 15
    mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(nb.cols)
  }
  
  #require(RColorBrewer)
  #require(ggplot2)
  plot_temp = ggplot(plotdf, aes(x = exp_env_factor, y = phen_corrected, group = gen_factor, colour = gen_factor))+
    geom_point()+geom_smooth(se = FALSE)+theme_classic()+ylim(-3,3)+
    scale_colour_manual(values = mycolors,
                       labels = plotspot.y)+
    #ggtitle(paste0("Indentifying row:",unique(plotdf$row))) +
    xlab("Environment") + ylab("Standardized Phenotype")+labs(colour = " ")+
    annotate("text", x = "E_1", y = 2.9, 
             label = paste0("Cov = ", unique(plotdf$cov_corrected),", p = ", round(unique(plotdf$cov_corrected_pvalue),2)), size = 8, hjust = 0)+
    annotate("text", x = "E_1", y = 2.4, 
             label = paste0("GxE = ",round(unique(plotdf$GxE_emm_estimate),2),", p = ",round(unique(plotdf$GxE_emm_pvalue),2)), size = 8, hjust = 0)+
    theme(plot.margin = margin(0, 0, 0, 0, "pt"))+
    scale_x_discrete(labels = plotspot.x)+
    theme(axis.text.x = element_text(size=16,colour = "black"),
          axis.title.x = element_text(size=20,face="bold")) +
    theme(axis.text.y = element_text(size=16,colour = "black"),
          axis.title.y = element_text(size=20,face="bold"))+
    theme(legend.text = element_text(size = 16))
  plot_temp

## Power Analysis
dat_csv$covtick <- NULL
for(i in 1:nrow(dat_csv)){
  if(dat_csv$cov_corrected_pvalue[i] > 0.05){dat_csv$covtick[i]=0}else{dat_csv$covtick[i]=1} 
}

dat_csv$gxetick <- NULL
for(i in 1:nrow(dat_csv)){
  if(dat_csv$GxE_emm_pvalue[i] > 0.05){dat_csv$gxetick[i]=0}else{dat_csv$gxetick[i]=1}
}


# Power of across all
require(tidyverse)
pow = dat_csv %>%
  group_by(delta_env, delta_gen, sample_size,n_pop,std_dev,interaction) %>%
  summarize("count_tick" = n(),
            "sum_covtick" = sum(covtick),
            "sum_gxetick" = sum(gxetick))
pow$covpower = pow$sum_covtick/pow$count_tick
pow$gxepower = pow$sum_gxetick/pow$count_tick

# Power at n_pop = 2 (Claim by Conover and Schultz)
claim = pow %>%
  filter(n_pop == 2) %>%
  group_by(sample_size,std_dev) %>%
  summarize(mean(covpower))
claim # Not worth plotting

# Heatmap of all values
totalpowlow = pow %>%
  filter(std_dev == 0.5)%>%
  group_by(n_pop, sample_size)%>%
  summarize("meancovpower" = mean(covpower),
            "meangxepower" = mean(gxepower))

totalpowhigh = pow %>%
  filter(std_dev == 1)%>%
  group_by(n_pop, sample_size)%>%
  summarize("meancovpower" = mean(covpower),
            "meangxepower" = mean(gxepower))

require(gridExtra)
covlow = ggplot(totalpowlow,aes(x = factor(sample_size), y = factor(n_pop), fill = meancovpower)) + 
  geom_tile() + scale_fill_gradient(low="grey", high="darkgreen") +
  xlab("Sample Size") + ylab("Number of Populations")+
  labs(fill = "Power")+
  theme_classic() + ggtitle("Covariance: Standard Deviation = 1")
covhigh = ggplot(totalpowhigh,aes(x = factor(sample_size), y = factor(n_pop), fill = meancovpower)) + 
  geom_tile() + scale_fill_gradient(low="grey", high="darkgreen") +
  xlab("Sample Size") + ylab("Number of Populations")+
  labs(fill = "Power")+
  theme_classic()+ggtitle("Covariance: Standard Deviation = 0.5")
grid.arrange(covlow,covhigh)

gxelow = ggplot(totalpowlow,aes(x = factor(sample_size), y = factor(n_pop), fill = meangxepower)) + 
  geom_tile() + scale_fill_gradient(low="grey", high="dodgerblue4") +
  xlab("Sample Size") + ylab("Number of Populations")+
  labs(fill = "Power")+
  theme_classic() + ggtitle("GxE: Standard Deviation = 1")
gxehigh = ggplot(totalpowhigh,aes(x = factor(sample_size), y = factor(n_pop), fill = meangxepower)) + 
  geom_tile() + scale_fill_gradient(low="grey", high="dodgerblue4") +
  xlab("Sample Size") + ylab("Number of Populations")+
  labs(fill = "Power")+
  theme_classic()+ggtitle("GxE: Standard Deviation = 0.5")
grid.arrange(gxehigh,gxelow,covhigh,covlow,ncol = 2)


# Filter to just values between 0.4 and 0.6 (this window captures full power)
covpow1 = dat_csv %>%
  filter(between(cov_corrected, 0.4,0.6))

gxepow1 = dat_csv %>%
  filter(between(GxE_emm_estimate,0.4,0.6))

covpow = covpow1 %>%
  group_by(delta_env, delta_gen, sample_size,n_pop,std_dev,interaction) %>%
  summarize("count_tick" = n(),
            "sum_covtick" = sum(covtick))
covpow$covpower = covpow$sum_covtick/covpow$count_tick

gxepow = gxepow1 %>%
  group_by(delta_env, delta_gen, sample_size,n_pop,std_dev,interaction) %>%
  summarize("count_tick" = n(),
            "sum_gxetick" = sum(gxetick))
gxepow$gxepower = gxepow$sum_gxetick/gxepow$count_tick

# Divide according to standard deviation
covpow_sd.5 = covpow %>%
  filter(std_dev == 0.5) %>%
  group_by(sample_size, n_pop) %>%
  summarize("meancovpower" = mean(covpower))

covpow_sd1 = covpow %>%
  filter(std_dev == 1) %>%
  group_by(sample_size, n_pop) %>%
  summarize("meancovpower" = mean(covpower))

gxepow_sd.5 = gxepow %>%
  filter(std_dev == 0.5) %>%
  group_by(sample_size, n_pop) %>%
  summarize("meangxepower" = mean(gxepower))
  
gxepow_sd1 = gxepow %>%
  filter(std_dev == 1) %>%
  group_by(sample_size, n_pop) %>%
  summarize("meangxepower" = mean(gxepower))

covpow_sd.5$pop_samp = paste(covpow_sd.5$n_pop,covpow_sd.5$sample_size)
covpow_sd1$pop_samp = paste(covpow_sd1$n_pop,covpow_sd1$sample_size)
gxepow_sd.5$pop_samp = paste(gxepow_sd.5$n_pop,gxepow_sd.5$sample_size)
gxepow_sd1$pop_samp = paste(gxepow_sd1$n_pop,gxepow_sd1$sample_size)

require(gridExtra)
csdlow = ggplot(covpow_sd.5,aes(x = factor(sample_size), y = factor(n_pop), fill = meancovpower)) + 
  geom_tile() + scale_fill_gradient(low="grey", high="darkgreen") +
  xlab("Sample Size") + ylab("Number of Populations")+
  labs(fill = "Power")+
  theme_classic() + ggtitle("Covariance: Standard Deviation = 1")
csdhigh = ggplot(covpow_sd1,aes(x = factor(sample_size), y = factor(n_pop), fill = meancovpower)) + 
  geom_tile() + scale_fill_gradient(low="grey", high="darkgreen") +
  xlab("Sample Size") + ylab("Number of Populations")+
  labs(fill = "Power")+
  theme_classic()+ggtitle("Covariance: Standard Deviation = 0.5")
grid.arrange(csdlow,csdhigh)

gsdlow = ggplot(gxepow_sd.5,aes(x = factor(sample_size), y = factor(n_pop), fill = meangxepower)) + 
  geom_tile() + scale_fill_gradient(low="grey", high="dodgerblue4") +
  xlab("Sample Size") + ylab("Number of Populations")+
  labs(fill = "Power")+
  theme_classic() + ggtitle("GxE: Standard Deviation = 1")
gsdhigh = ggplot(gxepow_sd1,aes(x = factor(sample_size), y = factor(n_pop), fill = meangxepower)) + 
  geom_tile() + scale_fill_gradient(low="grey", high="dodgerblue4") +
  xlab("Sample Size") + ylab("Number of Populations")+
  labs(fill = "Power")+
  theme_classic()+ggtitle("GxE: Standard Deviation = 0.5")
grid.arrange(gsdhigh,gsdlow,csdhigh,csdlow,ncol = 2)

## Goof around with sims with same total n of diff combos of sample size and n_pop

hundie1 = dat_csv %>%
  filter(sample_size == 10 & n_pop == 10)
hundie2 = dat_csv %>%
  filter(sample_size == 20 & n_pop == 5)
hundie = rbind(hundie1,hundie2)  

hundiepow = hundie %>%
  group_by(delta_env, delta_gen, sample_size,n_pop,std_dev,interaction) %>%
  summarize("count_tick" = n(),
            "sum_covtick" = sum(covtick),
            "sum_gxetick" = sum(gxetick))
hundiepow$covpower = hundiepow$sum_covtick/hundiepow$count_tick
hundiepow$gxepower = hundiepow$sum_gxetick/hundiepow$count_tick

fitty1 = dat_csv %>%
  filter(sample_size == 5 & n_pop == 10)
fitty2 = dat_csv %>%
  filter(sample_size == 10 & n_pop == 5)
fitty = rbind(fitty1,fitty2)  
fittypow = fitty %>%
  group_by(delta_env, delta_gen, sample_size,n_pop,std_dev,interaction) %>%
  summarize("count_tick" = n(),
            "sum_covtick" = sum(covtick),
            "sum_gxetick" = sum(gxetick))
fittypow$covpower = fittypow$sum_covtick/fittypow$count_tick
fittypow$gxepower = fittypow$sum_gxetick/fittypow$count_tick

## Corrected Covariance by GxE
ggplot(fitty, aes(x = cov_corrected, y = GxE_emm_estimate, group = factor(n_pop), alpha = 0.3,colour = col)) + 
  geom_point() + theme_classic() + ylim(0,3)+ xlim(-1,1)+
  xlab("Covariance Estimate") + ylab("GxE Estimate") +
  theme(legend.position = "none")+
  scale_colour_identity()+
  facet_grid(sample_size~n_pop)#,

# Heatmaps
hundielowsum = hundiepow %>%
  filter(std_dev == 0.5) %>%
  group_by(n_pop, sample_size)%>%
  summarize("meangxepower" = mean(gxepower),
            "meancovpower" = mean(covpower))

hundiehighsum = hundiepow %>%
  filter(std_dev == 1) %>%
  group_by(n_pop, sample_size)%>%
  summarize("meangxepower" = mean(gxepower),
            "meancovpower" = mean(covpower))

hundielowgxe = ggplot(hundielowsum,aes(x = factor(sample_size), y = factor(n_pop), fill = meangxepower)) + 
  geom_tile() + scale_fill_gradient(low="grey", high="dodgerblue4") +
  xlab("Sample Size") + ylab("Number of Populations")+
  labs(fill = "Power")+
  theme_classic() + ggtitle("GxE: Standard Deviation = 0.5")
hundielowcov = ggplot(hundielowsum,aes(x = factor(sample_size), y = factor(n_pop), fill = meancovpower)) + 
  geom_tile() + scale_fill_gradient(low="grey", high="darkgreen") +
  xlab("Sample Size") + ylab("Number of Populations")+
  labs(fill = "Power")+
  theme_classic() + ggtitle("CovGE: Standard Deviation = 0.5")
hundiehighgxe = ggplot(hundiehighsum,aes(x = factor(sample_size), y = factor(n_pop), fill = meangxepower)) + 
  geom_tile() + scale_fill_gradient(low="grey", high="dodgerblue4") +
  xlab("Sample Size") + ylab("Number of Populations")+
  labs(fill = "Power")+
  theme_classic() + ggtitle("GxE: Standard Deviation = 1")
hundiehighcov = ggplot(hundiehighsum,aes(x = factor(sample_size), y = factor(n_pop), fill = meancovpower)) + 
  geom_tile() + scale_fill_gradient(low="grey", high="darkgreen") +
  xlab("Sample Size") + ylab("Number of Populations")+
  labs(fill = "Power")+
  theme_classic() + ggtitle("CovGE: Standard Deviation = 1")
grid.arrange(hundielowgxe,hundiehighgxe,hundielowcov,hundiehighcov,ncol = 2)


fittylowsum = fittypow %>%
  filter(std_dev == 0.5) %>%
  group_by(n_pop, sample_size)%>%
  summarize("meangxepower" = mean(gxepower),
            "meancovpower" = mean(covpower))

fittyhighsum = fittypow %>%
  filter(std_dev == 1) %>%
  group_by(n_pop, sample_size)%>%
  summarize("meangxepower" = mean(gxepower),
            "meancovpower" = mean(covpower))

fittylowgxe = ggplot(fittylowsum,aes(x = factor(sample_size), y = factor(n_pop), fill = meangxepower)) + 
  geom_tile() + scale_fill_gradient(low="grey", high="dodgerblue4") +
  xlab("Sample Size") + ylab("Number of Populations")+
  labs(fill = "Power")+
  theme_classic() + ggtitle("GxE: Standard Deviation = 0.5")
fittylowcov = ggplot(fittylowsum,aes(x = factor(sample_size), y = factor(n_pop), fill = meancovpower)) + 
  geom_tile() + scale_fill_gradient(low="grey", high="darkgreen") +
  xlab("Sample Size") + ylab("Number of Populations")+
  labs(fill = "Power")+
  theme_classic() + ggtitle("CovGE: Standard Deviation = 0.5")
fittyhighgxe = ggplot(fittyhighsum,aes(x = factor(sample_size), y = factor(n_pop), fill = meangxepower)) + 
  geom_tile() + scale_fill_gradient(low="grey", high="dodgerblue4") +
  xlab("Sample Size") + ylab("Number of Populations")+
  labs(fill = "Power")+
  theme_classic() + ggtitle("GxE: Standard Deviation = 1")
fittyhighcov = ggplot(fittyhighsum,aes(x = factor(sample_size), y = factor(n_pop), fill = meancovpower)) + 
  geom_tile() + scale_fill_gradient(low="grey", high="darkgreen") +
  xlab("Sample Size") + ylab("Number of Populations")+
  labs(fill = "Power")+
  theme_classic() + ggtitle("CovGE: Standard Deviation = 1")
grid.arrange(fittylowgxe,fittyhighgxe,fittylowcov,fittyhighcov,ncol = 2)

