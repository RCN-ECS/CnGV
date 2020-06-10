#Another Day, Another Problem

The permutation code is not giving me accurate p-values for the GxE means values. 

I have tried everything I can think of (i've been debugging for 3 days straight) and I cannot find the problem. 

FOR EXAMPLE. Here is a 5 populations scenario that SHOULD have significant GxE across the board: 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/609_meanProbplot.png)

however, when I run my code, here are my GxE results: 


```{r}
GxE_emm_loop              0.65 # New way of estimating using for-loop
true_GxE_emm_loop         0.66
true_GxE_emm_loop_pvalue  0.00
GxE_emm_loop_lwrCI        0.64
GxE_emm_loop_uprCI        0.66
GxE_emm_loop_pvalue       0.00

true_omega_GxE            0.71 # Omega^2
true_omega_GxE_pvalue     0.00
GxE_omega_estimate        0.69
GxE_omega_lwrCI           0.68
GxE_omega_uprCI           0.71
GxE_omega_pvalue          0.00

true_eta_GxE              0.71 # Eta^2
true_eta_GxE_pvalue       0.00
GxE_eta                   0.70
GxE_eta_lwrCI             0.69
GxE_eta_uprCI             0.71
GxE_eta_pvalue            0.00

True_resid_variation      0.71 # Residual Variation
Resid_variation           0.70
GxE_residVar_lwrCI        0.69
GxE_residVar_uprCI        0.72

true_GxE_means            0.65 # Means data - HERE IS THE PROBLEM
true_GxE_means_pvalue     0.50
GxE_means                 0.64
GxE_means_lwrCI           0.63
GxE_means_uprCI           0.65
GxE_means_pvalue          0.50
```
As seen, the means are not showing up as significant even though they should. 

My entire permutation code is posted below. Full code is [here](https://github.com/RCN-ECS/CnGV/blob/master/src/Power_simulation_ClusterFormat_bonuscode.R) I cannot figure this out, and I really hope its not something really dumb. 

Starting Parameters:
```{r}
replicate   "5"            
delta_env   "1"            
delta_gen   "0.01"         
sample_size "5"            
n_pop       "5"            
std_dev     "0.5"          
interaction "5"      
```
```{permutation code}
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

# Gmeans - Permutation
G_matrix_perm = data.frame()
for(r in 1:length(unique(emm_GxE_perm$gen_factor))){
  g_perm <- filter(emm_GxE_perm, gen_factor == unique(emm_GxE_perm$gen_factor)[r])
  g_mean_perm <- sum(g_perm[,3])/length(unique(emm_GxE_perm$gen_factor))
  tempdat = data.frame("G_means" = g_mean_perm,
                       "gen_factor" = unique(g_perm$gen_factor))
  G_matrix_perm = rbind(G_matrix_perm,tempdat)
}

# Emeans - Permutation
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

#if(length(unique(perm_means$gen_factor))>2{
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

```
