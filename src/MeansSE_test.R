source("~/Documents/GitHub/CnGV/src/data_generation_function.R") # Generate data (either cat. or cont.)
source("~/Documents/GitHub/CnGV/src/sim_means_se.R") # Generate means and SE from raw sim. data (either cat. or cont.)

source("~/Documents/GitHub/CnGV/src/slope_generation_cont.R") # Bootstrapped slopes from means and std. error - cont.
source("~/Documents/GitHub/CnGV/src/slope_generation_cat.R") # Bootstrapped slopes from means and std. error - cat.

source("~/Documents/GitHub/CnGV/src/meansSE_boot_cont.R") # Gmeans, and Emeans from means and std. error - continuous
source("~/Documents/GitHub/CnGV/src/meansSE_boot_cat.R") # Gmeans, and Emeans from means and std. error - categorical

source("~/Documents/GitHub/CnGV/src/CovMatrix_sim_cont.R") # Covariance estimates from simulated raw data - continuous
source("~/Documents/GitHub/CnGV/src/Cov_matrix_sim_cat.R") # Covariance estimates from simulated raw data - categorical

# Continuous Starting parameters
Diff_means_cont <- list(
  "data_type" = c("continuous"), #continuous or categorical
  #"type" = NULL, #c("cogv","cngv","pure_GxE"),
  "intercept_G1" = c(4,5,6),
  "slope_G1" = 1, #seq(from = -1, to = 1, by = 0.5),
  "intercept_G2" = c(1,2,3),#c(0), #seq(from = -5, to = 5, by = 2),
  "slope_G2" = 1, #seq(from = -1, to = 1, by = 0.1),
  "sd" = 0.5, #seq(from = 0, to = 1, by = 0.5),
  "sample_size" = c(5)) #seq(from = 5, to = 10, by = 2),
  #"env_num" = 2, #c(2,5,10), 
  #"env" = NA,
  #"G1_env" = NA, 
  #"G2_env" = NA,
  #"true_covGE" = NA, 
  #"is.GxE" = NA, 
  #"slope_diff" = NA) 

# Create continuous data
cont_raw <- data_generation(Diff_means_cont) # raw
cont_mean <- sim_means_se(cont_raw) # means (This isnt set up to handle multiple indexes yet)

# Plot
print(ggplot(cont_raw, aes(x = env, y = phen, group = gen,colour = gen))) + geom_smooth(method = "gam") + theme_classic()

# Continous covariance and slope estimation
cont_raw_test <- Cov_matrix_sim_cont(cont_raw)
(cont_raw_test[[3]]) 
cont_mean_test <- meansSE_boot_cont(cont_mean, 100) # not set up for multiple indexes
(cont_mean_test[[3]]) 


# Categorical Starting parameters
Diff_means_cat <- list(
  "data_type" = c("categorical"), #continuous or categorical
  "intercept_G1" = seq(from = -5, to = 5, by = 2),
  "slope_G1" = seq(from = -1, to = 1, by = 0.5), #seq(from = -1, to = 1, by = 0.5),
  "intercept_G2" = seq(from = -5, to = 5, by = 2),#c(0), #seq(from = -5, to = 5, by = 2),
  "slope_G2" = seq(from = -1, to = 1, by = 0.5), #seq(from = -1, to = 1, by = 0.1),
  "sd" = 0.5, #seq(from = 0, to = 1, by = 0.5),
  "sample_size" = c(5)) #seq(from = 5, to = 10, by = 2),

# Generate categorical data
cat_raw <- data.frame(data_generation(Diff_means_cat)) # raw
cat_mean <- sim_means_se(cat_raw) # means (This isnt set up to handle multiple indexes yet)

# Categorical - covariance and slopes
cat_raw_test <- Cov_matrix_sim_cat(cat_raw)
(cont_raw_test[[3]])

cat_mean_test <- meansSE_boot_cat(cat_mean,100) # not set up for multiple indexes
(cat_mean_test[[3]])

#Big Poppa Test:

poppa_fun <- function(cont_raw,cat_raw){ #Continuous and categorical simulated raw data are inputs
  
  # Output
  cont_outdat = data.frame()
  cat_outdat = data.frame()

  for(i in 1:length(unique(cont_raw$index))){
  cont_level <- cont_raw$index[i]
  cont_dat <- filter(cont_raw, index == cont_level)
  
  cont_outdat. <- Cov_matrix_sim_cont(cont_raw)
  cont_outdat <- rbind(cont_outdat.[[3]],cont_outdat)
  }
  
  for(j in 1:length(unique(cat_raw$index))){
    cat_level <- cat_raw$index[i]
    cat_dat <- filter(cat_raw, index == cat_level)
    
    cat_outdat. <- Cov_matrix_sim_cat(cat_raw)
    cat_outdat <- rbind(cat_outdat.[[3]],cat_outdat)
  }
return(list(cat_outdat,cont_outdat))
  }

# Get data for plotting:
plot_dat <- poppa_fun(cont_raw,cat_raw)
cont_plot <- plot_dat[[2]]
cat_plot <- plot_dat[[1]]

#Continuous Plots: 

# EMM  
print(ggplot(cont_plot,aes(x = Covariance_est, y = GxE_emm))) + geom_jitter() + theme_classic()

# Lotterhos method 
print(ggplot(cont_plot,aes(x = Covariance_est, y = GxE_lot))) + geom_jitter() + theme_classic()

# Omega 
print(ggplot(cont_plot,aes(x = Covariance_est, y = GxE_omega))) + geom_jitter() + theme_classic()

# Categorical Plots: 
colorscheme = c( "check"="grey","cngv"="blue","cogv"= "red","pure_GxE"="black")
# EMM 
print(ggplot(cat_raw_test[[3]],aes(x = Covariance_est, y = GxE_emm))) + 
  geom_jitter() + theme_classic() + ylab("GxE magnitude - Emmeans") + xlab("Covariance Estimate") #+ scale_color_manual(values = colorscheme)

# Lotterhos method
print(ggplot(cat_raw_test[[3]],aes(x = Covariance_est, y = GxE_lot))) + geom_jitter() + theme_classic() + ylab("GxE magnitude - Emmeans-Lotterhos") + xlab("Covariance Estimate")

# Omega
print(ggplot(cat_raw_test[[3]],aes(x = Covariance_est, y = GxE_omega))) + geom_jitter() + theme_classic()+ ylab("GxE magnitude - Omega^2") + xlab("Covariance Estimate")


