# Power Analysis

source("~/Documents/GitHub/CnGV/src/data_generation_function.R") # Generate data (either cat. or cont.)
source("~/Documents/GitHub/CnGV/src/sim_means_se.R") # Generate means and SE from raw sim. data (either cat. or cont.)

source("~/Documents/GitHub/CnGV/src/slope_generation_cont.R") # Bootstrapped slopes from means and std. error - cont.
source("~/Documents/GitHub/CnGV/src/slope_generation_cat.R") # Bootstrapped slopes from means and std. error - cat.

source("~/Documents/GitHub/CnGV/src/meansSE_boot_cont.R") # Gmeans, and Emeans from means and std. error - continuous
source("~/Documents/GitHub/CnGV/src/meansSE_boot_cat.R") # Gmeans, and Emeans from means and std. error - categorical

source("~/Documents/GitHub/CnGV/src/CovMatrix_sim_cont.R") # Covariance estimates from simulated raw data - continuous
source("~/Documents/GitHub/CnGV/src/Cov_matrix_sim_cat.R") # Covariance estimates from simulated raw data - categorical


# Starting parameters
Diff_means_cat <- list(
  "data_type" = c("categorical"), 
  "intercept_G1" = seq(from = -5, to = 5, by = 2),
  "slope_G1" = seq(from = -1, to = 1, by = 0.5),
  "intercept_G2" = seq(from = -5, to = 5, by = 2),
  "slope_G2" = seq(from = -1, to = 1, by = 0.5), 
  "sd" = 0.05, #seq(from = 0, to = 1, by = 0.5),
  "sample_size" = c(5)) 

cat_raw <- data.frame(data_generation(Diff_means_cat)) # raw
cat_raw_test <- Cov_matrix_sim_cat(cat_raw) # Need some error or will divide by zero and 'splode.
(cat_raw_test[[3]])

#Re run but filter only where covariance = 0
options(scipen = 999)
df = cat_raw_test[[3]]
df %>% filter(dplyr::between(df$age, 5, 25))
filtdat = df %>% filter(between(df$Covariance_est,0,0.001))
df$cov_type_actual = NULL

for(i in 1:nrow(df)){
if(df$Covariance_est[i] < -0.01){df$cov_type_actual[i] = "CnGV"
  }else if(df$Covariance_est[i] > 0.01){df$cov_type_actual[i] = "CoGV"
  }else{df$cov_type_actual[i] = "pure_GxE"}
}

col3 = c("CnGV" = "blue","CoGV" = "green", "pure_GxE" = "black")
print(ggplot(df, aes(x = Covariance_est, y= abs(phendiff), colour = cov_type_actual))) + geom_point()+ 
  theme_classic() + scale_colour_manual(values = col3)+ #geom_smooth(colour = "black")+
  #geom_hline(yintercept =0) + 
  #xlim(0,1.0)+ylim(0,3)+
  #geom_ribbon(aes(ymin = 0, ymax = predict(loess(abs(phendiff) ~ GxE_emm))),alpha = .3)+
  ylab("Phenotypic Difference between G1E1 and G2E2") + xlab("Covariance Estimate")
  #annotate(geom="text", x=0.25, y=0.8, label="CnGV",color="black")+
  #annotate(geom="text", x=0.25, y=-1, label="CnGV",color="black") +
  #annotate(geom="text", x=0.75, y=2.5, label="CoGV",color="black")
  #annotate(geom="text", x=0.8, y=-3, label="CoGV",color="black")

