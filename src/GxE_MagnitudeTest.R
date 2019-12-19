source("~/Documents/GitHub/CnGV/src/data_generation_function.R") # Generate data (either cat. or cont.)
source("~/Documents/GitHub/CnGV/src/sim_means_se.R") # Generate means and SE from raw sim. data (either cat. or cont.)
source("~/Documents/GitHub/CnGV/src/meansSE_boot_cont.R") # Gmeans, and Emeans from means and std. error - continuous
source("~/Documents/GitHub/CnGV/src/meansSE_boot_cat.R") # Gmeans, and Emeans from means and std. error - categorical
source("~/Documents/GitHub/CnGV/src/CovMatrix_sim_cont.R") # Covariance estimates from simulated raw data - continuous
source("~/Documents/GitHub/CnGV/src/Cov_matrix_sim_cat.R") # Covariance estimates from simulated raw data - categorical

# Starting parameters
data <- list(
  "data_type" = c("continuous"), #continuous or categorical
  "type" = NULL, #c("cogv","cngv","pure_GxE"),
  "intercept_G1" = -1,
  "slope_G1" = c(1), #seq(from = -1, to = 1, by = 0.5),
  "intercept_G2" = 1,#c(0), #seq(from = -5, to = 5, by = 2),
  "slope_G2" = c(-1), #seq(from = -1, to = 1, by = 0.1),
  "sd" = 0, #seq(from = 0, to = 1, by = 0.5),
  "sample_size" = c(5), #seq(from = 5, to = 10, by = 2),
  "env_num" = 2, #c(2,5,10), 
  "env" = NA,
  "G1_env" = NA, 
  "G2_env" = NA,
  "true_covGE" = NA, 
  "is.GxE" = NA, 
  "slope_diff" = NA) 

repeatercont <- function(data){
  
  iterations <- 1 
  
  if(data$data_type == "continuous"){
    results = data.frame()
    for(i in 1:iterations){
      cont_raw <- data.frame(data_generation(data))
      cont_raw_test <- Cov_matrix_sim_cont(cont_raw)
      results <- rbind(results,cont_raw_test[[3]])
  }
  }else{
    results = data.frame()
    for(i in 1:iterations){
      cat_raw <- data.frame(data_generation(data))
      cat_raw_test <- Cov_matrix_sim_cat(cat_raw)
      results <- rbind(results,cat_raw_test[[3]])
    }}
    
    GxE_eta_CI =  quantile(results$GxE_eta, probs=c(0.025, 0.975), type=1)
    GxE_omega_CI = quantile(results$GxE_omega, probs=c(0.025, 0.975), type=1)
    GxE_lot_CI = quantile(results$GxE_lot, probs=c(0.025, 0.975), type=1)
    GxE_emm_CI = quantile(results$GxE_emm, probs=c(0.025, 0.975), type=1)
    
    resultframe = data.frame(
        "GxE_eta_mean" = mean(results$GxE_eta),
        "GxE_eta_lwrCI" = GxE_eta_CI[[1]],
        "GxE_eta_uprCI"  = GxE_eta_CI[[2]],
        "GxE_omega_mean" = mean(results$GxE_omega),
        "GxE_omega_lwrCI" = GxE_omega_CI[[1]],
        "GxE_omega_uprCI" = GxE_omega_CI[[2]],
        "GxE_lot_mean" = mean(results$GxE_lot),
        "GxE_lot_lwrCI" = GxE_lot_CI[[1]],
        "GxE_lot_uprCI" = GxE_lot_CI[[2]],
        "GxE_emm_mean" = mean(results$GxE_emm),
        "GxE_emm_lwrCI" = GxE_emm_CI[[1]],
        "GxE_emm_uprCI" = GxE_emm_CI[[2]])
  return(resultframe)
}

print(ggplot(cont_raw, aes(x = env, y = phen, group = gen,colour = gen))) + geom_smooth(method = "gam") + theme_classic()

hooloo = repeatercont(data)
hooloo

## Katie's original method: 

# Means of genotypes in each environmnent
G11 <- -0.25 #Env 1, G1
G12 <- 0.25   #Env 2, G1
G21 <- 0.25   #Env 1, G2
G22 <- -0.25 #Env 2, G2

calc_interaction <- function(G11, G12, G21, G22){
  
  plot(c(1,2), c(G11, G12), type="l", col="red", ylim=c(-2,2))
  points(c(1,2), c(G21, G22), type="l", col="blue", ylim=c(-2,2))
  
  
  (overall_mean <- mean(c(G11,G12,G21,G22)))
  G11m <- G11 - overall_mean
  G12m <- G12 - overall_mean
  G21m <- G21 - overall_mean
  G22m <- G22 - overall_mean
  (overall_mean <- mean(c(G11m,G12m,G21m,G22m)))
  print(overall_mean)
  
  # Marginal means
  (G1_mean <- mean(c(G11m,G12m)))
  (G2_mean <- mean(c(G22m,G21m)))
  (E1_mean <- mean(c(G11m,G21m)))
  (E2_mean <- mean(c(G12m,G22m)))
  
  ## Interaction effect for the i_th level of factor G (G2)
  ## and jth level of factor E (E2)
  print(c("Interaction effect for G22:",  
          overall_mean - G2_mean - E2_mean + G22))
  
  print(c("Interaction effect for G11:", 
          overall_mean - G1_mean - E1_mean + G11))
  
  print(c("Interaction effect for G12:",
          overall_mean - G1_mean - E2_mean + G12))
  
  print(c("Interaction effect for G21:", 
          overall_mean - G2_mean - E1_mean + G21))
}

calc_interaction(0.5, 1, 0, 0.5)
calc_interaction(0.5, 1, 0, 0.6)
t1 = calc_interaction(-1, 1, 1, -1)
t2 = calc_interaction(-0.25, 0.25, 0.25, -0.25)

