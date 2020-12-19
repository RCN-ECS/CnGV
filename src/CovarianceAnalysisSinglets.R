#######################################################################################
##                      Co/counter gradient Analysis                                 ##
##  Authors: Molly Albecker, Thais Bittar, Geoff Trussell, Katie Lotterhos           ##
#######################################################################################

# Thermal Optima Goofaround: 
gen_factor = c("G_1","G_2","G_3")
exp_env_factor = c("E_1","E_2","E_3","E_4","E_5","E_6")
nat_env_factor = c("E_3","E_4","E_5")
G1phen = c(1,3,5,4.5,4,3.5)
G2phen = c(0,2,4,5,4.5,4)
G3phen = c(0.5, 2, 3, 4,5, 4.5 )
phen_mean_SE = abs(rnorm(18, mean = 0, sd = 0.1))
thermoptcngv = data.frame("phen_data" = c(G1phen, G2phen,G3phen), 
                      "phen_mean_SE" = phen_mean_SE,
                      "gen_factor" = rep(gen_factor, each = 6), 
                      "exp_env_factor" = rep(exp_env_factor,3), 
                      "nat_env_factor" = rep(nat_env_factor, each = 6))
data.df= termoptcngv

G1phen = c(0,1,3,2.5,2,1.5)
G2phen = c(0,2,3,4,3,2)
G3phen = c(0,2,3,4,5,2)
phen_mean_SE = abs(rnorm(18, mean = 0, sd = 0.1))
thermoptcogv = data.frame("phen_data" = c(G1phen, G2phen,G3phen), 
                          "phen_mean_SE" = phen_mean_SE,
                          "gen_factor" = rep(gen_factor, each = 6), 
                          "exp_env_factor" = rep(exp_env_factor,3), 
                          "nat_env_factor" = rep(nat_env_factor, each = 6))
data.df= thermoptcogv


# Load packages
library("emmeans")
library("lme4")
library("dplyr")
library("gridExtra")

# 1. Load Data
setwd("~/Documents/GitHub/CnGV/CnGV/data/")
data.df1 <- read.csv("meta_df.csv")
data.df1$Data_file_name <- as.character(data.df1$Data_file_name)
data.df = filter(data.df1, data.df1$Data_file_name == "1064_Berggren_hatching_success")

# 2. Load Functions (CovarianceDataFunctions.R)
setwd("~/Documents/GitHub/CnGV/CnGV/src/")
source("CovarianceDataFunctions.R")

# 3. Format Data 
# Need: "gen_factor" column -as factor (each genotype/population = "G_1", "G_2", etc.)
# Need: "exp_env_factor" column - as factor (each experimental treatment = "E_1", "E_2", etc.)
# Need: "nat_env_factor" column - as factor (each genotype's native environment = "E_1", "E_2", etc.)
# Need: "phen_data" column - simply the column with focal phenotypic data
data.df$gen_factor = factor(data.df$gen_factor)
data.df$exp_env_factor = factor(data.df$exp_env_factor)
data.df$nat_env_factor = factor(data.df$nat_env_factor)

# 4. Establish starting conditions
n_boot <- 99
  
# 5. Is the data "raw" or "means" format? 
data_type = "means"

# 6. Standardize by standard deviation of group means
if(data_type == "raw"){ 
  data.df$group = paste(data.df$gen_factor,data.df$exp_env_factor,sep = "-")
  data.df$phen_corrected = (data.df$phen_data - mean(data.df$phen_data))/sd(tapply(data.df$phen_data, data.df$group, mean))
  }else{
  data.df$avg_phen_corrected = (data.df$phen_data - mean(data.df$phen_data))/sd(data.df$phen_data)
  }
  
# 7. Sanity Plots to form expectations
 (a = ggplot(data.df,aes(x=exp_env_factor,y=avg_phen_corrected, group = gen_factor, colour=nat_env_factor))+
   geom_point()+geom_line()+theme_classic() + ggtitle("Standardized Data") + labs(colour = "Native Environment")+
     geom_errorbar(aes(ymin = avg_phen_corrected-phen_mean_SE, ymax = avg_phen_corrected+phen_mean_SE),width = 0.1))
 b = ggplot(data.df,aes(x=exp_env_factor,y=phen_data, group = gen_factor, colour=nat_env_factor))+
   geom_point()+geom_line()+theme_classic() + ggtitle("Raw Data")
 grid.arrange(a,b)
 
# 8. Run the analysis
output = amarillo_armadillo(data.df, n_boot, data_type)
output
