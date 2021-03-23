#######################################################################################
##                      Co/counter gradient Analysis                                 ##
##  Authors: Molly Albecker, Thais Bittar, Geoff Trussell, Katie Lotterhos           ##
#######################################################################################

# Load packages
library("emmeans")
library("lme4")
library("dplyr")
library("gridExtra")

# 1. Load Data
setwd("~/Documents/GitHub/CnGV/CnGV/data/")
data.df1 <- read.csv("meta_df.csv")
data.df1$Data_file_name <- as.character(data.df1$Data_file_name)
data.df = filter(data.df1, data.df1$Data_file_name == "1632_Wrangle_shell_strength")

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
n_boot <- 999
  
# 5. Is the data "raw" or "means" format? 
data_type = unique(data.df$data_type)

# 6. Standardize by standard deviation of group means
if(data_type == "raw"){ 
  data.df$group = paste(data.df$gen_factor,data.df$exp_env_factor,sep = "-")
  data.df$phen_corrected = (data.df$phen_data - mean(data.df$phen_data))/sd(tapply(data.df$phen_data, data.df$group, mean))
  }else{
  data.df$avg_phen_corrected = (data.df$phen_data - mean(data.df$phen_data))/sd(data.df$phen_data)
  }
  
# 7. Sanity Plots to form expectations
 ggplot(data.df,aes(x=exp_env_factor,y=avg_phen_corrected, group = gen_factor, colour=nat_env_factor))+
   geom_point()+geom_line()+theme_classic() + ggtitle("Standardized Data") + labs(colour = "Native Environment")+
     geom_errorbar(aes(ymin = avg_phen_corrected-phen_mean_SE, ymax = avg_phen_corrected+phen_mean_SE),width = 0.1)
 ggplot(data.df,aes(x=exp_env_factor,y=phen_data, group = gen_factor, colour=nat_env_factor))+
   geom_point()+geom_smooth()+theme_classic() + ggtitle("Raw Data")

# 8. Run the analysis
output = amarillo_armadillo(data.df, n_boot, data_type)
output
