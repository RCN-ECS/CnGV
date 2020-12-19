#################################################################################################
##              Functions for co/counter gradient meta-analysis                                ##
##        Authors: Molly Albecker, Thais Bittar, Geoff Trussell, Katie Lotterhos               ##
#################################################################################################

# Load packages
library("emmeans")
library("lme4")
library("dplyr")

# Load Compiled datasets
setwd("~/Documents/GitHub/CnGV/CnGV/data/")
info.df <- read.csv("Extraction_Initialize.csv")
data.df1 <- read.csv("meta_df.csv")
output <- data.frame()

# Data wrangle
data.df <- data.df1[-which(is.na(data.df1$Data_file_name)),] # Rid of empty cells
data.df$Data_file_name <- as.character(data.df$Data_file_name)
info.df$Data_file_name <- as.character(info.df$Data_file_name)

# Test data
# testies1 = c("1064_Berggren_hatching_success","126_Debecker_GST_content","1904_Firman_rapid_sperm","3095_Smith_bulk_density","936_Hamann_inflorescence","945_Lucey_living_individuals_pluseggs")
# testies2 = c("652_Urban_outside_leaf_refugia","4493_Smith_maturation_size","303_Brancalion_leaf_mass", "543_Pepino_mouth_width","45_Verheyen_SOD_activity","1129_Cattano_length_hatchiling")
# data.df. = filter(data.df, data.df$Data_file_name %in% testies2)
# data.df = data.df.

# Load functions
 setwd("~/Documents/GitHub/CnGV/CnGV/src/")
source("CovarianceDataFunctions.R")

# Run analysis
for(x in 1:length(unique(data.df$Data_file_name))){
  
  # Pull out data for each study 
  temp = filter(data.df, Data_file_name == unique(data.df$Data_file_name)[x])
  ref = filter(info.df, Data_file_name == unique(temp$Data_file_name))
  
  # Format data 
  temp$gen_factor = factor(temp$gen_factor)
  temp$exp_env_factor = factor(temp$exp_env_factor)
  
  # Establish starting conditions
  skip_to_next <- FALSE
  n_boot <- 999
  
  # Rename native environments 
  temp$nat_env_factor = gsub("N_", "E_", temp$nat_env_factor)
  
  # Is the data "raw" or "means" format? 
  data_type = NA
  if(unique(temp$Raw_data_available) == "raw data"){data_type = "raw"}else{data_type = "means"}
  
  # Standardize by standard deviation of group means
  if(data_type == "raw"){ 
    temp$group = paste(temp$gen_factor,temp$exp_env_factor,sep = "-")
    temp$phen_corrected = (temp$phen_data - mean(temp$phen_data))/sd(tapply(temp$phen_data, temp$group, mean))
  }else{
    temp$avg_phen_corrected = (temp$phen_data - mean(temp$phen_data))/sd(temp$phen_data)
  }
  
  tryCatch(
    {temp2 = amarillo_armadillo(temp, n_boot, data_type)
  
    # Compile Results 
    temp_out <- data.frame("Index" = x, 
                           "Data_file_name" = ref$Data_file_name,
                           "First.Author" = ref$First_Author,
                           "Data.type" = data_type,
                           "Phylum" = ref$phylum_division,
                           "Genus" = ref$genus,
                           "Species" = ref$species,
                           "Phenotype" = ref$phenotype,
                           "Phenotype.unit" = ref$phenotype_unit,
                           "Experimental_comparison" = ref$experimental_comparison,
                           "Mention.CGV?" = ref$Mention.CGV,
                           "Sig.GxE?" = ref$GxE_sig,
                           "Covariance" = temp2$Covariance.Estimate,
                           "Covariance_LCI" = temp2$Covariance.Lower.CI,
                           "Covariance_UCI" = temp2$Covariance.Upper.CI,
                           "Covariance_Pvalue" = temp2$Covariance.p.value,
                           "GxE_Estimate" = temp2$GxE.Estimate,
                           "GxE_LCI" = temp2$GxE.Lower.CI,
                           "GxE_UCI" = temp2$GxE.Upper.CI,
                           "GxE_Pvalue" = temp2$GxE.p.value,
                           "Omega2" = temp2$Omega2,
                           "Omega2_LCI" = temp2$Omega2.Lower.CI,
                           "Omega2_UCI" = temp2$Omega2.Upper.CI,
                           "Omega2_Pvalue" = temp2$Omega2.p.value)}, error = function(e) {skip_to_next <<- TRUE})
  
  if(skip_to_next) {
    exit_df = data.frame("index" = unique(temp$Data_file_name))
    #write.csv(exit_df, paste0("~/Desktop/fail_",x,".csv")) 
    write.csv(exit_df, paste0("/scratch/albecker/Power_analysis/fails/fail_",x,".csv")) 
    next }else{
    output = rbind(output,temp_out)
    }
}

write.csv(output, "Meta_analysis_results.csv")

##################################  
##    Result Visualization      ##
##################################  
library(ggplot2)
library(viridis)

res <- read.csv("~/Desktop/Meta_analysis_results.csv")
dim(res)
unique(res$Phylum)
phyShape <- c("Chordata"=15,"Cnidaria"=16,"Arthropoda"=17,"Dinoflagellata"=18,"Bryophyta"=0,"Tracheophyta"=1,"Mollusca"=2,"Annelida"=5)

# CovGE vs. GxE
res$color = NULL
res = res %>%
   mutate(color = ifelse(Covariance_Pvalue <= 0.025 & GxE_Pvalue <= 0.05, "3",
                         ifelse(GxE_Pvalue <= 0.05, "2", #"#481F70FF",\
                                ifelse(Covariance_Pvalue <= 0.025, "1","4")))) #"#20A486FF",
                                      # ifelse(is.na(Covariance_Pvalue) == TRUE & is.na(GxE_Pvalue) == TRUE, NA, "4"))))) #"#FDE725FF", "white")))))

ggplot(filter(res,is.na(color)), aes(x = Covariance, y = GxE_Estimate, fill = color))+
  labs(fill = "Significance")+
  geom_vline(aes(xintercept = 0))+
  geom_point(shape = 21, size = 4) + 
  scale_fill_manual(values = c("1" = "#20A486FF","2" = "#481F70FF", "3" = "#FDE725FF", "4" = "white"),
    labels = c("1" = "Covariance Significant", "2" = "GxE Significant", "3" = "Both Significant","4"="None Signficant"))+
  xlab(expression("Cov"["GE"]))+
  ylab(expression(bar(Delta)*""["GxE"]))+
  theme(legend.position="bottom")  +
  theme_classic(base_family = "Times", base_size = 20)+ facet_wrap(~factor(Phylum),ncol=2)



