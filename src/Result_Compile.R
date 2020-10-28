################################################
###             Package Results Code         ###
################################################

## Load packages
library(readr)
library(plyr)
library(gridExtra)

## Compile Data
setwd("/scratch/albecker/Power_analysis/power_output/")
Sim1 <- list.files(pattern= "*Results_") 
sim2 <- plyr::ldply(Sim1, read_csv)
write.csv(sim2, "/scratch/albecker/Power_analysis/Power_output_results.csv")

Var1 <- list.files(pattern= "*variance_")
var2 <- plyr::ldply(Var1, read_csv)
write.csv(var2, "/scratch/albecker/Power_analysis/Variance_output_results.csv")

setwd("/scratch/albecker/Power_analysis/PL_output/")
PL <- list.files(pattern= "*PL_") 
PL1 = plyr::ldply(PL, read_csv)
write.csv(PL1, "/scratch/albecker/Power_analysis/PL_output_results.csv")
 
setwd("/scratch/albecker/Power_analysis/phenotype_output/")
Phen <- list.files(pattern= "*Phenotype_") 
Phen1 = plyr::ldply(Phen, read_csv)
write.csv(Phen1, "/scratch/albecker/Power_analysis/phenotype_output_results.csv")

setwd("/scratch/albecker/Power_analysis/permutation_output/")
perm <- list.files(pattern= "*Phenotype_") 
perm1 = plyr::ldply(perm, read_csv)
write.csv(perm1, "/scratch/albecker/Power_analysis/permutation_output_results.csv")

setwd("/scratch/albecker/Power_analysis/GEmeans_output/")
cov <- list.files(pattern= "*covmatrix_") 
cov1 = plyr::ldply(cov, read_csv)
write.csv(cov1, "/scratch/albecker/Power_analysis/covmatrix_output_results.csv")

setwd("/scratch/albecker/Power_analysis/Anova_output/")
ssq <- list.files(pattern= "*model_SSQ_") 
ssq1 = plyr::ldply(ssq, read_csv)
write.csv(ssq1, "/scratch/albecker/Power_analysis/model_SSQ_output_results.csv")

setwd("/scratch/albecker/Power_analysis/Anova_output/")
coef <- list.files(pattern= "*anova_coef_") 
coef1 = plyr::ldply(coef, read_csv)
write.csv(coef1, "/scratch/albecker/Power_analysis/anova_coef_output_results.csv")