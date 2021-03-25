##########################################################
###             Simulation Visualization Code          ###
##########################################################

## Load packages
library(ggplot2)
library(grid)
library(readr)
library(tidyverse)
library(gridExtra)
library(ggthemes)
library(viridis)

# Load Data compiled on cluster
setwd("~/Documents/GitHub/CnGV/CnGV/results/Sim_3.10.21/")
start_params = read.csv("df_rerun.csv")                      # Input Parameters
start_df = read.csv("Power_output_results.csv")              # Results of simulations 

phen_data = read.csv("~/Desktop/phenotype_output_results.csv")   # Data for phenotype plots

# Split up into two experimental designs
dat_csv = start_df %>% filter(env_scenario == 1) %>% droplevels() # Reciprocal Transplant
dat_dub = start_df %>% filter(env_scenario == 2) %>% droplevels() # Common Garden

####### Check Average time for longest sims ##########
(timecheck = start_df %>%
  #filter(n_pop == max(n_pop)) %>%
  #filter(std_dev == max(std_dev)) %>%
  #filter(sample_size == max(sample_size)) %>%
  summarize(average_time = mean(Sim_time)))

####### Check 0s for FPR and FNR ##########
sizecheck = dat_csv %>%
   #filter(true_cov == 0) 
   filter(true_GxE_emm >= 0.75)
   #filter(abs(true_cov) >=0.75) 
dim(sizecheck) # Should be around 1000 (or 100 if 1 replicate)

####### Check # for each type ##########
(sizecheck = start_df %>%
    filter(n_pop == min(n_pop)) %>%
    filter(sample_size == min(sample_size)) %>%
    #filter(n_pop == max(n_pop)) %>%
    #filter(sample_size == max(sample_size)) %>%
    summarize(size = n()))

####### Check for missing rows ##########
`%notin%` <- Negate(`%in%`)

h = start_df %>% filter(env_scenario == 2) %>% filter(n_pop == 2)  
dim(h)

(missing_rows = anti_join(start_df, start_params, by = "row"))
start_df = filter(start_df, row %notin% missing_rows$row)
#write.csv(start_df1[,-1], "~/Desktop/rerun.csv")



#####################################
##          Parameter Coverage      ##
######################################
dat_csv$ss_chr <- NA
dat_csv = dat_csv %>% mutate(ss_chr = ifelse(sample_size == 2, "2 Samples",
                                  ifelse(sample_size == 4, "4 Samples",
                                         ifelse(sample_size == 8, "8 Samples", "16 Samples"))))
dat_csv$np_chr <- NA
dat_csv = dat_csv %>% mutate(np_chr = ifelse(n_pop == 2, "2 Genotypes",
                                                         ifelse(n_pop == 4, "4 Genotypes","8 Genotypes")))
dat_csv$np_f = factor(dat_csv$np_chr, levels=c("2 Genotypes","4 Genotypes","8 Genotypes")) 
dat_csv$ss_f = factor(dat_csv$ss_chr, levels=c("2 Samples","4 Samples","8 Samples","16 Samples"))                                               

(hexy = ggplot(dat_csv, aes(x = true_cov, y = true_GxE_emm)) + 
   geom_hex()+ 
   ylab(expression(""*bar(Delta)*""["GxE"]*" of population"))+xlab(expression("Cov"["GE"]*" of population"))+
   ggtitle("A   Full Reciprocal Transplant") + facet_grid(ss_f~np_f) +
  #facet_grid(sample_size~n_pop)+
   theme_classic(base_family = "Times"))

dat_dub$ss_chr <- NA
dat_dub = dat_dub %>% mutate(ss_chr = ifelse(sample_size == 2, "2 Samples",
                                             ifelse(sample_size == 4, "4 Samples",
                                                    ifelse(sample_size == 8, "8 Samples", "16 Samples"))))
dat_dub$np_chr <- NA
dat_dub = dat_dub %>% mutate(np_chr = ifelse(n_pop == 2, "2 Genotypes",
                                             ifelse(n_pop == 4, "4 Genotypes",
                                                    ifelse(n_pop == 8, "8 Genotypes","16 Genotypes"))))
dat_dub$np_f = factor(dat_dub$np_chr, levels=c("2 Genotypes","4 Genotypes","8 Genotypes","16 Genotypes")) 
dat_dub$ss_f = factor(dat_dub$ss_chr, levels=c("2 Samples","4 Samples","8 Samples","16 Samples"))    
dat_dub = filter(dat_dub, n_pop != 2)
(hexy2 = ggplot(dat_dub, aes(x = true_cov, y = true_GxE_emm)) + 
    geom_hex()+ 
    ylab(expression(""*bar(Delta)*""["GxE"]*" of population"))+xlab(expression("Cov"["GE"]*" of population"))+
    ggtitle("B   Paired Common Garden") + facet_grid(ss_f~np_f) +
    theme_classic(base_family = "Times"))

grid.arrange(hexy,hexy2) 



#########################################################
##          Error Matrices  -- Recip Transplant        ##
#########################################################

## Covariance Permutation 
dat_csv$Covconfintperm = rep("NA",nrow(dat_csv))
for(i in 1:nrow(dat_csv)){
    if(dat_csv$true_cov[i] != 0 && dat_csv$covariance_pvalue[i] <= 0.05){dat_csv$Covconfintperm[i] = "True Positive"
    }else if(dat_csv$true_cov[i] == 0 & dat_csv$covariance_pvalue[i] <=  0.05){dat_csv$Covconfintperm[i] = "False Positive"
    }else if(dat_csv$true_cov[i]!= 0 & dat_csv$covariance_pvalue[i] >  0.05){dat_csv$Covconfintperm[i] = "False Negative"
    }else if(dat_csv$true_cov[i] == 0 & dat_csv$covariance_pvalue[i] >  0.05){dat_csv$Covconfintperm[i] = "True Negative"
    }else{dat_csv$Covconfintperm[i] = "None"}
}

## Covariance Bootstrap 
dat_csv$Covconfintboot = rep(NA, nrow(dat_csv))
for(i in 1:nrow(dat_csv)){
  if(dat_csv$true_cov[i] != 0 &&
     dat_csv$covariance_lwrCI[i] < 0 &&
     dat_csv$covariance_uprCI[i] < 0 
  ){dat_csv$Covconfintboot[i] = "True Positive"
  }else if(dat_csv$true_cov[i] != 0 &&
           dat_csv$covariance_lwrCI[i] > 0 &&
           dat_csv$covariance_uprCI[i] > 0
  ){dat_csv$Covconfintboot[i] = "True Positive"
  }else if(dat_csv$true_cov[i] == 0 &&
           dat_csv$covariance_lwrCI[i] < 0 &&
           dat_csv$covariance_uprCI[i] < 0
  ){dat_csv$Covconfintboot[i] = "False Positive"
  }else if(dat_csv$true_cov[i] == 0 &&
           dat_csv$covariance_lwrCI[i] > 0 &&
           dat_csv$covariance_uprCI[i] > 0
  ){dat_csv$Covconfintboot[i] = "False Positive"
  }else if(dat_csv$true_cov[i] != 0 && 
           dat_csv$covariance_lwrCI[i] <= 0 && 
           dat_csv$covariance_uprCI[i] >= 0
  ){dat_csv$Covconfintboot[i] = "False Negative"
  }else if(dat_csv$true_cov[i]== 0 && 
           dat_csv$covariance_lwrCI[i] <= 0 && 
           dat_csv$covariance_uprCI[i] >= 0
  ){dat_csv$Covconfintboot[i] = "True Negative"
  }else{dat_csv$Covconfintboot[i] = "None"}
}

## GxE Permutation 
dat_csv$GxEconfintperm = rep("NA", nrow(dat_csv))
for(i in 1:nrow(dat_csv)){
  if(dat_csv$true_GxE_emm[i] != 0 & round(dat_csv$GxE_emm_pvalue[i],2) <= 0.05){dat_csv$GxEconfintperm[i] = "True Positive"
  }else if(dat_csv$true_GxE_emm[i] == 0 & round(dat_csv$GxE_emm_pvalue[i],2) <= 0.05){dat_csv$GxEconfintperm[i] = "False Positive"
  }else if(dat_csv$true_GxE_emm[i] != 0 & round(dat_csv$GxE_emm_pvalue[i],2) > 0.05){dat_csv$GxEconfintperm[i] = "False Negative"
  }else if(dat_csv$true_GxE_emm[i] == 0 & round(dat_csv$GxE_emm_pvalue[i],2) > 0.05){dat_csv$GxEconfintperm[i] = "True Negative"
  }else{dat_csv$GxEconfintperm == "None"}
}

## GxE Bootstrap 
dat_csv$GxEconfintboot = rep("NA", nrow(dat_csv))
for(i in 1:nrow(dat_csv)){
  if(dat_csv$true_GxE_emm[i] != 0 & round(dat_csv$GxE_emm_lwrCI[i],2) > 0)
  {dat_csv$GxEconfintboot[i] = "True Positive"
  }else if(dat_csv$true_GxE_emm[i] == 0 & round(dat_csv$GxE_emm_lwrCI[i],2) > 0)
  {dat_csv$GxEconfintboot[i] = "False Positive"
  }else if(dat_csv$true_GxE_emm[i] != 0 & round(dat_csv$GxE_emm_lwrCI[i],2) == 0)
  {dat_csv$GxEconfintboot[i] = "False Negative"
  }else if(dat_csv$true_GxE_emm[i] == 0 & round(dat_csv$GxE_emm_lwrCI[i],2) == 0)
  {dat_csv$GxEconfintboot[i] = "True Negative"
  }else{dat_csv$GxEconfintboot[i] = "None"}
}

## GxE Anova 
dat_csv$GxEanova_conf = rep("NA", nrow(dat_csv))
for(i in 1:nrow(dat_csv)){
  if(dat_csv$true_GxE_emm[i] != 0 & round(dat_csv$GxE_Anova[i],2) <= 0.05){dat_csv$GxEanova_conf[i] = "True Positive"
  }else if(dat_csv$true_GxE_emm[i] == 0 & round(dat_csv$GxE_Anova[i],2) <= 0.05){dat_csv$GxEanova_conf[i] = "False Positive"
  }else if(dat_csv$true_GxE_emm[i] != 0 & round(dat_csv$GxE_Anova[i],2) > 0.05){dat_csv$GxEanova_conf[i] = "False Negative"
  }else if(dat_csv$true_GxE_emm[i] == 0 & round(dat_csv$GxE_Anova[i],2) > 0.05){dat_csv$GxEanova_conf[i] = "True Negative"
  }else{dat_csv$GxEanova_conf == "None"}
}

## Filter Data

dat_csv1 <- dat_csv %>% 
  filter(std_dev == 1) %>%
  filter(total_samples > 17)

dat_csv2 <- dat_csv %>% 
  filter(std_dev == 1)%>%
  filter(total_samples > 17)


####### Reciprocal Transplant - Overall FPR/FNR  ######
cov_perm_table = dat_csv1 %>%
    group_by("name" = Covconfintperm) %>%
    summarize("n" = n())
#fpr.fnr(cov_perm_table, divided = FALSE, scenario = 1)

cov_boot_table = dat_csv1 %>%
    group_by("name" = Covconfintboot) %>%
  summarize("n" = n())
#fpr.fnr(cov_boot_table, divided = FALSE, scenario = 1)

gxe_anova_table = dat_csv2 %>%
  group_by("name" =GxEanova_conf) %>%
  summarize("n" = n())
#fpr.fnr(gxe_anova_table, divided = FALSE, scenario = 1)

gxe_perm_table = dat_csv2 %>%
    group_by("name" =GxEconfintperm) %>%
    summarize("n" = n())
#fpr.fnr(gxe_perm_table, divided = FALSE, scenario = 1)

gxe_boot_table = dat_csv2 %>%
  group_by("name" =GxEconfintboot) %>%
  summarize("n" = n())
#fpr.fnr(gxe_boot_table, divided = FALSE, scenario = 1)

####### Reciprocal Transplant - False Positive Tables ########
(raw_confusion_hmap1 = dat_csv1 %>%
    group_by(sample_size, n_pop, "name" = Covconfintperm) %>%
    summarize("n" = n()))
raw_conf1 = fpr.fnr(raw_confusion_hmap1, divided = TRUE, scenario = 1)
#raw_conf_plot1 <- heatmap_fun(raw_conf1,"rate") # Can also do "percent"

(raw_confusion_hmap2 = dat_csv1 %>%
    group_by(sample_size, n_pop, "name" =Covconfintboot) %>%
    summarize("n" = n()))
raw_conf2 = fpr.fnr(raw_confusion_hmap2, divided = TRUE, scenario = 1)
#raw_conf_plot2 <- heatmap_fun(raw_conf2,"rate")

(raw_confusion_hmap3 = dat_csv2 %>%
    group_by(sample_size, n_pop,"name" = GxEconfintboot) %>%
    summarize("n" = n()))
raw_conf3 = fpr.fnr(raw_confusion_hmap3, divided = TRUE, scenario = 1)
#raw_conf_plot3 <- heatmap_fun(raw_conf3,"rate")

raw_confusion_hmap4 = dat_csv2 %>%
    group_by(sample_size, n_pop,"name" = GxEconfintperm) %>%
    summarize("n" = n())
raw_conf4 = fpr.fnr(raw_confusion_hmap4, divided = TRUE, scenario = 1)
#raw_conf_plot4 <- heatmap_fun(raw_conf4,"rate")

raw_confusion_hmap5 = dat_csv2 %>%
  group_by(sample_size, n_pop,"name" = GxEanova_conf) %>%
  summarize("n" = n())
raw_conf5 = fpr.fnr(raw_confusion_hmap5, divided = TRUE, scenario = 1)
#raw_conf_plot5 <- heatmap_fun(raw_conf5,"rate")

# Compile FPs for plot
raw_conf1$ID = rep("Cov_Perm", nrow(raw_conf1))
raw_conf2$ID = rep("Cov_Boot", nrow(raw_conf2))
raw_conf3$ID = rep("GxE_Boot", nrow(raw_conf3))
raw_conf4$ID = rep("GxE_Perm", nrow(raw_conf4))
raw_conf5$ID = rep("GxE_Anova", nrow(raw_conf5))

fpdf = rbind(raw_conf1,raw_conf2,raw_conf3,raw_conf4,raw_conf5)

gxeFPR = rbind(raw_conf4,raw_conf5)
gxeFPR = gxeFPR[gxeFPR$name == "False Positive",]

covFPR = rbind(raw_conf1,raw_conf2)
covFPR = covFPR[covFPR$name == "False Positive",]

## False Positive Rates
fpdf$npop_plot = NA
for(i in 1:nrow(fpdf)){
  if(fpdf$n_pop[i] == 2){fpdf$npop_plot[i] = "2 Genotypes"
  }else if(fpdf$n_pop[i] == 4){fpdf$npop_plot[i] = "4 Genotypes"
  }else{fpdf$npop_plot[i] = "8 Genotypes"}
}
fpdf1 = fpdf[fpdf$name == "False Positive",]

####### Reciprocal Transplant - False Negative Tables ########

covperm1 = fnr.effsize(dat_csv1, metric = "cov", data.type = "raw", analysis ="perm",resolution = "fine")
covperm1$ID = rep("Cov_Perm",nrow(covperm1))
covboot1 = fnr.effsize(dat_csv1, metric = "cov", data.type = "raw", analysis ="boot",resolution = "fine")
covboot1$ID = rep("Cov_Boot",nrow(covboot1))
gxeperm1 = fnr.effsize(dat_csv1, metric = "gxe", data.type = "raw", analysis ="perm",resolution = "fine")
gxeperm1$ID = rep("GxE_Perm",nrow(gxeperm1))
gxeboot1 = fnr.effsize(dat_csv1, metric = "gxe", data.type = "raw", analysis = "boot",resolution = "fine")
gxeboot1$ID = rep("GxE_Boot",nrow(gxeboot1))
gxeanova1 = fnr.effsize(dat_csv1, metric = "gxe", data.type = "raw", analysis = "anova",resolution = "fine")
gxeanova1$ID = rep("GxE_Anova",nrow(gxeanova1))

fndf = rbind(covperm1,covboot1,gxeperm1,gxeboot1,gxeanova1)
gxePow = rbind(gxeperm1,gxeanova1)
covPow = rbind(covperm1,covboot1)


########################################################
##          Error Matrices  -- Common Garden          ##
########################################################

# CovGE -- Permutation - Scenario 2
dat_dub$Covconfintperm = rep("NA",nrow(dat_dub))

for(i in 1:nrow(dat_dub)){
  if(dat_dub$true_cov[i] != 0 && dat_dub$covariance_pvalue[i] <= 0.05){dat_dub$Covconfintperm[i] = "True Positive"
  }else if(dat_dub$true_cov[i] == 0 & dat_dub$covariance_pvalue[i] <= 0.05){dat_dub$Covconfintperm[i] = "False Positive"
  }else if(dat_dub$true_cov[i]!= 0 & dat_dub$covariance_pvalue[i] > 0.05){dat_dub$Covconfintperm[i] = "False Negative"
  }else if(dat_dub$true_cov[i] == 0 & dat_dub$covariance_pvalue[i] > 0.05){dat_dub$Covconfintperm[i] = "True Negative"
  }else{dat_dub$Covconfintperm[i] = "None"}
}

# Cov Boot check
dat_dub$Covconfintboot = rep(NA, nrow(dat_dub))

for(i in 1:nrow(dat_dub)){
  
  if(dat_dub$true_cov[i] != 0 &&
     dat_dub$covariance_lwrCI[i] < 0 &&
     dat_dub$covariance_uprCI[i] < 0
  ){dat_dub$Covconfintboot[i] = "True Positive"
  }else if(dat_dub$true_cov[i] != 0 &&
           dat_dub$covariance_lwrCI[i] > 0 &&
           dat_dub$covariance_uprCI[i] > 0
  ){dat_dub$Covconfintboot[i] = "True Positive"
  }else if(dat_dub$true_cov[i] == 0 &&
           dat_dub$covariance_lwrCI[i] < 0 &&
           dat_dub$covariance_uprCI[i] < 0
  ){dat_dub$Covconfintboot[i] = "False Positive"
  }else if(dat_dub$true_cov[i] == 0 &&
           dat_dub$covariance_lwrCI[i] > 0 &&
           dat_dub$covariance_uprCI[i] > 0
  ){dat_dub$Covconfintboot[i] = "False Positive"
  }else if(dat_dub$true_cov[i] != 0 && 
           dat_dub$covariance_lwrCI[i] <= 0 && 
           dat_dub$covariance_uprCI[i] >= 0
  ){dat_dub$Covconfintboot[i] = "False Negative"
  }else if(dat_dub$true_cov[i]== 0 && 
           dat_dub$covariance_lwrCI[i] <= 0 && 
           dat_dub$covariance_uprCI[i] >= 0
  ){dat_dub$Covconfintboot[i] = "True Negative"
  }else{dat_dub$Covconfintboot[i] = "None"}
}

# GxE Anova check
dat_dub$GxEanova_conf = rep("NA", nrow(dat_dub))

for(i in 1:nrow(dat_dub)){
  if(dat_dub$true_GxE_emm[i] != 0 & round(dat_dub$GxE_Anova[i],2) <= 0.05){dat_dub$GxEanova_conf[i] = "True Positive"
  }else if(dat_dub$true_GxE_emm[i] == 0 & round(dat_dub$GxE_Anova[i],2) <= 0.05){dat_dub$GxEanova_conf[i] = "False Positive"
  }else if(dat_dub$true_GxE_emm[i] != 0 & round(dat_dub$GxE_Anova[i],2) > 0.05){dat_dub$GxEanova_conf[i] = "False Negative"
  }else if(dat_dub$true_GxE_emm[i] == 0 & round(dat_dub$GxE_Anova[i],2) > 0.05){dat_dub$GxEanova_conf[i] = "True Negative"
  }else{dat_dub$GxEanova_conf == "None"}
}

# GxE Perm check
dat_dub$GxEconfintperm = rep("NA", nrow(dat_dub))

for(i in 1:nrow(dat_dub)){
  if(dat_dub$true_GxE_emm[i] != 0 & round(dat_dub$GxE_emm_pvalue[i],2) <= 0.05){dat_dub$GxEconfintperm[i] = "True Positive"
  }else if(dat_dub$true_GxE_emm[i] == 0 & round(dat_dub$GxE_emm_pvalue[i],2) <= 0.05){dat_dub$GxEconfintperm[i] = "False Positive"
  }else if(dat_dub$true_GxE_emm[i] != 0 & round(dat_dub$GxE_emm_pvalue[i],2) > 0.05){dat_dub$GxEconfintperm[i] = "False Negative"
  }else if(dat_dub$true_GxE_emm[i] == 0 & round(dat_dub$GxE_emm_pvalue[i],2) > 0.05){dat_dub$GxEconfintperm[i] = "True Negative"
  }else{dat_dub$GxEconfintperm == "None"}
}

# GxE Boot check
dat_dub$GxEconfintboot = rep("NA", nrow(dat_dub))

for(i in 1:nrow(dat_dub)){
  if(dat_dub$true_GxE_emm[i] != 0 & round(dat_dub$GxE_emm_lwrCI[i],2) > 0)
  {dat_dub$GxEconfintboot[i] = "True Positive"
  }else if(dat_dub$true_GxE_emm[i] == 0 & round(dat_dub$GxE_emm_lwrCI[i],2) > 0)
  {dat_dub$GxEconfintboot[i] = "False Positive"
  }else if(dat_dub$true_GxE_emm[i] != 0 & round(dat_dub$GxE_emm_lwrCI[i],2) == 0)
  {dat_dub$GxEconfintboot[i] = "False Negative"
  }else if(dat_dub$true_GxE_emm[i] == 0 & round(dat_dub$GxE_emm_lwrCI[i],2) == 0)
  {dat_dub$GxEconfintboot[i] = "True Negative"
  }else{dat_dub$GxEconfintboot[i] = "None"}
}

dat_dub1 <- dat_dub %>% 
  filter(std_dev == 1) %>%
  #filter(between(abs(covariance),0.2,0.6))%>%
  filter(total_samples > 17)
  #filter(total_samples > 128)
dat_dub2 <- dat_dub %>% 
  filter(std_dev == 1) %>%
 # filter(between(GxE_emm,0.3,0.6)) %>%
  filter(total_samples > 17)
 #filter(total_samples > 128)


####### Common Garden - Overall FPR/FNR   ######
cov_perm_table = dat_dub1 %>%
  group_by("name" = Covconfintperm) %>%
  summarize("n" = n())
#fpr.fnr(cov_perm_table, divided = FALSE, scenario = 2)

cov_boot_table = dat_dub1 %>%
  group_by("name" =Covconfintboot) %>%
  summarize("n" = n())
#fpr.fnr(cov_boot_table, divided = FALSE, scenario = 2)

gxe_anova_table = dat_dub2 %>%
  group_by("name" = GxEanova_conf) %>%
  summarize("n" = n())
#fpr.fnr(gxe_anova_table, divided = FALSE, scenario = 2)

gxe_perm_table = dat_dub2 %>%
  group_by("name" =GxEconfintperm) %>%
  summarize("n" = n())
#fpr.fnr(gxe_perm_table, divided = FALSE, scenario = 2)

gxe_boot_table = dat_dub2 %>%
  group_by("name" =GxEconfintboot) %>%
  summarize("n" = n())
#fpr.fnr(gxe_boot_table, divided = FALSE, scenario = 2)

####### Common Garden - False Positive Tables  #####
(raw_confusion_hmap1 = dat_dub1 %>%
    group_by(sample_size, n_pop, "name" =Covconfintperm) %>%
    summarize("n" = n()))
dub_raw_conf1 = fpr.fnr(raw_confusion_hmap1, divided = TRUE, scenario = 2)
#raw_conf_plot1 <- heatmap_fun(dub_raw_conf1,"rate") # Can also do "percent"

(raw_confusion_hmap2 = dat_dub1 %>%
    group_by(sample_size, n_pop, "name" =Covconfintboot) %>%
    summarize("n" = n()))
dub_raw_conf2 = fpr.fnr(raw_confusion_hmap2, divided = TRUE, scenario = 2)
#raw_conf_plot2 <- heatmap_fun(dub_raw_conf2,"rate")

(raw_confusion_hmap3 = dat_dub2 %>%
    group_by(sample_size, n_pop,"name" = GxEconfintboot) %>%
    summarize("n" = n()))
dub_raw_conf3 = fpr.fnr(raw_confusion_hmap3, divided = TRUE, scenario = 2)
#raw_conf_plot3 <- heatmap_fun(dub_raw_conf3,"rate")

raw_confusion_hmap4 = dat_dub2 %>%
  group_by(sample_size, n_pop,"name" = GxEconfintperm) %>%
  summarize("n" = n())
dub_raw_conf4 = fpr.fnr(raw_confusion_hmap4, divided = TRUE, scenario = 2)
#raw_conf_plot4 <- heatmap_fun(dub_raw_conf4,"rate")

raw_confusion_hmap5 = dat_dub2 %>%
  group_by(sample_size, n_pop,"name" = GxEanova_conf) %>%
  summarize("n" = n())
dub_raw_conf5 = fpr.fnr(raw_confusion_hmap5, divided = TRUE, scenario = 2)
#raw_conf_plot5 <- heatmap_fun(dub_raw_conf5,"rate")

# Compile FPs for plot
dub_raw_conf1$ID = rep("Cov_Perm", nrow(dub_raw_conf1))
dub_raw_conf2$ID = rep("Cov_Boot", nrow(dub_raw_conf2))
dub_raw_conf3$ID = rep("GxE_Boot", nrow(dub_raw_conf3))
dub_raw_conf4$ID = rep("GxE_Perm", nrow(dub_raw_conf4))
dub_raw_conf5$ID = rep("GxE_Anova", nrow(dub_raw_conf5))
dub_fpdf = rbind(dub_raw_conf1,dub_raw_conf2,dub_raw_conf3,dub_raw_conf4,dub_raw_conf5)
dub_fpdf = dub_fpdf[dub_fpdf$name == "False Positive",]

dub_gxe = rbind(dub_raw_conf3,dub_raw_conf4,dub_raw_conf5)
dub_gxe = dub_gxe %>% filter(name == "False Positive") 
dub_cov = rbind(dub_raw_conf1,dub_raw_conf2)
dub_cov = dub_cov[dub_cov$name == "False Positive",]

dub_raw_conf3[is.nan(dub_raw_conf3)] <- 0
dub_raw_conf3 %>%
  filter(name == "False Positive") %>%
  summarize(mean(rate))

dub_fpdf$npop_plot = NA
for(i in 1:nrow(dub_fpdf)){
  if(dub_fpdf$n_pop[i] == 4){dub_fpdf$npop_plot[i] = "4 Genotypes"
  }else if(dub_fpdf$n_pop[i] == 8){dub_fpdf$npop_plot[i] = "8 Genotypes"
  }else{dub_fpdf$npop_plot[i] = "16 Genotypes"}
}
dub_fpdf$npop_plot = factor(dub_fpdf$npop_plot, levels=c("4 Genotypes","8 Genotypes","16 Genotypes")) 


####### Common Garden - False Negative Tables  #####
covperm2 = fnr.effsize(dat_dub1, metric = "cov", data.type = "raw", analysis ="perm",scenario = 2,resolution = "fine")
covperm2$ID = rep("Cov_Perm",nrow(covperm2))
covboot2 = fnr.effsize(dat_dub1, metric = "cov", data.type = "raw", analysis ="boot",scenario = 2,resolution = "fine")
covboot2$ID = rep("Cov_Boot",nrow(covboot2))
gxeperm2 = fnr.effsize(dat_dub1, metric = "gxe", data.type = "raw", analysis ="perm",scenario = 2,resolution = "fine")
gxeperm2$ID = rep("GxE_Perm",nrow(gxeperm2))
gxeboot2 = fnr.effsize(dat_dub1, metric = "gxe", data.type = "raw", analysis = "boot",scenario = 2,resolution = "fine")
gxeboot2$ID = rep("GxE_Boot",nrow(gxeboot2))
gxeanova2 = fnr.effsize(dat_dub1, metric = "gxe", data.type = "raw", analysis = "anova",scenario = 2,resolution = "fine")
gxeanova2$ID = rep("GxE_Anova",nrow(gxeanova2))

fndfdub = rbind(covperm2,covboot2,gxeperm2,gxeboot2,gxeanova2)
gxePowdub = rbind(gxeperm2,gxeanova2)
covPowdub = rbind(covperm2,covboot2)




###### Reciprocal Transplant - CovGE - False Negative Heatmap  ######
covPow[is.nan(covPow)] <- 0
covPow1 <- 
  covPow %>% 
  filter(between(bin,0.2,0.6))%>%
  group_by(sample_size,n_pop,ID) %>%
  summarize("FNR" = mean(fnr))
covPow1$n_pop<-as.factor(covPow1$n_pop)
covPow2<-
  covPow1 %>%
  mutate(#xadj = ifelse(ID %in% c("Cov_Perm"), -.45/2,  .45/2),
    yadj = ifelse(ID %in% c("Cov_Perm"), -.475/2,  .475/2),
    xpos = as.numeric(factor(sample_size)),
    ypos = as.numeric(n_pop) + yadj,
    ID2 = ifelse(ID %in% c("Cov_Perm"), "Perm.", "Boot."))

(covFNFRT = covPow2 %>% 
    ggplot(aes(as.factor(sample_size), ypos, fill = FNR)) + 
    geom_tile(height = 0.475, width = 0.95, color= "white") +
    geom_text(aes(label = paste0(ID2,"\n",round(FNR, 2))), size = 4, family = "Times", colour = "white", show.legend = F) +
    #ggtitle("False Negative Rates: CovGE")+
    scale_fill_viridis(breaks=seq(0,1,0.25), #breaks in the scale bar
                       limits=c(0,1))+
    scale_x_discrete(name = "Sample Size",
                     labels = unique(covPow1$sample_size)) +
    scale_y_continuous(breaks = 1:3, 
                       name = "Number of Genotypes",
                       labels = c(2,4,8)) +
    labs(fill = "False Negative Rate")+
    #  ggtitle(expression("Reciprocal Transplant: Cov"["GE"]*" False Negative Rates"))+
    theme_classic(base_size = 18, base_family = "Times")+
    theme(axis.text = element_text(colour = "black"))+
    theme(legend.position = "none")+
    ggtitle("FRT: False Negative Rates")+
    theme(plot.title = element_text(size = 18, face = "bold")))

###### Reciprocal Transplant - CovGE - False Positive Heatmap #######
covFPR[is.nan(covFPR)] <- 0
covFPR1 <- 
  covFPR %>% 
  mutate(#xadj = ifelse(ID %in% c("Cov_Perm"), -.45/2,  .45/2),
    yadj = ifelse(ID %in% c("Cov_Perm"), -.475/2,  .475/2),
    xpos = as.numeric(factor(sample_size)),
    ypos = as.numeric(factor(n_pop)) + yadj,
    ID2 = ifelse(ID %in% c("Cov_Perm"), "Perm.", "Boot."))

(covFPFRT = covFPR1 %>% 
    ggplot(aes(as.factor(sample_size), ypos, fill = rate)) + 
    geom_tile(height = 0.475, width = 0.95, color= "white") +
    geom_text(aes(label = paste0(ID2,"\n",round(rate, 2))), size = 4, family = "Times", colour = "white", show.legend = F) +
    #ggtitle("False Positive Rates: CovGE")+
    scale_fill_viridis(breaks=seq(0,.5,0.1), #breaks in the scale bar
                       limits=c(0,.4))+
    scale_x_discrete(name = "Sample Size",
                     labels = unique(covFPR1$sample_size)) +
    scale_y_continuous(breaks = 1:3, name = "Number of Genotypes",
                       labels = c(2,4,8)) +
    labs(fill = "False Positive Rate")+
    theme_classic(base_size = 18, base_family = "Times")+
    theme(axis.text = element_text(colour = "black"))+
    theme(legend.position = "none")+
    ggtitle(expression("A    Cov"["GE"]*": Reciprocal Transplant"))+
    theme(plot.title = element_text(size = 18, face = "bold"))) 

###### Reciprocal Transplant - GxE - False Positive Heatmap  #######
gxeFPR[is.nan(gxeFPR)] <- 0
gxeFPR1 <- 
  gxeFPR %>% 
  mutate(#xadj = ifelse(ID == "GxE_Anova", -.45/3,ifelse(ID == "GxE_Perm", 0.0, .45/3)),
    yadj = ifelse(ID == "GxE_Anova", -.475/2,  .475/2),
    #ifelse(ID == "GxE_Perm", 0.0, 0.9/3)),
    xpos = as.numeric(factor(sample_size)),
    ypos = as.numeric(factor(n_pop))+yadj,
    col = ifelse(rate >= 0.5, "black","white"),
    ID2 = ifelse(ID == "GxE_Anova", "Anova",  
                 ifelse(ID == "GxE_Perm", "Perm.", "Boot.")))

(gxeFPFRT = gxeFPR1 %>% 
    ggplot(aes(as.factor(sample_size), ypos, fill = rate)) + 
    geom_tile(height = .475, width = 0.95, color= "white") +
    geom_text(aes(label = paste0(ID2,"\n",round(rate, 2)), colour = col), 
              size = 4, family = "Times", show.legend = F) +
    scale_colour_identity()+
    scale_fill_viridis(breaks=seq(0,.5,0.1), #breaks in the scale bar
                       limits=c(0,.4))+
    scale_x_discrete(name = "Sample Size",
                     labels = unique(gxeFPR1$sample_size)) +
    scale_y_continuous(breaks = 1:3, name = "Number of Genotypes",
                       labels = c(2,4,8)) +
    labs(fill = "False Positive Rate")+
    theme_classic(base_size = 18, base_family = "Times")+
    theme(axis.text = element_text(colour = "black"))+
    theme(legend.position = "none")+
    ggtitle(expression("C   "*bar(Delta)*""["GxE"]*": Reciprocal Transplant"))+
    theme(plot.title = element_text(size = 18, face = "bold")))


###### Reciprocal Transplant - GxE - False Negative Heatmap #####
gxePow[is.nan(gxePow)] <- 0
gxeFNR1 = gxePow %>%
  filter(between(bin,0.3,0.6))%>%  
  
  group_by(sample_size,n_pop,ID) %>%
  summarize("fnr" = mean(fnr))
gxeFNR1$n_pop<-as.factor(gxeFNR1$n_pop)
gxeFNR2 <- 
  gxeFNR1 %>% 
  mutate(#xadj = ifelse(ID == "GxE_Anova", -.45/3,ifelse(ID == "GxE_Perm", 0.0, .45/3)),
    yadj = ifelse(ID == "GxE_Anova", -.475/2,  .475/2),
    # ifelse(ID == "GxE_Perm", 0.0, 0.9/3)),
    xpos = as.numeric(factor(sample_size)),
    ypos = as.numeric(n_pop)+yadj,
    col = ifelse(fnr >= 0.5, "black","white"),
    ID2 = ifelse(ID == "GxE_Anova", "Anova",  
                 ifelse(ID == "GxE_Perm", "Perm.", "Boot.")))

(gxeFNFRT = gxeFNR2 %>% 
    ggplot(aes(as.factor(sample_size), ypos, fill = fnr)) + 
    geom_tile(height = 0.475, width = 0.95, color= "white") +
    geom_text(aes(label = paste0(ID2,"\n",round(fnr, 2)), colour = col), 
              size = 4, family = "Times", show.legend = F) +
    scale_colour_identity()+
    scale_fill_viridis(breaks=seq(0,1,0.25), #breaks in the scale bar
                       limits=c(0,1))+
    scale_x_discrete(name = "Sample Size",
                     labels = unique(gxeFNR2$sample_size)) +
    scale_y_continuous(breaks = 1:3, name = "Number of Genotypes",
                       labels = c(2,4,8)) +
    # ggtitle(expression("Reciprocal Transplant: "*bar(Delta)*""["GxE"]*" False Positive Rates"))+  
    labs(fill = "False Negative Rate")+
    theme_classic(base_size = 18, base_family = "Times")+
    theme(axis.text = element_text(colour = "black"))+
    theme(legend.position = "none")+
    ggtitle("FRT: False Negative Rates")+
    theme(plot.title = element_text(size = 18, face = "bold")))

# False Positive Panel
# grid.arrange(covFPFRT, covFPCG, gxeFPFRT, gxeFPCG, ncol = 2)


###### Common Garden - CovGE - False Negative Heatmap #####
covPowdub[is.nan(covPowdub)] <- 0
covPowdub1 <- 
  covPowdub %>% 
  filter(between(bin,0.2,0.6))%>%
  group_by(sample_size,n_pop,ID) %>%
  summarize("FNR" = mean(fnr))
covPowdub1$n_pop<-as.factor(covPowdub1$n_pop)
covPowdub2<-
  covPowdub1 %>%
  mutate(#xadj = ifelse(ID %in% c("Cov_Perm"), -.45/2,  .45/2),
    yadj = ifelse(ID %in% c("Cov_Perm"), -.475/2,  .475/2),
    xpos = as.numeric(factor(sample_size)),
    ypos = as.numeric(n_pop) + yadj,
    ID2 = ifelse(ID %in% c("Cov_Perm"), "Perm.", "Boot."))

(covFNCG = covPowdub2 %>% 
    ggplot(aes(as.factor(sample_size), ypos, fill = FNR)) + 
    geom_tile(height = 0.475, width = 0.95, color= "white") +
    geom_text(aes(label = paste0(ID2,"\n",round(FNR, 2))), size = 4, family = "Times", colour = "white", show.legend = F) +
    scale_fill_viridis(breaks=seq(0,1,0.25), #breaks in the scale bar
                       limits=c(0,1))+
    scale_x_discrete(name = "Sample Size",
                     labels = unique(covPowdub2$sample_size)) +
    scale_y_continuous(breaks = 1:3, 
                       name = "Number of Genotypes",
                       labels = c(2,4,8)) +
    labs(fill = "False Negative Rate")+
    #  ggtitle(expression("Common Garden: Cov"["GE"]*" False Negative Rates"))+
    theme_classic(base_size = 18, base_family = "Times")+
    theme(axis.text = element_text(colour = "black"))+
    theme(legend.position = "none")+
    ggtitle("CG: False Negative Rates")+
    theme(plot.title = element_text(size = 18, face = "bold")))

###### Common Garden - CovGE - False Positive Heatmap #####
dub_cov[is.nan(dub_cov)] <- 0
dub_cov1 <- 
  dub_cov %>% 
  mutate(#xadj = ifelse(ID %in% c("Cov_Perm"), -.45/2,  .45/2),
    yadj = ifelse(ID %in% c("Cov_Perm"), -.475/2,  .475/2),
    xpos = as.numeric(factor(sample_size)),
    ypos = as.numeric(factor(n_pop)) + yadj,
    col = ifelse(rate >= 0.29, "black","white"),
    ID2 = ifelse(ID %in% c("Cov_Perm"), "Perm.", "Boot."))

(covFPCG = dub_cov1 %>% 
  ggplot(aes(as.factor(sample_size), ypos, fill = rate)) + 
  geom_tile(height = 0.475, width = 0.95, color= "white") +
  geom_text(aes(label = paste0(ID2,"\n",round(rate, 2)), colour = col), size =4, family = "Times", show.legend = F) +
  #ggtitle("False Positive Rates: CovGE")+
    scale_fill_viridis(breaks=seq(0,.5,0.1), #breaks in the scale bar
                       limits=c(0,.4))+
    scale_colour_identity()+
  scale_x_discrete(name = "Sample Size",
                   labels = unique(dub_cov1$sample_size)) +
  scale_y_continuous(breaks = 1:3, name = "Number of Genotypes",
                     labels = c(4,8,16)) +
   # ggtitle(expression("Common Garden: Cov"["GE"]*" False Positive Rates"))+
  labs(fill = "False Positive Rate")+
  theme_classic(base_size = 18, base_family = "Times")+
  theme(axis.text = element_text(colour = "black"))+
  theme(legend.position = "none")+
  ggtitle(expression("B    Cov"["GE"]*": Common Garden"))+
  theme(plot.title = element_text(size = 18, face = "bold")))


###### Common Garden - GxE - False Positive Heatmap #####
dub_gxe[is.nan(dub_gxe)] <- 0
dub_gxe1 <- 
  dub_gxe %>% 
  mutate(#xadj = ifelse(ID == "GxE_Anova", -.45/3,ifelse(ID == "GxE_Perm", 0.0, .45/3)),
    yadj = ifelse(ID == "GxE_Anova",  -.475/2,  .475/2),
                  #ifelse(ID == "GxE_Perm", 0.0, 0.9/3)),
    xpos = as.numeric(factor(sample_size)),
    ypos = as.numeric(factor(n_pop))+yadj,
    col = ifelse(rate >= 0.5, "black","white"),
    ID2 = ifelse(ID == "GxE_Anova", "Anova",  
                 ifelse(ID == "GxE_Perm", "Perm.", "Boot.")))

(gxeFPCG = dub_gxe1 %>% 
  ggplot(aes(as.factor(sample_size), ypos, fill = rate)) + 
  geom_tile(height = .475, width = 0.95, color= "white") +
  geom_text(aes(label = paste0(ID2,"\n",round(rate, 2)), colour = col), 
            size = 4, family = "Times", show.legend = F) +
  scale_colour_identity()+
  scale_fill_viridis(breaks=seq(0,.5,0.1), #breaks in the scale bar
                     limits=c(0,.4))+
  scale_x_discrete(name = "Sample Size",
                   labels = unique(gxeFPR1$sample_size)) +
  scale_y_continuous(breaks = 1:3, name = "Number of Genotypes",
                     labels = c(4,8,16)) +
    #ggtitle(expression("Common Garden: "*bar(Delta)*""["GxE"]*" False Positive Rates"))+
  labs(fill = "False Positive Rate")+
  theme_classic(base_size = 18, base_family = "Times")+
  theme(axis.text = element_text(colour = "black"))+
  theme(legend.position = "none")+
  ggtitle(expression("D    "*bar(Delta)*""["GxE"]*": Common Garden"))+
  theme(plot.title = element_text(size = 18, face = "bold"))) 




###### Recip. Transplant - CovGE - Power Heatmap ######
covPow[is.nan(covPow)] <- 0
covPow1 <- #
  covPow %>% 
  filter(between(bin,0.3,0.6))%>%
  group_by(sample_size,n_pop,ID) %>%
  summarize("FNR" = mean(fnr))
covPow1$totals = covPow1$n_pop * covPow1$n_pop * covPow1$sample_size
covPow1$power = 1-covPow1$FNR

covPow1$n_pop<-as.factor(covPow1$n_pop)
covPow2<-
  covPow1 %>%
  mutate(#xadj = ifelse(ID %in% c("Cov_Perm"), -.45/2,  .45/2),
    yadj = ifelse(ID %in% c("Cov_Perm"), -.475/2,  .475/2),
    xpos = as.numeric(factor(sample_size)),
    ypos = as.numeric(n_pop) + yadj,
    col = ifelse(power >= 0.5, "black","white"),
    ID2 = ifelse(ID %in% c("Cov_Perm"), "Perm.", "Boot."))

(covFRT_power = covPow2 %>% 
    ggplot(aes(as.factor(sample_size), ypos, fill = power)) + 
    geom_tile(height = 0.475, width = 0.95, color= "white") +
    geom_text(aes(label = paste0(ID2,"\n",covPow2$totals, "\n",round(power, 2))), size = 4, family = "Times", colour = covPow2$col, show.legend = F) +
    scale_colour_identity()+    #ggtitle("False Negative Rates: CovGE")+
    scale_fill_viridis(breaks=seq(0,1,0.25), #breaks in the scale bar
                       limits=c(0,1))+
    scale_x_discrete(name = "Sample Size",
                     labels = unique(covPow1$sample_size)) +
    scale_y_continuous(breaks = 1:3, 
                       name = "Number of Genotypes",
                       labels = c(2,4,8)) +
    labs(fill = "False Negative Rate")+
    #  ggtitle(expression("Reciprocal Transplant: Cov"["GE"]*" False Negative Rates"))+
    theme_classic(base_size = 18, base_family = "Times")+
    theme(axis.text = element_text(colour = "black"))+
    theme(legend.position = "none")+
    #ggtitle("C    FRT: Power")+
    theme(plot.title = element_text(size = 18, face = "bold")))

###### Common Garden- CovGE - Power Heatmap ######
covPowdub[is.nan(covPowdub)] <- 0
covPowdub1 <- 
  covPowdub %>% 
  filter(between(bin,0.2,0.6))%>%
  group_by(sample_size,n_pop,ID) %>%
  summarize("FNR" = mean(fnr))
covPowdub1$power = 1-covPowdub1$FNR
covPowdub1$totals = covPowdub1$n_pop * 2 * covPowdub1$sample_size

covPowdub1$n_pop<-as.factor(covPowdub1$n_pop)
covPowdub2<-
  covPowdub1 %>%
  mutate(#xadj = ifelse(ID %in% c("Cov_Perm"), -.45/2,  .45/2),
    yadj = ifelse(ID %in% c("Cov_Perm"), -.475/2,  .475/2),
    xpos = as.numeric(factor(sample_size)),
    ypos = as.numeric(n_pop) + yadj,
    col = ifelse(power >= 0.5, "black","white"),
    ID2 = ifelse(ID %in% c("Cov_Perm"), "Perm.", "Boot."))

(covCGPower = covPowdub2 %>% 
    ggplot(aes(as.factor(sample_size), ypos, fill = power)) + 
    geom_tile(height = 0.475, width = 0.95, color= "white") +
    geom_text(aes(label = paste0(ID2,"\n",covPowdub2$totals, "\n",round(power, 2))), size = 4, family = "Times", colour = covPowdub2$col, show.legend = F) +
    scale_colour_identity()+
    scale_fill_viridis(breaks=seq(0,1,0.25), #breaks in the scale bar
                       limits=c(0,1))+
    scale_x_discrete(name = "Sample Size",
                     labels = unique(covPowdub2$sample_size)) +
    scale_y_continuous(breaks = 1:3, 
                       name = "Number of Genotypes",
                       labels = c(4,8,16)) +
    theme(axis.text.y = element_text(lineheight = .5, 
                                     size = 6))+
    labs(fill = "Power")+
    #  ggtitle(expression("Common Garden: Cov"["GE"]*" False Negative Rates"))+
    theme_classic(base_size = 18, base_family = "Times")+
    theme(axis.text = element_text(colour = "black"))+
    theme(legend.position = "none")+
    ggtitle("D    CG: Power")+
    theme(plot.title = element_text(size = 18, face = "bold")))


###### Reciprocal Transplant - GxE - Power Heatmap ######
gxePow[is.nan(gxePow)] <- 0
gxeFNR1 = gxePow %>%
  filter(between(bin,0.2,0.5))%>%  
  group_by(sample_size,n_pop,ID) %>%
  summarize("power" = mean(power))
#gxeFNR1$power = 1- gxeFNR1$fnr
gxeFNR1$totals = gxeFNR1$sample_size * gxeFNR1$n_pop * gxeFNR1$n_pop
gxeFNR1$n_pop<-as.factor(gxeFNR1$n_pop)

gxeFNR2 <- 
  gxeFNR1 %>% 
  mutate(#xadj = ifelse(ID == "GxE_Anova", -.45/3,ifelse(ID == "GxE_Perm", 0.0, .45/3)),
    yadj = ifelse(ID == "GxE_Anova", -.475/2,  .475/2),
    # ifelse(ID == "GxE_Perm", 0.0, 0.9/3)),
    xpos = as.numeric(factor(sample_size)),
    ypos = as.numeric(n_pop)+yadj,
    col = ifelse(power >= 0.5, "black","white"),
    ID2 = ifelse(ID == "GxE_Anova", "Anova",  
                 ifelse(ID == "GxE_Perm", "Perm.", "Boot.")))

(gxepowerFRT = gxeFNR2 %>% 
    ggplot(aes(as.factor(sample_size), ypos, fill = power)) + 
    geom_tile(height = 0.475, width = 0.95, color= "white") +
    geom_text(aes(label = paste0(ID2,"\n",gxeFNR2$totals, "\n",round(power, 2)), colour = col), 
              size = 4, family = "Times", show.legend = F) +
    scale_colour_identity()+
    scale_fill_viridis(breaks=seq(0,1,0.25), #breaks in the scale bar
                       limits=c(0,1))+
    scale_x_discrete(name = "Sample Size",
                     labels = unique(gxeFNR2$sample_size)) +
    scale_y_continuous(breaks = 1:3, name = "Number of Genotypes",
                       labels = c(2,4,8)) +
    # ggtitle(expression("Reciprocal Transplant: "*bar(Delta)*""["GxE"]*" False Positive Rates"))+  
    labs(fill = "False Negative Rate")+
    theme_classic(base_size = 18, base_family = "Times")+
    theme(axis.text = element_text(colour = "black"))+
    theme(legend.position = "none")+
    ggtitle("C   FRT: Power")+
    theme(plot.title = element_text(size = 18, face = "bold")))

###### Common Garden - GxE - Power Heatmap ######
gxePowdub[is.nan(gxePowdub)] <- 0
gxePowdub1 = gxePowdub %>%
  filter(between(bin,0.3,0.6))%>%
  group_by(sample_size,n_pop,ID) %>%
  summarize("fnr" = mean(fnr))
gxePowdub1$power = 1-gxePowdub1$fnr
gxePowdub1$totals = gxePowdub1$n_pop * 2 * gxePowdub1$sample_size
gxePowdub1$n_pop<-as.factor(gxePowdub1$n_pop)
gxePowdub2 <- 
  gxePowdub1 %>% 
  mutate(#xadj = ifelse(ID == "GxE_Anova", -.45/3,ifelse(ID == "GxE_Perm", 0.0, .45/3)),
    yadj = ifelse(ID == "GxE_Anova",  -.475/2,  .475/2),
    # ifelse(ID == "GxE_Perm", 0.0, 0.9/2)),
    xpos = as.numeric(factor(sample_size)),
    ypos = as.numeric(n_pop)+yadj,
    col = ifelse(power >= 0.5, "black","white"),
    ID2 = ifelse(ID == "GxE_Anova", "Anova",  
                 ifelse(ID == "GxE_Perm", "Perm.", "Boot.")))

(gxepower_dub  = gxePowdub2 %>% 
    ggplot(aes(as.factor(sample_size), ypos, fill = power)) + 
    geom_tile(height = 0.475, width = 0.95, color= "white") +
    geom_text(aes(label = paste0(ID2,"\n",gxePowdub2$totals, "\n",round(power, 2)), colour = col), 
              size = 4, family = "Times", show.legend = F) +
    scale_colour_identity()+
    scale_fill_viridis(breaks=seq(0,1,0.25), #breaks in the scale bar
                       limits=c(0,1))+
    scale_x_discrete(name = "Sample Size",
                     labels = unique(gxePowdub2$sample_size)) +
    scale_y_continuous(breaks = 1:3, name = "Number of Genotypes",
                       labels = c(4,8,16)) +
    #ggtitle(expression("Reciprocal Transplant: "*bar(Delta)*""["GxE"]*" False Positive Rates"))+  
    labs(fill = "Power")+
    theme_classic(base_size = 18, base_family = "Times")+
    theme(axis.text = element_text(colour = "black"))+
    theme(legend.position = "none")+
    ggtitle("D   CG: Power")+
    theme(plot.title = element_text(size = 18, face = "bold")))

#grid.arrange(eff_var1, gxeFPFRT, gxeFPCG, gxepowerFRT,gxepower_dub, #ncol = 4,
 #            layout_matrix = cbind(c(1,2,4),c(1,3,5)))




###### Reciprocal Transplant - CovGE & GxE - Power ~ Effect Size #####
frtGxE1 = fnr.effsize(dat_csv1, metric = "GxE", scenario = 1,  data.type = "raw", analysis ="perm",resolution = "fine")
frtGxE1$ID = rep("GxE_Perm", nrow(frtGxE1))
frtGxE2 = fnr.effsize(dat_csv1, metric = "GxE", scenario = 1,  data.type = "raw", analysis ="anova",resolution = "fine")
frtGxE2$ID = rep("GxE_Anova", nrow(frtGxE2))
frtcov1 = fnr.effsize(dat_csv1, metric = "Cov", scenario = 1,   data.type = "raw", analysis ="perm",resolution = "fine")
frtcov1$ID = rep("Cov_Perm", nrow(frtcov1))
frtcov2 = fnr.effsize(dat_csv1, metric = "Cov", scenario = 1,   data.type = "raw", analysis ="boot",resolution = "fine")
frtcov2$ID = rep("Cov_Boot", nrow(frtcov2))

trade2_frt = rbind(frtGxE1,frtGxE2, frtcov1,frtcov2)
trade2_frt$totals = trade2_frt$sample_size*trade2_frt$n_pop *trade2_frt$n_pop
#trade3_frt = trade2_frt %>%
#  group_by(ID,totals,bin)%>%
#  summarize("avg_power" = mean(power))
trade2_frt$grp = paste(trade2_frt$bin,trade2_frt$totals,trade2_frt$ID,sep = "_")
trade2_frt$linegrp = paste(trade2_frt$totals,trade2_frt$ID,sep = "_")

(Power256_FRT = ggplot(trade2_frt %>% filter(totals == 256) %>% filter(bin != 0.0),aes(x = bin, y = power)) + 
    geom_rect(aes(xmin="0.2", xmax="0.5", ymin=-0.05, ymax=1.00), alpha=0.002) +
    geom_smooth(aes(group = linegrp, colour = factor(ID),linetype = factor(ID)), se = F)+
    geom_point(aes(shape = factor(ID), group = grp, fill = factor(ID)),position = position_dodge(width = 0.4),size = 4,alpha = 0.7)+
    geom_hline(aes(yintercept = 0.8),linetype = "dashed")+
    theme_classic(base_size = 22, base_family = "Times")+
    scale_colour_manual(values = c("Cov_Boot" = "#3CBB75FF","Cov_Perm" = "#3CBB75FF",
                                   "GxE_Anova" = "#453781FF","GxE_Perm" = "#453781FF"),#"#56C667FF","#FDE725FF"),
                        labels = c("Cov_Boot" = expression("Cov"["GE"]*": Bootstrap"),"Cov_Perm" = expression("Cov"["GE"]*": Permutation"),
                                   "GxE_Anova" = expression(bar(Delta)*""["GxE"]*": Anova"),"GxE_Perm" = expression(bar(Delta)*""["GxE"]*": Permutation")))+
    scale_fill_manual(values = c("Cov_Boot" = "#3CBB75FF","Cov_Perm" = "#3CBB75FF","GxE_Anova" = "#453781FF","GxE_Perm" = "#453781FF"),
                      labels = c("Cov_Boot" = expression("Cov"["GE"]*": Bootstrap"),"Cov_Perm" = expression("Cov"["GE"]*": Permutation"),
                                 "GxE_Anova" = expression(bar(Delta)*""["GxE"]*": Anova"),"GxE_Perm" = expression(bar(Delta)*""["GxE"]*": Permutation")))+
    scale_shape_manual(values = c("Cov_Boot" =21, "Cov_Perm" = 22,"GxE_Anova"=23, "GxE_Perm"=24),
                       labels = c("Cov_Boot" = expression("Cov"["GE"]*": Bootstrap"),"Cov_Perm" = expression("Cov"["GE"]*": Permutation"),
                                  "GxE_Anova" = expression(bar(Delta)*""["GxE"]*": Anova"),"GxE_Perm" = expression(bar(Delta)*""["GxE"]*": Permutation")))+
    scale_linetype_manual(values = c("Cov_Boot" = "dashed", "Cov_Perm" = "solid","GxE_Anova"= "twodash", "GxE_Perm"= "solid"),
                          labels = c("Cov_Boot" = expression("Cov"["GE"]*": Bootstrap"),"Cov_Perm" = expression("Cov"["GE"]*": Permutation"),
                                     "GxE_Anova" = expression(bar(Delta)*""["GxE"]*": Anova"),"GxE_Perm" = expression(bar(Delta)*""["GxE"]*": Permutation")))+
    xlab("Effect Size")+ylab("Power")+ theme(legend.position = "none")+
    ggtitle("C    FRT: Power according to effect size")+
    theme(plot.title = element_text(size=18))+
    theme(axis.text = element_text(colour = "black")))
   #theme(legend.text.align = 0)+
   #labs(colour = " ", fill = " ", shape = " ", linetype = " "))

###### Reciprocal Transplant - CovGE & GxE - Power ~ Total Sample Size #####
trade2_frt$linegrp2 = paste(trade2_frt$ID,sep = "_")
(Power_FRT = ggplot(trade2_frt %>% filter(between(bin,0.2,0.5)),aes(x = as.factor(totals), y = power)) + 
    #geom_rect(aes(xmin="0.2", xmax="0.5", ymin=-0.05, ymax=1.00), alpha=0.002) +
    geom_smooth(aes(group = linegrp2, colour = factor(ID),linetype = factor(ID)), se = F)+
    geom_point(aes(shape = factor(ID), group = grp, fill = factor(ID)),position = position_dodge(width = 0.4),size = 4,alpha = 0.7)+
    geom_hline(aes(yintercept = 0.8),linetype = "dashed")+
    theme_classic(base_size = 22, base_family = "Times")+
    scale_colour_manual(values = c("Cov_Boot" = "#3CBB75FF","Cov_Perm" = "#3CBB75FF",
                                   "GxE_Anova" = "#453781FF","GxE_Perm" = "#453781FF"),#"#56C667FF","#FDE725FF"),
                        labels = c("Cov_Boot" = expression("Cov"["GE"]*": Bootstrap"),"Cov_Perm" = expression("Cov"["GE"]*": Permutation"),
                                   "GxE_Anova" = expression(bar(Delta)*""["GxE"]*": Anova"),"GxE_Perm" = expression(bar(Delta)*""["GxE"]*": Permutation")))+
    scale_fill_manual(values = c("Cov_Boot" = "#3CBB75FF","Cov_Perm" = "#3CBB75FF","GxE_Anova" = "#453781FF","GxE_Perm" = "#453781FF"),
                      labels = c("Cov_Boot" = expression("Cov"["GE"]*": Bootstrap"),"Cov_Perm" = expression("Cov"["GE"]*": Permutation"),
                                 "GxE_Anova" = expression(bar(Delta)*""["GxE"]*": Anova"),"GxE_Perm" = expression(bar(Delta)*""["GxE"]*": Permutation")))+
    scale_shape_manual(values = c("Cov_Boot" =21, "Cov_Perm" = 22,"GxE_Anova"=23, "GxE_Perm"=24),
                       labels = c("Cov_Boot" = expression("Cov"["GE"]*": Bootstrap"),"Cov_Perm" = expression("Cov"["GE"]*": Permutation"),
                                  "GxE_Anova" = expression(bar(Delta)*""["GxE"]*": Anova"),"GxE_Perm" = expression(bar(Delta)*""["GxE"]*": Permutation")))+
    scale_linetype_manual(values = c("Cov_Boot" = "dashed", "Cov_Perm" = "solid","GxE_Anova"= "twodash", "GxE_Perm"= "solid"),
                          labels = c("Cov_Boot" = expression("Cov"["GE"]*": Bootstrap"),"Cov_Perm" = expression("Cov"["GE"]*": Permutation"),
                                     "GxE_Anova" = expression(bar(Delta)*""["GxE"]*": Anova"),"GxE_Perm" = expression(bar(Delta)*""["GxE"]*": Permutation")))+
    xlab("Total Sample Size")+ylab("Power")+ theme(legend.position = "none")+
    ggtitle("A    FRT: Power according to total sample size")+
    theme(plot.title = element_text(size=18))+
    theme(axis.text = element_text(colour = "black")))

###### Common Garden - CovGE & GxE - Power ~ Effect Size #####
sigGxE1 = fnr.effsize(dat_dub1, metric = "GxE", scenario = 2,  data.type = "raw", analysis ="perm",resolution = "fine")
sigGxE1$ID = rep("GxE_Perm", nrow(sigGxE1))
sigGxE2 = fnr.effsize(dat_dub1, metric = "GxE", scenario = 2,  data.type = "raw", analysis ="anova",resolution = "fine")
sigGxE2$ID = rep("GxE_Anova", nrow(sigGxE2))
sigcov1 = fnr.effsize(dat_dub1, metric = "Cov", scenario = 2,   data.type = "raw", analysis ="perm",resolution = "fine")
sigcov1$ID = rep("Cov_Perm", nrow(sigcov1))
sigcov2 = fnr.effsize(dat_dub1, metric = "Cov", scenario = 2,   data.type = "raw", analysis ="boot",resolution = "fine")
sigcov2$ID = rep("Cov_Boot", nrow(sigcov2))

trade2 = rbind(sigGxE1,sigGxE2, sigcov1,sigcov2)
trade2$totals = trade2$sample_size*2*trade2$n_pop
#trade3 = trade2 %>%
#  group_by(ID,totals,bin)%>%
#  summarize("avg_power" = mean(power))
trade2$grp = paste(trade2$bin,trade2$totals,trade2$ID,sep = "_")
trade2$linegrp = paste(trade2$totals,trade2$ID,sep = "_")

(Power256_CG = ggplot(trade2 %>% filter(totals == 256) %>% filter(bin != 0.0),aes(x = bin, y = power)) + 
    geom_rect(aes(xmin="0.2", xmax="0.5", ymin=0.00, ymax=1.00), alpha=0.002) +
    geom_smooth(aes(group = linegrp, colour = factor(ID),linetype = factor(ID)), se = F)+
    geom_point(aes(shape = factor(ID), group = grp, fill = factor(ID)),position = position_dodge(width = 0.4),size = 4,alpha = 0.7)+
    geom_hline(aes(yintercept = 0.8),linetype = "dashed")+
    theme_classic(base_size = 22, base_family = "Times")+
    scale_colour_manual(values = c("Cov_Boot" = "#3CBB75FF","Cov_Perm" = "#3CBB75FF",
                                   "GxE_Anova" = "#453781FF","GxE_Perm" = "#453781FF"),#"#56C667FF","#FDE725FF"),
                        labels = c("Cov_Boot" = expression("Cov"["GE"]*": Bootstrap"),"Cov_Perm" = expression("Cov"["GE"]*": Permutation"),
                                 "GxE_Anova" = expression(bar(Delta)*""["GxE"]*": Anova"),"GxE_Perm" = expression(bar(Delta)*""["GxE"]*": Permutation")))+
    scale_fill_manual(values = c("Cov_Boot" = "#3CBB75FF","Cov_Perm" = "#3CBB75FF","GxE_Anova" = "#453781FF","GxE_Perm" = "#453781FF"),
                      labels = c("Cov_Boot" = expression("Cov"["GE"]*": Bootstrap"),"Cov_Perm" = expression("Cov"["GE"]*": Permutation"),
                                 "GxE_Anova" = expression(bar(Delta)*""["GxE"]*": Anova"),"GxE_Perm" = expression(bar(Delta)*""["GxE"]*": Permutation")))+
    scale_shape_manual(values = c("Cov_Boot" =21, "Cov_Perm" = 22,"GxE_Anova"=23, "GxE_Perm"=24),
                       labels = c("Cov_Boot" = expression("Cov"["GE"]*": Bootstrap"),"Cov_Perm" = expression("Cov"["GE"]*": Permutation"),
                                 "GxE_Anova" = expression(bar(Delta)*""["GxE"]*": Anova"),"GxE_Perm" = expression(bar(Delta)*""["GxE"]*": Permutation")))+
    scale_linetype_manual(values = c("Cov_Boot" = "dashed", "Cov_Perm" = "solid","GxE_Anova"= "twodash", "GxE_Perm"= "solid"),
                          labels = c("Cov_Boot" = expression("Cov"["GE"]*": Bootstrap"),"Cov_Perm" = expression("Cov"["GE"]*": Permutation"),
                                     "GxE_Anova" = expression(bar(Delta)*""["GxE"]*": Anova"),"GxE_Perm" = expression(bar(Delta)*""["GxE"]*": Permutation")))+
    xlab("Effect Size")+ylab("Power")+
    theme(axis.text = element_text(colour = "black"))+ 
    ggtitle("D    CG: Power according to effect size")+
    theme(plot.title = element_text(size=18))+
    theme(legend.position = "none"))
    #theme(legend.text.align = 0)+
    #labs(colour = " ", fill = " ", shape = " ", linetype = " ")+
    #facet_wrap(~totals,ncol = 4))

###### Reciprocal Transplant - CovGE & GxE - Power ~ Total Sample Size #####

trade2$linegrp2 = paste(trade2$ID,sep = "_")

(Power_CG = ggplot(trade2 %>% filter(between(bin,0.2,0.5)),aes(x = as.factor(totals), y = power)) + 
    #geom_rect(aes(xmin="0.2", xmax="0.5", ymin=-0.05, ymax=1.00), alpha=0.002) +
    geom_point(aes(shape = factor(ID), group = grp, fill = factor(ID)),position = position_dodge(width = 0.4),size = 4,alpha = 0.7)+
    geom_smooth(aes(group = linegrp2, colour = factor(ID),linetype = factor(ID)), se = F)+
    geom_hline(aes(yintercept = 0.8),linetype = "dashed")+
    theme_classic(base_size = 22, base_family = "Times")+
    scale_colour_manual(values = c("Cov_Boot" = "#3CBB75FF","Cov_Perm" = "#3CBB75FF",
                                   "GxE_Anova" = "#453781FF","GxE_Perm" = "#453781FF"),#"#56C667FF","#FDE725FF"),
                        labels = c("Cov_Boot" = expression("Cov"["GE"]*": Bootstrap"),"Cov_Perm" = expression("Cov"["GE"]*": Permutation"),
                                   "GxE_Anova" = expression(bar(Delta)*""["GxE"]*": Anova"),"GxE_Perm" = expression(bar(Delta)*""["GxE"]*": Permutation")))+
    scale_fill_manual(values = c("Cov_Boot" = "#3CBB75FF","Cov_Perm" = "#3CBB75FF","GxE_Anova" = "#453781FF","GxE_Perm" = "#453781FF"),
                      labels = c("Cov_Boot" = expression("Cov"["GE"]*": Bootstrap"),"Cov_Perm" = expression("Cov"["GE"]*": Permutation"),
                                 "GxE_Anova" = expression(bar(Delta)*""["GxE"]*": Anova"),"GxE_Perm" = expression(bar(Delta)*""["GxE"]*": Permutation")))+
    scale_shape_manual(values = c("Cov_Boot" =21, "Cov_Perm" = 22,"GxE_Anova"=23, "GxE_Perm"=24),
                       labels = c("Cov_Boot" = expression("Cov"["GE"]*": Bootstrap"),"Cov_Perm" = expression("Cov"["GE"]*": Permutation"),
                                  "GxE_Anova" = expression(bar(Delta)*""["GxE"]*": Anova"),"GxE_Perm" = expression(bar(Delta)*""["GxE"]*": Permutation")))+
    scale_linetype_manual(values = c("Cov_Boot" = "dashed", "Cov_Perm" = "solid","GxE_Anova"= "twodash", "GxE_Perm"= "solid"),
                          labels = c("Cov_Boot" = expression("Cov"["GE"]*": Bootstrap"),"Cov_Perm" = expression("Cov"["GE"]*": Permutation"),
                                     "GxE_Anova" = expression(bar(Delta)*""["GxE"]*": Anova"),"GxE_Perm" = expression(bar(Delta)*""["GxE"]*": Permutation")))+
    xlab("Total Sample Size")+ylab("Power")+ theme(legend.position = "none")+
    ggtitle("B    CG: Power according to total sample size")+
    theme(plot.title = element_text(size=18))+
    theme(axis.text = element_text(colour = "black")))


###### Arrangement - Power Analysis Panels  #####

grid.arrange(Power_FRT,Power_CG,Power256_FRT, Power256_CG)


##################################
###  Confusion Matrix on Means  ###
###################################

# MEANS -- CovGE -- Permutation
dat_csv1$meansCovconfintperm = rep("NA",nrow(dat_csv1))
for(i in 1:nrow(dat_csv1)){
  if(dat_csv1$true_cov_means[i] != 0 && dat_csv1$cov_means_pvalue[i] <= 0.05){dat_csv1$meansCovconfintperm[i] = "True Positive"
  }else if(dat_csv1$true_cov_means[i] == 0 & dat_csv1$cov_means_pvalue[i] <= 0.05){dat_csv1$meansCovconfintperm[i] = "False Positive"
  }else if(dat_csv1$true_cov_means[i]!= 0 & dat_csv1$cov_means_pvalue[i] > 0.05){dat_csv1$meansCovconfintperm[i] = "False Negative"
  }else if(dat_csv1$true_cov_means[i] == 0 & dat_csv1$cov_means_pvalue[i] > 0.05){dat_csv1$meansCovconfintperm[i] = "True Negative"
  }else{dat_csv1$meansCovconfintperm[i] = "None"}
}

# MEANS CovGE -- Bootstrap
dat_csv1$MeansCovconfintboot = rep(NA, nrow(dat_csv1))

for(i in 1:nrow(dat_csv1)){
  
  if(dat_csv1$true_cov_means[i] != 0 &&
     dat_csv1$cov_means_lwrCI[i] < 0 &&
     dat_csv1$cov_means_uprCI[i] < 0
  ){dat_csv1$MeansCovconfintboot[i] = "True Positive"
  }else if(dat_csv1$true_cov_means[i] != 0 &&
           dat_csv1$cov_means_lwrCI[i] > 0 &&
           dat_csv1$cov_means_uprCI[i] > 0
  ){dat_csv1$MeansCovconfintboot[i] = "True Positive"
  }else if(dat_csv1$true_cov_means[i] == 0 &&
           dat_csv1$cov_means_lwrCI[i] < 0 &&
           dat_csv1$cov_means_uprCI[i] < 0
  ){dat_csv1$MeansCovconfintboot[i] = "False Positive"
  }else if(dat_csv1$true_cov_means[i] == 0 &&
           dat_csv1$cov_means_lwrCI[i] > 0 &&
           dat_csv1$cov_means_uprCI[i] > 0
  ){dat_csv1$MeansCovconfintboot[i] = "False Positive"
  }else if(dat_csv1$true_cov_means[i] != 0 && 
           dat_csv1$cov_means_lwrCI[i] <= 0 && 
           dat_csv1$cov_means_uprCI[i] >= 0
  ){dat_csv1$MeansCovconfintboot[i] = "False Negative"
  }else if(dat_csv1$true_cov_means[i]== 0 && 
           dat_csv1$cov_means_lwrCI[i] <= 0 && 
           dat_csv1$cov_means_uprCI[i] >= 0
  ){dat_csv1$MeansCovconfintboot[i] = "True Negative"
  }else{dat_csv1$MeansCovconfintboot[i] = "None"}
}

# MEANS GxE -- Bootstrap

dat_csv1$MeanGxEconfintboot = rep("NA", nrow(dat_csv1))
for(i in 1:nrow(dat_csv1)){
  if(dat_csv1$true_GxE_means[i] != 0 & round(dat_csv1$GxE_means_lwrCI[i],2) > 0)
  {dat_csv1$MeanGxEconfintboot[i] = "True Positive"
  }else if(dat_csv1$true_GxE_means[i] == 0 & round(dat_csv1$GxE_means_lwrCI[i],2) > 0)
  {dat_csv1$MeanGxEconfintboot[i] = "False Positive"
  }else if(dat_csv1$true_GxE_means[i] != 0 & round(dat_csv1$GxE_means_lwrCI[i],2) <= 0)
  {dat_csv1$MeanGxEconfintboot[i] = "False Negative"
  }else if(dat_csv1$true_GxE_means[i] == 0 & round(dat_csv1$GxE_means_lwrCI[i],2) <= 0)
  {dat_csv1$MeanGxEconfintboot[i] = "True Negative"
  }else{dat_csv1$MeanGxEconfintboot[i] = "None"}
}

# MEANS GxE -- Permutation
dat_csv1$meansGxEconfintperm = rep("NA", nrow(dat_csv1))

for(i in 1:nrow(dat_csv1)){
  if(dat_csv1$true_GxE_means[i] != 0 & dat_csv1$GxE_means_pvalue[i] <= 0.05){dat_csv1$meansGxEconfintperm[i] = "True Positive"
  }else if(dat_csv1$true_GxE_means[i] == 0 & dat_csv1$GxE_means_pvalue[i] <= 0.05){dat_csv1$meansGxEconfintperm[i] = "False Positive"
  }else if(dat_csv1$true_GxE_means[i] != 0 & dat_csv1$GxE_means_pvalue[i] > 0.05){dat_csv1$meansGxEconfintperm[i] = "False Negative"
  }else if(dat_csv1$true_GxE_means[i] == 0 & dat_csv1$GxE_means_pvalue[i] > 0.05){dat_csv1$meansGxEconfintperm[i] = "True Negative"
  }else{dat_csv1$meansGxEconfintperm == "None"}
}


####### Reciprocal Transplant - Overall FPR/FNR  - Means ######

means_cov_perm_table = dat_csv1 %>%
  group_by("name" = meansCovconfintperm) %>%
  summarize("n" = n())
fpr.fnr(means_cov_perm_table, divided = FALSE, scenario = 1)

means_cov_boot_table = dat_csv1 %>%
  group_by("name" =MeansCovconfintboot) %>%
  summarize("n" = n())
fpr.fnr(means_cov_boot_table, divided = FALSE, scenario = 1)

means_gxe_boot_table = dat_csv1 %>%
  group_by("name" =MeanGxEconfintboot) %>%
  summarize("n" = n())
fpr.fnr(means_gxe_boot_table, divided = FALSE, scenario = 1)

means_gxe_perm_table = dat_csv1 %>%
  group_by("name" =meansGxEconfintperm) %>%
  summarize("n" = n())
fpr.fnr(means_gxe_perm_table, divided = FALSE, scenario = 1)

####### Reciprocal Transplant - False Positive Tables  - Means ######
(means_confusion_hmap1 = dat_csv1 %>%
    group_by(sample_size, n_pop, "name" =meansCovconfintperm) %>%
    summarize("n" = n()))
means_conf1 = fpr.fnr(means_confusion_hmap1, divided = TRUE, scenario = 1)
#means_conf_plot1 <- heatmap_fun(means_conf1,"rate") # Can also do "percent"

(means_confusion_hmap2 = dat_csv1 %>%
    group_by(sample_size, n_pop, "name" =MeansCovconfintboot) %>%
    summarize("n" = n()))
means_conf2 = fpr.fnr(means_confusion_hmap2, divided = TRUE, scenario = 1)
#means_conf_plot2 <- heatmap_fun(means_conf2,"rate")

(means_confusion_hmap3 = dat_csv1 %>%
    group_by(sample_size, n_pop,"name" = MeanGxEconfintboot) %>%
    summarize("n" = n()))
means_conf3 = fpr.fnr(means_confusion_hmap3, divided = TRUE, scenario = 1)
#means_conf_plot3 <- heatmap_fun(means_conf3,"rate")

means_confusion_hmap4 = dat_csv1 %>%
  group_by(sample_size, n_pop,"name" = meansGxEconfintperm) %>%
  summarize("n" = n())
means_conf4 = fpr.fnr(means_confusion_hmap4, divided = TRUE, scenario = 1)
#means_conf_plot4 <- heatmap_fun(means_conf4,"rate")

# Compile FPs for plot
means_conf1$ID = rep("Cov_Perm", nrow(means_conf1))
means_conf2$ID = rep("Cov_Boot", nrow(means_conf2))
means_conf3$ID = rep("GxE_Boot", nrow(means_conf3))
means_conf4$ID = rep("GxE_Perm", nrow(means_conf4))
fpdf_means = rbind(means_conf1,means_conf2,means_conf3,means_conf4)
fpdf_means = fpdf_means[fpdf_means$name == "False Positive",]

## False Positive Rates
fpdf_means$npop_plot = NA
for(i in 1:nrow(fpdf_means)){
  if(fpdf_means$n_pop[i] == 2){fpdf_means$npop_plot[i] = "2 Populations"
  }else if(fpdf_means$n_pop[i] == 4){fpdf_means$npop_plot[i] = "4 Populations"
  }else{fpdf_means$npop_plot[i] = "8 Populations"}
}

gxeFPRmean = rbind(means_conf3,means_conf4)
gxeFPRmean = gxeFPRmean[gxeFPRmean$name == "False Positive",]

covFPRmean = rbind(means_conf1,means_conf2)
covFPRmean = covFPRmean[covFPRmean$name == "False Positive",]


####### Reciprocal Transplant - False Negative Tables  - Means ######
covperm1mean = fnr.effsize(dat_csv1, metric = "cov", data.type = "means", analysis ="perm",resolution = "fine",scenario = 1)
covperm1mean$ID = rep("Cov_Perm",nrow(covperm1mean))
covboot1mean = fnr.effsize(dat_csv1, metric = "cov", data.type = "means",  analysis ="boot",resolution = "fine")
covboot1mean$ID = rep("Cov_Boot",nrow(covboot1mean))
gxeperm1mean = fnr.effsize(dat_csv1, metric = "gxe", data.type = "means",  analysis ="perm",resolution = "fine")
gxeperm1mean$ID = rep("GxE_Perm",nrow(gxeperm1mean))
gxeboot1mean = fnr.effsize(dat_csv1, metric = "gxe", data.type = "means", analysis = "boot",resolution = "fine")
gxeboot1mean$ID = rep("GxE_Boot",nrow(gxeboot1mean))

fndfmean = rbind(covperm1mean,covboot1mean,gxeperm1mean,gxeboot1mean)
gxePowmean = rbind(gxeperm1mean,gxeboot1mean)
covPowmean = rbind(covperm1mean,covboot1mean)


###############################################
###  Confusion Matrix on Means - Env Scen 2 ###
###############################################

dat_dub$meansCovconfintperm = rep("NA",nrow(dat_dub))

for(i in 1:nrow(dat_dub)){
  if(dat_dub$true_cov_means[i] != 0 && dat_dub$cov_means_pvalue[i] <= 0.05){dat_dub$meansCovconfintperm[i] = "True Positive"
  }else if(dat_dub$true_cov_means[i] == 0 & dat_dub$cov_means_pvalue[i] <= 0.05){dat_dub$meansCovconfintperm[i] = "False Positive"
  }else if(dat_dub$true_cov_means[i]!= 0 & dat_dub$cov_means_pvalue[i] > 0.05){dat_dub$meansCovconfintperm[i] = "False Negative"
  }else if(dat_dub$true_cov_means[i] == 0 & dat_dub$cov_means_pvalue[i] > 0.05){dat_dub$meansCovconfintperm[i] = "True Negative"
  }else{dat_dub$meansCovconfintperm[i] = "None"}
}

# MEANS Cov Boot check
dat_dub$MeansCovconfintboot = rep(NA, nrow(dat_dub))

for(i in 1:nrow(dat_dub)){
  
  if(dat_dub$true_cov_means[i] != 0 &&
     dat_dub$cov_means_lwrCI[i] < 0 &&
     dat_dub$cov_means_uprCI[i] < 0
  ){dat_dub$MeansCovconfintboot[i] = "True Positive"
  }else if(dat_dub$true_cov_means[i] != 0 &&
           dat_dub$cov_means_lwrCI[i] > 0 &&
           dat_dub$cov_means_uprCI[i] > 0
  ){dat_dub$MeansCovconfintboot[i] = "True Positive"
  }else if(dat_dub$true_cov_means[i] == 0 &&
           dat_dub$cov_means_lwrCI[i] < 0 &&
           dat_dub$cov_means_uprCI[i] < 0
  ){dat_dub$MeansCovconfintboot[i] = "False Positive"
  }else if(dat_dub$true_cov_means[i] == 0 &&
           dat_dub$cov_means_lwrCI[i] > 0 &&
           dat_dub$cov_means_uprCI[i] > 0
  ){dat_dub$MeansCovconfintboot[i] = "False Positive"
  }else if(dat_dub$true_cov_means[i] != 0 && 
           dat_dub$cov_means_lwrCI[i] <= 0 && 
           dat_dub$cov_means_uprCI[i] >= 0
  ){dat_dub$MeansCovconfintboot[i] = "False Negative"
  }else if(dat_dub$true_cov_means[i]== 0 && 
           dat_dub$cov_means_lwrCI[i] <= 0 && 
           dat_dub$cov_means_uprCI[i] >= 0
  ){dat_dub$MeansCovconfintboot[i] = "True Negative"
  }else{dat_dub$MeansCovconfintboot[i] = "None"}
}

# MEANS -- GxE --  Bootstrap - Scen 2
dat_dub$MeanGxEconfintboot = rep("NA", nrow(dat_dub))

for(i in 1:nrow(dat_dub)){
  if(dat_dub$true_GxE_means[i] != 0 & dat_dub$GxE_means_lwrCI[i] > 0)
  {dat_dub$MeanGxEconfintboot[i] = "True Positive"
  }else if(dat_dub$true_GxE_means[i] == 0 & dat_dub$GxE_means_lwrCI[i] > 0)
  {dat_dub$MeanGxEconfintboot[i] = "False Positive"
  }else if(dat_dub$true_GxE_means[i] != 0 & dat_dub$GxE_means_lwrCI[i] <= 0)
  {dat_dub$MeanGxEconfintboot[i] = "False Negative"
  }else if(dat_dub$true_GxE_means[i] == 0 & dat_dub$GxE_means_lwrCI[i] <= 0)
  {dat_dub$MeanGxEconfintboot[i] = "True Negative"
  }else{dat_dub$MeanGxEconfintboot[i] = "None"}
}

# MEANS GxE Perm check
dat_dub$meansGxEconfintperm = rep("NA", nrow(dat_dub))

for(i in 1:nrow(dat_dub)){
  if(dat_dub$true_GxE_means[i] != 0 & dat_dub$GxE_means_pvalue[i] <= 0.05){dat_dub$meansGxEconfintperm[i] = "True Positive"
  }else if(dat_dub$true_GxE_means[i] == 0 & dat_dub$GxE_means_pvalue[i] <= 0.05){dat_dub$meansGxEconfintperm[i] = "False Positive"
  }else if(dat_dub$true_GxE_means[i] != 0 & dat_dub$GxE_means_pvalue[i] > 0.05){dat_dub$meansGxEconfintperm[i] = "False Negative"
  }else if(dat_dub$true_GxE_means[i] == 0 & dat_dub$GxE_means_pvalue[i] > 0.05){dat_dub$meansGxEconfintperm[i] = "True Negative"
  }else{dat_dub$meansGxEconfintperm == "None"}
}

## Counts for table
dat_dub1 <- dat_dub %>% 
  filter(std_dev == 1) %>%
  #filter(between(abs(covariance),0.2,0.6))%>%
  filter(total_samples > 17)
#filter(total_samples > 128)
dat_dub2 <- dat_dub %>% 
  filter(std_dev == 1) %>%
  # filter(between(GxE_emm,0.3,0.6)) %>%
  filter(total_samples > 17)
#filter(total_samples > 128)
 

####### Common Garden - Overall FPR/FNR  - Means ########

means_cov_perm_table = dat_dub1 %>%
  group_by("name" = meansCovconfintperm) %>%
  summarize("n" = n())
fpr.fnr(means_cov_perm_table, divided = FALSE, scenario = 2)

means_cov_boot_table = dat_dub1 %>%
  group_by("name" =MeansCovconfintboot) %>%
  summarize("n" = n())
fpr.fnr(means_cov_boot_table, divided = FALSE, scenario = 2)

means_gxe_boot_table = dat_dub2 %>%
  group_by("name" =MeanGxEconfintboot) %>%
  summarize("n" = n())
fpr.fnr(means_gxe_boot_table, divided = FALSE, scenario = 2)

means_gxe_perm_table = dat_dub2 %>%
  group_by("name" =meansGxEconfintperm) %>%
  summarize("n" = n())
fpr.fnr(means_gxe_perm_table, divided = FALSE, scenario = 2)

####### Common Garden - False Positive Tables  - Means ########
(means_confusion_hmap1 = dat_dub1 %>%
    group_by(sample_size, n_pop, "name" =meansCovconfintperm) %>%
    summarize("n" = n()))
dub_means_conf1 = fpr.fnr(means_confusion_hmap1, divided = TRUE, scenario = 2)
#means_conf_plot1 <- heatmap_fun(dub_means_conf1,"rate") # Can also do "percent"

(means_confusion_hmap2 = dat_dub1 %>%
    group_by(sample_size, n_pop, "name" =MeansCovconfintboot) %>%
    summarize("n" = n()))
dub_means_conf2 = fpr.fnr(means_confusion_hmap2, divided = TRUE, scenario = 2)
#means_conf_plot2 <- heatmap_fun(dub_means_conf2,"rate")

(means_confusion_hmap3 = dat_dub2 %>%
    group_by(sample_size, n_pop,"name" = MeanGxEconfintboot) %>%
    summarize("n" = n()))
dub_means_conf3 = fpr.fnr(means_confusion_hmap3, divided = TRUE, scenario = 2)
#means_conf_plot3 <- heatmap_fun(dub_means_conf3,"rate")

means_confusion_hmap4 = dat_dub2 %>%
  group_by(sample_size, n_pop,"name" = meansGxEconfintperm) %>%
  summarize("n" = n())
dub_means_conf4 = fpr.fnr(raw_confusion_hmap4, divided = TRUE, scenario = 2)
#means_conf_plot4 <- heatmap_fun(dub_means_conf4,"rate")

# Compile FPs for plot
dub_means_conf1$ID = rep("Cov_Perm", nrow(dub_means_conf1))
dub_means_conf2$ID = rep("Cov_Boot", nrow(dub_means_conf2))
dub_means_conf3$ID = rep("GxE_Boot", nrow(dub_means_conf3))
dub_means_conf4$ID = rep("GxE_Perm", nrow(dub_means_conf4))
dub_fpdf_means = rbind(dub_means_conf1,dub_means_conf2,dub_means_conf3,dub_means_conf4)
dub_fpdf_means = dub_fpdf_means[dub_fpdf_means$name == "False Positive",]

## False Positive Rates
dub_fpdf_means$npop_plot = NA
for(i in 1:nrow(dub_fpdf_means)){
  if(dub_fpdf_means$n_pop[i] == 2){dub_fpdf_means$npop_plot[i] = "2 Populations"
  }else if(dub_fpdf_means$n_pop[i] == 4){dub_fpdf_means$npop_plot[i] = "4 Populations"
  }else{dub_fpdf_means$npop_plot[i] = "8 Populations"}
}

gxeFPRmeandub = rbind(dub_means_conf3,dub_means_conf4)
gxeFPRmeandub = gxeFPRmeandub[gxeFPRmeandub$name == "False Positive",]

covFPRmeandub = rbind(dub_means_conf1,dub_means_conf2)
covFPRmeandub = covFPRmeandub[covFPRmeandub$name == "False Positive",]


####### Common Garden - False Negative Tables  - Means ########
covperm1mean_dub = fnr.effsize(dat_dub1, metric = "cov", data.type = "means", analysis ="perm",scenario = 2,resolution = "fine")
covperm1mean_dub$ID = rep("Cov_Perm",nrow(covperm1mean_dub))
covboot1mean_dub = fnr.effsize(dat_dub1, metric = "cov", data.type = "means",  analysis ="boot",scenario = 2,resolution = "fine")
covboot1mean_dub$ID = rep("Cov_Boot",nrow(covboot1mean_dub))
gxeperm1mean_dub = fnr.effsize(dat_dub2, metric = "gxe", data.type = "means",  analysis ="perm",scenario = 2,resolution = "fine")
gxeperm1mean_dub$ID = rep("GxE_Perm",nrow(gxeperm1mean_dub))
gxeboot1mean_dub = fnr.effsize(dat_dub2, metric = "gxe", data.type = "means", analysis = "boot",scenario = 2,resolution = "fine")
gxeboot1mean_dub$ID = rep("GxE_Boot",nrow(gxeboot1mean_dub))

fndfmean_dub = rbind(covperm1mean_dub,covboot1mean_dub,gxeperm1mean_dub,gxeboot1mean_dub)
gxePowmean_dub = rbind(gxeperm1mean_dub,gxeboot1mean_dub)
covPowmean_dub = rbind(covperm1mean_dub,covboot1mean_dub)


####### Reciprocal Transplant - Means - CovGE - False Negative Heatmap ######
covPowmean[is.nan(covPowmean)] <- 0
covPowmean1 <- 
  covPowmean %>% 
  filter(between(bin,0.2,0.6))%>%
  group_by(sample_size,n_pop,ID) %>%
  summarize("FNR" = mean(fnr))
covPowmean1$n_pop<-as.factor(covPowmean1$n_pop)
covPowmean2<-
  covPowmean1 %>%
  mutate(#xadj = ifelse(ID %in% c("Cov_Perm"), -.45/2,  .45/2),
    yadj = ifelse(ID %in% c("Cov_Perm"), -.475/2,  .475/2),
    xpos = as.numeric(factor(sample_size)),
    ypos = as.numeric(n_pop) + yadj,
    ID2 = ifelse(ID %in% c("Cov_Perm"), "Perm.", "Boot."))

(covFNFRTmean = covPowmean2 %>% 
    ggplot(aes(as.factor(sample_size), ypos, fill = FNR)) + 
    geom_tile(height = 0.475, width = 0.95, color= "white") +
    geom_text(aes(label = paste0(ID2,"\n",round(FNR, 2))), size = 4, family = "Times", colour = "white", show.legend = F) +
    #ggtitle("False Negative Rates: CovGE")+
    scale_fill_viridis(breaks=seq(0,1,0.25), #breaks in the scale bar
                       limits=c(0,1))+
    scale_x_discrete(name = "Sample Size",
                     labels = unique(covPowmean2$sample_size)) +
    scale_y_continuous(breaks = 1:3, 
                       name = "Number of Genotypes",
                       labels = c(2,4,8)) +
    labs(fill = "False Negative Rate")+
    #  ggtitle(expression("Reciprocal Transplant: Cov"["GE"]*" False Negative Rates"))+
    theme_classic(base_size = 18, base_family = "Times")+
    theme(axis.text = element_text(colour = "black"))+
    theme(legend.position = "none")+
    ggtitle("FRT: False Negative Rates for Means")+
    theme(plot.title = element_text(size = 18, face = "bold")))

####### Reciprocal Transplant - Means - CovGE - False Postive Heatmap ######
covFPRmean[is.nan(covFPRmean)] <- 0
covFPRmean1 <- 
  covFPRmean %>% 
  mutate(#xadj = ifelse(ID %in% c("Cov_Perm"), -.45/2,  .45/2),
    yadj = ifelse(ID %in% c("Cov_Perm"), -.475/2,  .475/2),
    xpos = as.numeric(factor(sample_size)),
    ypos = as.numeric(factor(n_pop)) + yadj,
    ID2 = ifelse(ID %in% c("Cov_Perm"), "Perm.", "Boot."))

(covFPFRTmean = covFPRmean1 %>% 
    ggplot(aes(as.factor(sample_size), ypos, fill = rate)) + 
    geom_tile(height = 0.475, width = 0.95, color= "white") +
    geom_text(aes(label = paste0(ID2,"\n",round(rate, 2))), size = 4, family = "Times", colour = "white", show.legend = F) +
    #ggtitle("False Positive Rates: CovGE")+
    scale_fill_viridis(breaks=seq(0,0.4,0.1), #breaks in the scale bar
                       limits=c(0,.4))+
    scale_x_discrete(name = "Sample Size",
                     labels = unique(covFPRmean1$sample_size)) +
    scale_y_continuous(breaks = 1:3, name = "Number of Genotypes",
                       labels = c(2,4,8)) +
    labs(fill = "False Positive Rate")+
    #  ggtitle(expression("Reciprocal Transplant: Cov"["GE"]*" False Positive Rates"))+
    
    theme_classic(base_size = 18, base_family = "Times")+
    theme(axis.text = element_text(colour = "black"))+
    theme(legend.position = "none")+
    ggtitle("FRT: False Positive Rates for Means")+
    theme(plot.title = element_text(size = 18, face = "bold"))) 

####### Reciprocal Transplant - Means - GxE - False Positive Heatmap ######
gxeFPRmean[is.nan(gxeFPRmean)] <- 0
gxeFPRmean1 <- 
  gxeFPRmean %>% 
  mutate(
    yadj = ifelse(ID == "GxE_Perm", -.475/2,  .475/2),
    #ifelse(ID == "GxE_Perm", 0.0, 0.9/3)),
    xpos = as.numeric(factor(sample_size)),
    ypos = as.numeric(factor(n_pop))+yadj,
    col = ifelse(rate >= 0.5, "black","white"),
    ID2 = ifelse(ID == "GxE_Boot", "Boot.", "Perm."))

(gxeFPFRTmean = filter(gxeFPRmean1,ID2 == "Perm.") %>% 
    ggplot(aes(as.factor(sample_size), ypos, fill = rate)) + 
    geom_tile(height = .475, width = 0.95, color= "white") +
    geom_text(aes(label = paste0(ID2,"\n",round(rate, 2)), colour = col), 
              size = 4, family = "Times", show.legend = F) +
    scale_colour_identity()+
    scale_fill_viridis(breaks=seq(0,1,0.25), #breaks in the scale bar
                       limits=c(0,1))+
    scale_x_discrete(name = "Sample Size",
                     labels = unique(gxeFPRmean$sample_size)) +
    scale_y_continuous(breaks = 1:3, name = "Number of Genotypes",
                       labels = c(2,4,8)) +
    # ggtitle(expression("Reciprocal Transplant: "*bar(Delta)*""["GxE"]*" False Positive Rates"))+  
    labs(fill = "False Positive Rate")+
    theme_classic(base_size = 18, base_family = "Times")+
    theme(axis.text = element_text(colour = "black"))+
    theme(legend.position = "none")+
    ggtitle("FRT: False Positive Rates for Means")+
    theme(plot.title = element_text(size = 18, face = "bold")))

####### Reciprocal Transplant - Means - GxE - False Negative Heatmap ######
gxePowmean[is.nan(gxePowmean)] <- 0
gxePowmean1 = gxePowmean %>%
  filter(between(bin,0.3,0.6))%>%  
  group_by(sample_size,n_pop,ID) %>%
  summarize("fnr" = mean(fnr))
gxePowmean1$n_pop<-as.factor(gxePowmean1$n_pop)
gxePowmean2 <- 
  gxePowmean1 %>% 
  mutate(#xadj = ifelse(ID == "GxE_Anova", -.45/3,ifelse(ID == "GxE_Perm", 0.0, .45/3)),
    yadj = ifelse(ID == "GxE_Perm", -.475/2,  .475/2),
    # ifelse(ID == "GxE_Perm", 0.0, 0.9/3)),
    xpos = as.numeric(factor(sample_size)),
    ypos = as.numeric(n_pop)+yadj,
    col = ifelse(fnr >= 0.5, "black","white"),
    ID2 = ifelse(ID == "GxE_Perm", "Perm.", "Boot."))

(gxeFNFRTmean = filter(gxePowmean2, ID2 == "Perm.") %>% 
    ggplot(aes(as.factor(sample_size), ypos, fill = fnr)) + 
    geom_tile(height = 0.475, width = 0.95, color= "white") +
    geom_text(aes(label = paste0(ID2,"\n",round(fnr, 2)), colour = col), 
              size = 4, family = "Times", show.legend = F) +
    scale_colour_identity()+
    scale_fill_viridis(breaks=seq(0,1,0.25), #breaks in the scale bar
                       limits=c(0,1))+
    scale_x_discrete(name = "Sample Size",
                     labels = unique(gxePowmean2$sample_size)) +
    scale_y_continuous(breaks = 1:3, name = "Number of Genotypes",
                       labels = c(2,4,8)) +
    # ggtitle(expression("Reciprocal Transplant: "*bar(Delta)*""["GxE"]*" False Positive Rates"))+  
    labs(fill = "False Negative Rate")+
    theme_classic(base_size = 18, base_family = "Times")+
    theme(axis.text = element_text(colour = "black"))+
    theme(legend.position = "none")+
    ggtitle("FRT: False Negative Rates for Means")+
    theme(plot.title = element_text(size = 18, face = "bold")))




####### Common Garden - Means - CovGE - False Negative Heatmap #####
covPowmean_dub[is.nan(covPowmean_dub)] <- 0
covPowmean_dub1 <- 
  covPowmean_dub %>% 
  filter(between(bin,0.2,0.6))%>%
  group_by(sample_size,n_pop,ID) %>%
  summarize("FNR" = mean(fnr))
covPowmean_dub1$n_pop<-as.factor(covPowmean_dub1$n_pop)
covPowmean_dub2<-
  covPowmean_dub1 %>%
  mutate(#xadj = ifelse(ID %in% c("Cov_Perm"), -.45/2,  .45/2),
    yadj = ifelse(ID %in% c("Cov_Perm"), -.475/2,  .475/2),
    xpos = as.numeric(factor(sample_size)),
    ypos = as.numeric(n_pop) + yadj,
    ID2 = ifelse(ID %in% c("Cov_Perm"), "Perm.", "Boot."))

(covFNFRTmean_dub = covPowmean_dub2 %>% 
    ggplot(aes(as.factor(sample_size), ypos, fill = FNR)) + 
    geom_tile(height = 0.475, width = 0.95, color= "white") +
    geom_text(aes(label = paste0(ID2,"\n",round(FNR, 2))), size = 4, family = "Times", colour = "white", show.legend = F) +
    #ggtitle("False Negative Rates: CovGE")+
    scale_fill_viridis(breaks=seq(0,1,0.25), #breaks in the scale bar
                       limits=c(0,1))+
    scale_x_discrete(name = "Sample Size",
                     labels = unique(covPowmean_dub2$sample_size)) +
    scale_y_continuous(breaks = 1:3, 
                       name = "Number of Genotypes",
                       labels = c(2,4,8)) +
    labs(fill = "False Negative Rate")+
    #  ggtitle(expression("Reciprocal Transplant: Cov"["GE"]*" False Negative Rates"))+
    theme_classic(base_size = 18, base_family = "Times")+
    theme(axis.text = element_text(colour = "black"))+
    theme(legend.position = "none")+
    ggtitle("CG: False Negative Rates for Means")+
    theme(plot.title = element_text(size = 18, face = "bold")))

####### Common Garden - Means - CovGE - False Positive Heatmap #####
covFPRmeandub[is.nan(covFPRmeandub)] <- 0
covFPRmeandub1 <- 
  covFPRmeandub %>% 
  mutate(#xadj = ifelse(ID %in% c("Cov_Perm"), -.45/2,  .45/2),
    yadj = ifelse(ID %in% c("Cov_Perm"), -.475/2,  .475/2),
    xpos = as.numeric(factor(sample_size)),
    ypos = as.numeric(factor(n_pop)) + yadj,
    ID2 = ifelse(ID %in% c("Cov_Perm"), "Perm.", "Boot."))


(covFPRmean_dub = covFPRmeandub1 %>% 
    ggplot(aes(as.factor(sample_size), ypos, fill = rate)) + 
    geom_tile(height = 0.475, width = 0.95, color= "white") +
    geom_text(aes(label = paste0(ID2,"\n",round(rate, 2))), size = 4, family = "Times", colour = "white", show.legend = F) +
    #ggtitle("False Positive Rates: CovGE")+
    scale_fill_viridis(breaks=seq(0,1,0.25), #breaks in the scale bar
                       limits=c(0,1))+
    scale_x_discrete(name = "Sample Size",
                     labels = unique(covFPRmeandub$sample_size)) +
    scale_y_continuous(breaks = 1:3, name = "Number of Genotypes",
                       labels = c(2,4,8)) +
    labs(fill = "False Positive Rate")+
    #  ggtitle(expression("Reciprocal Transplant: Cov"["GE"]*" False Positive Rates"))+
    
    theme_classic(base_size = 18, base_family = "Times")+
    theme(axis.text = element_text(colour = "black"))+
    theme(legend.position = "none")+
    ggtitle("CG: False Positive Rates for Means")+
    theme(plot.title = element_text(size = 18, face = "bold"))) 

####### Common Garden - Means - GxE - False Positive Heatmap #####
gxeFPRmeandub[is.nan(gxeFPRmeandub)] <- 0
gxeFPRmeandub1 <- 
  gxeFPRmeandub %>% 
  mutate(
    yadj = ifelse(ID == "GxE_Perm", -.475/2,  .475/2),
    #ifelse(ID == "GxE_Perm", 0.0, 0.9/3)),
    xpos = as.numeric(factor(sample_size)),
    ypos = as.numeric(factor(n_pop))+yadj,
    col = ifelse(rate >= 0.5, "black","white"),
    ID2 = ifelse(ID == "GxE_Boot", "Boot.", "Perm."))

(gxeFPCGmean_dub = filter(gxeFPRmeandub1, ID2 == "Perm.") %>% 
    ggplot(aes(as.factor(sample_size), ypos, fill = rate)) + 
    geom_tile(height = .475, width = 0.95, color= "white") +
    geom_text(aes(label = paste0(ID2,"\n",round(rate, 2)), colour = col), 
              size = 4, family = "Times", show.legend = F) +
    scale_colour_identity()+
    scale_fill_viridis(breaks=seq(0,1,0.25), #breaks in the scale bar
                       limits=c(0,1))+
    scale_x_discrete(name = "Sample Size",
                     labels = unique(gxeFPRmeandub1$sample_size)) +
    scale_y_continuous(breaks = 1:3, name = "Number of Genotypes",
                       labels = c(2,4,8)) +
    # ggtitle(expression("Reciprocal Transplant: "*bar(Delta)*""["GxE"]*" False Positive Rates"))+  
    labs(fill = "False Positive Rate")+
    theme_classic(base_size = 18, base_family = "Times")+
    theme(axis.text = element_text(colour = "black"))+
    theme(legend.position = "none")+
    ggtitle("CG: False Positive Rates for Means")+
    theme(plot.title = element_text(size = 18, face = "bold")))

####### Common Garden - Means - GxE - False Negative Heatmap #####
gxePowmean_dub[is.nan(gxePowmean_dub)] <- 0
gxePowmean_dub1 = gxePowmean_dub %>%
  filter(between(bin,0.3,0.6))%>%  
  
  group_by(sample_size,n_pop,ID) %>%
  summarize("fnr" = mean(fnr))
gxePowmean_dub1$n_pop<-as.factor(gxePowmean_dub1$n_pop)
gxePowmean_dub2 <- 
  gxePowmean_dub1 %>% 
  mutate(#xadj = ifelse(ID == "GxE_Anova", -.45/3,ifelse(ID == "GxE_Perm", 0.0, .45/3)),
    yadj = ifelse(ID == "GxE_Perm", -.475/2,  .475/2),
    # ifelse(ID == "GxE_Perm", 0.0, 0.9/3)),
    xpos = as.numeric(factor(sample_size)),
    ypos = as.numeric(n_pop)+yadj,
    col = ifelse(fnr >= 0.5, "black","white"),
    ID2 = ifelse(ID == "GxE_Perm", "Perm.", "Boot."))

(gxeFNCGmean_dub = filter(gxePowmean_dub2, ID2 == "Perm.") %>% 
    ggplot(aes(as.factor(sample_size), ypos, fill = fnr)) + 
    geom_tile(height = 0.475, width = 0.95, color= "white") +
    geom_text(aes(label = paste0(ID2,"\n",round(fnr, 2)), colour = col), 
              size = 4, family = "Times", show.legend = F) +
    scale_colour_identity()+
    scale_fill_viridis(breaks=seq(0,1,0.25), #breaks in the scale bar
                       limits=c(0,1))+
    scale_x_discrete(name = "Sample Size",
                     labels = unique(gxePowmean_dub2$sample_size)) +
    scale_y_continuous(breaks = 1:3, name = "Number of Genotypes",
                       labels = c(2,4,8)) +
    # ggtitle(expression("Reciprocal Transplant: "*bar(Delta)*""["GxE"]*" False Positive Rates"))+  
    labs(fill = "False Negative Rate")+
    theme_classic(base_size = 18, base_family = "Times")+
    theme(axis.text = element_text(colour = "black"))+
    theme(legend.position = "none")+
    ggtitle("CG: False Negative Rates for Means")+
    theme(plot.title = element_text(size = 18, face = "bold")))

## False Positives
#grid.arrange(covFPFRT,covFPFRTmean,covFPCG,covFPRmean_dub, ncol=2, nrow=2,
  #           top = textGrob("CovGE False Positive Rates",gp=gpar(fontsize=20)))

#grid.arrange(gxeFPFRT, gxeFPFRTmean, gxeFPCG, gxeFPCGmean_dub, ncol=2, nrow=2,
#             top = textGrob("GxE False Positive Rates",gp=gpar(fontsize=20)))

## False Negatives
#grid.arrange(covFNFRT, covFNFRTmean, covFNCG, covFNFRTmean_dub, ncol=2, nrow=2,
   #          top = textGrob("CovGE False Negative Rates",gp=gpar(fontsize=20)))

#grid.arrange(gxeFNFRT, gxeFNFRTmean, gxeFNCG_dub, gxeFNCGmean_dub, ncol=2, nrow=2,
        #     top = textGrob("GxE False Negative Rates",gp=gpar(fontsize=20)))



######################################
## Tradeoff with GxE and Covariance ##
######################################

# Assign 1's for significant outcomes and 0's for nonsignificant outcomes

dat_csv$covtick <- NULL
for(i in 1:nrow(dat_csv)){
  if(dat_csv$covariance_lwrCI[i] > 0 & dat_csv$covariance_uprCI[i] > 0){dat_csv$covtick[i]=1
  }else if(dat_csv$covariance_lwrCI[i] < 0 & dat_csv$covariance_uprCI[i] < 0){dat_csv$covtick[i]=1
  }else{dat_csv$covtick[i]=0} 
}

dat_dub$covtick <- NULL
for(i in 1:nrow(dat_dub)){
  if(dat_dub$covariance_lwrCI[i] > 0 & dat_dub$covariance_uprCI[i] > 0){dat_dub$covtick[i]=1
  }else if(dat_dub$covariance_lwrCI[i] < 0 & dat_dub$covariance_uprCI[i] < 0){dat_dub$covtick[i]=1
  }else{dat_dub$covtick[i]=0} 
}

dat_csv$gxetick <- NULL
for(i in 1:nrow(dat_csv)){
  if(dat_csv$GxE_emm_pvalue[i] > 0.05){dat_csv$gxetick[i]=0}else{dat_csv$gxetick[i]=1}
}

dat_dub$gxetick <- NULL
for(i in 1:nrow(dat_dub)){
  if(dat_dub$GxE_emm_pvalue[i] > 0.05){dat_dub$gxetick[i]=0}else{dat_dub$gxetick[i]=1}
}

sigGxE = dat_csv %>% 
  #filter(sig ==TRUE) %>% # filter out false positives as potential solution to weed out messiness.
  filter(true_GxE_emm != 0) %>%
  filter(true_cov != 0) %>%
  filter(Covconfintboot != "false positive") %>%
  filter(GxEconfintperm != "false positive")


sigGxE2 = dat_dub %>% 
  #filter(sig ==TRUE) %>% # filter out false positives as potential solution to weed out messiness.
  filter(true_GxE_emm != 0) %>%
  filter(true_cov != 0) %>%
  filter(Covconfintboot != "false positive") %>%
  filter(GxEconfintperm != "false positive")

sigGxEdub = fnr.effsize(sigGxE2, metric = "GxE", scenario = 2, data.type = "raw", analysis ="perm",resolution = "fine")
sigGxEdub$ID = rep("GxE_Perm", nrow(sigGxEdub))
sigcovdub = fnr.effsize(sigGxE2, metric = "Cov", scenario = 2,   data.type = "raw", analysis ="perm",resolution = "fine")
sigcovdub$ID = rep("Cov_Perm", nrow(sigcovdub))
trade_dub = rbind(sigcovdub,sigGxEdub)

(bin = ggplot(sigGxE, aes(x = true_GxE_emm, y = covtick))+
    geom_smooth(method = "glm",method.args = list(family = "binomial"),se = T,colour = "black") + 
    geom_point()+
    xlab(expression(""*bar(Delta)*""["GxE"]*": Population"))+ylab(expression("Proportion significant Cov"["GE"]))+
    theme_bw(base_size = 18, base_family = "Times")+
    theme(axis.text.x = element_text(colour = "black"))+
    theme(axis.text.y = element_text(colour = "black"))+
    ggtitle(expression("C   FRT: Significant Cov"["GE"]*" vs. "*bar(Delta)*""["GxE"]))+
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 2)))


(lin = ggplot(filter(sigGxE,replicate ==1), aes(x = true_GxE_emm, y = abs(true_cov)))+
    geom_point(alpha = 0.15)+
    geom_smooth(method = "glm", colour = "black", size = 1.5)+
    xlab(expression(bar(Delta)*""["GxE"]*": Population"))+ylab(expression("|Cov"["GE"]*" |: Population"))+
    theme_bw(base_size = 18, base_family = "Times")+
    theme(axis.text.x = element_text(colour = "black"))+
    theme(axis.text.y = element_text(colour = "black"))+
    ggtitle(expression("A   FRT: |Cov"["GE"]*"| vs. "*bar(Delta)*""["GxE"]))+
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 2)))

(bin2 = ggplot(sigGxE2, aes(x = true_GxE_emm, y = covtick))+
    geom_smooth(method = "glm",method.args = list(family = "binomial"),se = T,colour = "black") + 
    geom_point()+
    xlab(expression(bar(Delta)*""["GxE"]*": Population"))+ylab(expression("Proportion of significant Cov"["GE"]))+
    theme_bw(base_size = 18, base_family = "Times")+
    theme(axis.text.x = element_text(colour = "black"))+
    theme(axis.text.y = element_text(colour = "black"))+
    ggtitle(expression("D   CG: Significant Cov"["GE"]*" vs. "*bar(Delta)*""["GxE"]))+
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 2)))


(lin2 = ggplot(filter(sigGxE2, replicate == 1), aes(x = true_GxE_emm, y = abs(true_cov)))+
    geom_point(alpha = 0.15)+
    geom_smooth(method = "glm", colour = "black", size = 1.5)+
    xlab(expression(bar(Delta)*""["GxE"]*": Population"))+ylab(expression("| Cov"["GE"]*" |: Population"))+
    theme_bw(base_size = 18, base_family = "Times")+
    theme(axis.text.x = element_text(colour = "black"))+
    theme(axis.text.y = element_text(colour = "black"))+
    ggtitle(expression("B   CG: |Cov"["GE"]*"| vs. "*bar(Delta)*""["GxE"]))+
    theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2)))

grid.arrange(lin, lin2, bin, bin2,fndf2,fndfdub, ncol = 2)

### Relationship with Omega^2
dat_csv$sig = NULL
for(i in 1:nrow(dat_csv)){ # Use only if one or the other is significant
  if(dat_csv$GxE_omega_pvalue[i] <= 0.05 | dat_csv$covariance_lwrCI[i] > 0 & dat_csv$covariance_uprCI[i] > 0){dat_csv$sig[i] = TRUE
  }else if(dat_csv$GxE_omega_pvalue[i] <= 0.05 | dat_csv$covariance_lwrCI[i] < 0 & dat_csv$covariance_uprCI[i] < 0){dat_csv$sig[i] = TRUE
  }else{dat_csv$sig[i]=FALSE} 
}

sigGxE = filter(dat_csv, sig ==TRUE)

ggplot(sigGxE, aes(x = GxE_omega, y = covtick))+
  geom_smooth(method = "glm",method.args = list(family = "binomial"),se = T,colour = "black") + 
  geom_point()+
  xlab("Magnitude of GxE - Omega Squared")+ylab("Proportion of significant CovGE values")+
  theme_bw(base_size = 18, base_family = "Times")+
  theme(axis.text.x = element_text(size=16,colour = "black"),
        axis.title.x = element_text(size=18,face="bold")) +
  theme(axis.text.y = element_text(size=16,colour = "black"),
        axis.title.y = element_text(size=18,face="bold")) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2))


###############################################
###      ANOVA vs. GxE_EMM comparison       ###
###############################################
devcol = c("None Significant" = "white", 
           "All Significant" = "#440154FF", 
           "ANOVA P-value < 0.05" = "#39568CFF",
           "Permutation P-value < 0.05" = "#73D055FF")
devshape = c("None Significant" = 21, 
           "All Significant" = 22, 
           "ANOVA P-value < 0.05" = 23,
           "Permutation P-value < 0.05" = 24)
dat_csv$sigsig = NA
dat_csv = dat_csv %>% mutate(sigsig = ifelse(GxE_emm_pvalue <= 0.05 & GxE_omega_pvalue <= 0.05, "All Significant",
                                             ifelse(GxE_emm_pvalue > 0.05 & GxE_omega_pvalue > 0.05, "None Significant",
                                                ifelse(GxE_emm_pvalue > 0.05 & GxE_omega_pvalue <= 0.05, "ANOVA P-value < 0.05", 
                                                   "Permutation P-value < 0.05"))))
oe = dat_csv %>% filter(replicate == 1) #%>% filter(std_dev == 1)
                                        
(eff_var1 = ggplot(oe, aes(x = true_GxE_emm, y = GxE_omega,linetype = factor(std_dev), group = std_dev)) + 
    geom_point(aes(shape = factor(sigsig), fill = factor(std_dev)), alpha = 0.65, size = 4)+
    scale_fill_manual(values = c("0.5" = "#39568CFF", "1" = "#73D055FF"),
                      labels = c("0.5" = "Low (0.5)", "1" = "High (1.0)"))+
    scale_shape_manual(values = devshape)+
    geom_smooth(colour = "yellow1",method = "lm", formula = y ~ splines::bs(x, 3),se = F, size = 1.5) + 
    xlab(expression(""*bar(Delta)*""["GxE"]*" of population"))+
    scale_linetype_manual(values = c("0.5" = "dotdash", "1" = "solid"),
                          labels = c("0.5" = "Low (0.5)", "1" = "High (1.0)"))+
    ylab(expression(""*omega^2*" of population"))+
    labs(fill = "Standard Deviation", shape =  "Significance", linetype = "Standard Deviation")+
    guides(fill = guide_legend(override.aes=list(shape = 21, linetype = c("dotdash", "solid"))))+
    theme_classic(base_size = 18, base_family = "Times")+
    theme(axis.text =element_text(colour = "black")))


fdat_dub$sigsig = NA
dat_dub = dat_dub %>% mutate(sigsig = ifelse(GxE_emm_pvalue <= 0.05 & GxE_omega_pvalue <= 0.05, "Both Significant",
                                             ifelse(GxE_emm_pvalue > 0.05 & GxE_omega_pvalue > 0.05, "None Significant",
                                                    ifelse(GxE_emm_pvalue > 0.05 & GxE_omega_pvalue <= 0.05, "Omega Significant", 
                                                           "E.M.M Significant"))))

(eff_var2 = ggplot(filter(dat_dub,replicate == 1), aes(x = true_GxE_emm, y = GxE_omega, group = factor(std_dev))) + 
    geom_point(aes(shape = factor(sigsig), fill = factor(sigsig)), alpha = 1, size = 4)+
    scale_fill_manual(values = devcol)+
    scale_shape_manual(values = devshape)+
    geom_smooth(aes(linetype = factor(std_dev)),method = "lm", formula = y ~ splines::bs(x, 3),se = F,colour = "black",size = 1.5) + 
    xlab(expression("Actual "*bar(Delta)*""["GxE"]))+
    scale_linetype_manual(values = c("0.5" = "dotted", "1" = "solid"))+
    ylab(expression(""*omega^2))+
    labs(fill = "Significance", shape =  "Significance", linetype = "Standard Deviation")+
    guides(shape = guide_legend(override.aes = list(size = 6)))+
    ggtitle("Paired Common Garden Design")+
    theme_classic(base_size = 24, base_family = "Times")+
    theme(axis.text =element_text(colour = "black")))
grid.arrange(eff_var1,eff_var2)

# Do confidence intervals for each overlap with their true values?  
dat_csv$probgxe = NULL
dat_csv$probcov = NULL

for(i in 1:nrow(dat_csv)){
  if((dat_csv$true_GxE_emm[i] < dat_csv$GxE_emm_lwrCI[i]) | (dat_csv$true_GxE_emm[i] >dat_csv$GxE_emm_uprCI[i])){dat_csv$probgxe[i] = TRUE
  }else{dat_csv$probgxe[i] = FALSE}
  if((dat_csv$true_cov[i] < dat_csv$covariance_lwrCI[i]) | (dat_csv$true_cov[i] >dat_csv$covariance_uprCI[i])){dat_csv$probcov[i] = TRUE
  }else{dat_csv$probcov[i] = FALSE}
}

GxEanoms = dat_csv %>%
  filter(probgxe == TRUE)
dim(GxEanoms)
length(which(GxEanoms$GxE_emm_pvalue<=0.05)) 

cov_anom = dat_csv %>%
  filter(probcov ==TRUE)
dim(cov_anom)
length(which(cov_anom$covariance_pvalue<=0.05))

(overlapper_cov = ggplot(cov_anom[cov_anom$replicate==1,],aes(x = row))+ theme_classic() + 
    geom_errorbar(aes(ymin = covariance_lwrCI,ymax = covariance_uprCI),color = "black")+
    geom_point(aes(y = true_cov),size = 2, color = "red",alpha = 0.5)+
    geom_point(aes(y = covariance,colour = factor(std_dev)), size =2, shape = 4)+
    xlab("Unique Parameter Set (row)") + ylab("Covariance"))

(overlapper_GxE = ggplot(GxEanoms[GxEanoms$replicate==1,],aes(x = row))+ theme_classic() + 
    geom_errorbar(aes(ymin = GxE_emm_lwrCI,ymax = GxE_emm_uprCI),color = "black")+
    geom_point(aes(y = true_GxE_emm),size = 2, color = "red",alpha = 0.5)+
    geom_point(aes(y = GxE_emm,colour = factor(std_dev)), size =2, shape = 4)+
    xlab("Unique Parameter Set (row)") + ylab("GxE - Estimated Marginal Mean"))

#############################
##    Means vs. Raw       ###
#############################

dat_csv$meancoverror = abs(dat_csv$cov_means_uprCI - dat_csv$cov_means_lwrCI)
dat_dub$meancoverror = abs(dat_dub$cov_means_uprCI - dat_dub$cov_means_lwrCI)
dat_csv$coverror = abs(dat_csv$covariance_uprCI - dat_csv$covariance_lwrCI)
dat_dub$coverror = abs(dat_dub$covariance_uprCI - dat_dub$covariance_lwrCI)
dat_csv$meangxeerror = dat_csv$GxE_means_uprCI - dat_csv$GxE_means_lwrCI
dat_dub$meangxeerror = dat_dub$GxE_means_uprCI - dat_dub$GxE_means_lwrCI
dat_csv$gxeerror = dat_csv$GxE_emm_uprCI - dat_csv$GxE_emm_lwrCI
dat_dub$gxeerror = dat_dub$GxE_emm_uprCI - dat_dub$GxE_emm_lwrCI


# Covariance
dat_csv_2 = filter(dat_csv, sample_size != 2)
(covpopcheck = ggplot(dat_csv_2,aes(x = true_cov, y = covariance))+
    geom_point(aes(colour = factor(n_pop),shape = factor(std_dev)),alpha = 0.5)+
    ylab(expression("Cov"["GE"]*": Sample Estimate"))+xlab(expression("Cov"["GE"]*": Population"))+
    geom_abline(slope = 1, intercept = 0,colour = "black",size = 1.5)+
    geom_smooth(method= "glm",colour="red")+
    scale_colour_viridis(discrete = TRUE)+
    theme_classic(base_size = 20, base_family = "Times")+
    ggtitle(expression("A   FRT: Cov"["GE"]))+
    theme(legend.position="none")+
    labs(colour = "Number of genotypes", shape = "Residual error")+
    theme(axis.text = element_text(colour = "black")))

(covpopcheckdub = ggplot(dat_dub,aes(x = true_cov, y = covariance))+
    geom_point(aes(colour = factor(n_pop),shape = factor(std_dev)),alpha = 0.5)+
    ylab(expression("Cov"["GE"]*": Sample Estimate"))+xlab(expression("Cov"["GE"]*": Population"))+
    geom_abline(slope = 1, intercept = 0,colour = "black",size = 1.5)+
    geom_smooth(method= "glm",colour="red")+
    scale_colour_viridis(discrete = TRUE)+
    theme_classic(base_size = 20, base_family = "Times")+
    ggtitle(expression("B   CG: Cov"["GE"]))+
    labs(colour = "Number of genotypes", shape = "Residual error")+
    theme(legend.position="none")+
    theme(axis.text = element_text(colour = "black")))

(covmeancheck = ggplot(dat_csv,aes(x = covariance, y = cov_means))+
    geom_point()+ylab(expression("Cov"["GE"]*": Group Means"))+xlab(expression("Cov"["GE"]*": Raw data"))+
    geom_abline(slope = 1, intercept = 0,colour = "red")+
    theme_classic(base_size = 20, base_family = "Times")+
    ggtitle("Full Reciprocal Transplant Design")+
    theme(axis.text = element_text(colour = "black")))

(coverrorcheck = ggplot(dat_csv, aes(x = coverror,y = meancoverror)) + 
    geom_point(alpha = 0.5)+ylab(expression("Length of Cov"["GE"]*"CI: Group Means"))+xlab(expression("Length of Cov"["GE"]*"CI: Raw data "))+
    geom_abline(slope = 1, intercept = 0,size = 1, colour = "red")+
    theme_classic(base_size = 20, base_family = "Times")+
    ggtitle("")+
    theme(axis.text = element_text(colour = "black")))

grid.arrange(covmeancheck, coverrorcheck)

# Lwr CI
(coverrorcheck_lwr = ggplot(dat_csv, aes(x = covariance_lwrCI, y = cov_means_lwrCI)) + 
    geom_point(alpha = 0.5)+ylab(expression("Lower limit 95% CI for Cov"["GE"]*": Group Means"))+xlab(expression("Lower Limit of 95% CI for Cov"["GE"]*": Raw data "))+
    geom_abline(slope = 1, intercept = 0,size = 1, colour = "red")+
    theme_classic(base_size = 20, base_family = "Times")+
    theme(axis.text = element_text(colour = "black")))

# Upr CI
(coverrorcheck_upr = ggplot(dat_csv, aes(x = covariance_uprCI, y = cov_means_uprCI  )) + 
    geom_point(alpha = 0.5)+ylab(expression("Upper limit 95% CI for Cov"["GE"]*": Group Means"))+xlab(expression("Upper Limit of 95% CI for Cov"["GE"]*": Raw data "))+
    geom_abline(slope = 1, intercept = 0,size = 1, colour = "red")+
    theme_classic(base_size = 20, base_family = "Times")+
    theme(axis.text = element_text(colour = "black")))

grid.arrange(covmeancheck, coverrorcheck_lwr, coverrorcheck_upr)

(covmeancheck2 = ggplot(dat_dub,aes(x = covariance, y = cov_means))+
    geom_point()+ylab(expression("Cov"["GE"]*": Group Means"))+xlab(expression("Cov"["GE"]*": Raw data"))+
    geom_abline(slope = 1, intercept = 0,colour = "red")+
    theme_classic(base_size = 20, base_family = "Times")+
    ggtitle("Paired Common Garden Design")+
    theme(axis.text = element_text(colour = "black")))

(coverrorcheck2 = ggplot(dat_dub, aes(x = coverror,y = meancoverror)) + 
    geom_point(alpha = 0.5)+ylab(expression("Length of Cov"["GE"]*"CI: Group Means"))+xlab(expression("Length of Cov"["GE"]*"CI: Raw data "))+
    geom_abline(slope = 1, intercept = 0,size = 1, colour = "red")+
    theme_classic(base_size = 20, base_family = "Times")+
    ggtitle("")+
    theme(axis.text = element_text(colour = "black")))

grid.arrange(covmeancheck, coverrorcheck,covmeancheck2, coverrorcheck2,ncol = 2)

# GxE
(gxepopcheck = ggplot(dat_csv_2,aes(x = true_GxE_emm, y = GxE_emm))+
    geom_point(aes(colour = factor(n_pop),shape = factor(std_dev)),alpha = 0.5)+
    ylab(expression(""*bar(Delta)*""["GxE"]*": Sample Estimate"))+xlab(expression(""*bar(Delta)*""["GxE"]*": Population"))+
    geom_abline(slope = 1, intercept = 0,colour = "black",size = 1.5)+
    geom_smooth(method= "glm",colour="red")+
    scale_colour_viridis(discrete = TRUE)+
    theme_classic(base_size = 20, base_family = "Times")+
    ggtitle(expression("C   FRT: "*bar(Delta)*""["GxE"]))+
    labs(colour = "Number of genotypes", shape = "Residual error")+
    theme(legend.position="none")+
    theme(axis.text = element_text(colour = "black")))

(gxepopcheck_dub = ggplot(dat_dub,aes(x = true_GxE_emm, y = GxE_emm))+
    geom_point(aes(colour = factor(n_pop),shape = factor(std_dev)),alpha = 0.5)+
    ylab(expression(""*bar(Delta)*""["GxE"]*": Sample Estimate"))+xlab(expression(""*bar(Delta)*""["GxE"]*": Population"))+
    geom_abline(slope = 1, intercept = 0,colour = "black",size = 1.5)+
    geom_smooth(method= "glm",colour="red")+
    scale_colour_viridis(discrete = TRUE)+
    theme_classic(base_size = 20, base_family = "Times")+
    ggtitle(expression("D   CG: "*bar(Delta)*""["GxE"]))+
    labs(colour = "Number of genotypes", shape = "Residual error")+
    theme(axis.text = element_text(colour = "black")))

grid.arrange(covpopcheck, covpopcheckdub, gxepopcheck, gxepopcheck_dub,ncol = 2)

(gxemeancheck = ggplot(dat_csv,aes(x = true_GxE_emm, y = true_GxE_means))+
    geom_point()+ylab(expression(bar(Delta)*""["GxE"]*": Group Means"))+xlab(expression(bar(Delta)*""["GxE"]*": Raw data"))+
    geom_abline(slope = 1, intercept = 0,colour = "red")+
    theme_classic(base_size = 20, base_family = "Times")+
    theme(axis.text = element_text(colour = "black")))

(gxeerrorcheck = ggplot(dat_csv, aes(x = gxeerror,y = meangxeerror)) + 
    geom_point(alpha = 0.5)+ylab(expression("Length of 95% CI for"*bar(Delta)*""["GxE"]*": Group Means"))+xlab(expression("Length of 95% CI for"*bar(Delta)*""["GxE"]*": Raw data"))+
    geom_abline(slope = 1, intercept = 0,size = 1, colour = "red")+
    theme_classic(base_size = 20, base_family = "Times")+
    theme(axis.text = element_text(colour = "black")))

grid.arrange(gxemeancheck, gxeerrorcheck)

(gxeerrorcheck_lwr = ggplot(dat_csv, aes(x = GxE_emm_lwrCI,y = GxE_means_lwrCI)) + 
    geom_point(alpha = 0.5)+ylab(expression("Lower limit of 95% CI for"*bar(Delta)*""["GxE"]*": Group Means"))+xlab(expression("Lower limit of 95% CI for"*bar(Delta)*""["GxE"]*": Raw data"))+
    geom_abline(slope = 1, intercept = 0,size = 1, colour = "red")+
    theme_classic(base_size = 20, base_family = "Times")+
    theme(axis.text = element_text(colour = "black")))

(gxeerrorcheck_upr = ggplot(dat_csv, aes(x = GxE_emm_uprCI,y = GxE_means_uprCI)) + 
    geom_point(alpha = 0.5)+ylab(expression("Upper limit of 95% CI for"*bar(Delta)*""["GxE"]*": Group Means"))+xlab(expression("Upper limit of 95% CI for"*bar(Delta)*""["GxE"]*": Raw data"))+
    geom_abline(slope = 1, intercept = 0,size = 1, colour = "red")+
    theme_classic(base_size = 20, base_family = "Times")+
    theme(axis.text = element_text(colour = "black")))

grid.arrange(covmeancheck, gxemeancheck, coverrorcheck_lwr, gxeerrorcheck_lwr, coverrorcheck_upr,  gxeerrorcheck_upr, ncol = 2)

(gxemeancheck2 = ggplot(dat_dub,aes(x = true_GxE_emm, y = true_GxE_means))+
    geom_point()+theme_classic()+ylab(expression(bar(Delta)*""["GxE"]*": Group Means"))+xlab(expression(bar(Delta)*""["GxE"]*": Raw data"))+
    geom_abline(slope = 1, intercept = 0,colour = "red")+
    theme_classic(base_size = 20, base_family = "Times")+
    theme(axis.text = element_text(colour = "black")))

(gxeerrorcheck2 = ggplot(dat_dub, aes(x = gxeerror,y = meangxeerror)) + 
    geom_point(alpha = 0.5)+ylab(expression("Length of 95% CI for"*bar(Delta)*""["GxE"]*": Group Means"))+xlab(expression("Length of 95% CI for"*bar(Delta)*""["GxE"]*": Raw data"))+
    geom_abline(slope = 1, intercept = 0,size = 1, colour = "red")+
    theme_classic(base_size = 20, base_family = "Times")+
    theme(axis.text = element_text(colour = "black")))

grid.arrange(gxemeancheck2, gxeerrorcheck2)

# Do GxE pvalues from permutation match pvalues from anova? 

(pval1 = ggplot(filter(dat_csv, replicate ==1),aes(x = GxE_Anova, y = GxE_emm_pvalue,group = factor(total_samples),colour = factor(sample_size),shape = factor(n_pop)))+ # means data
    geom_point(size = 2, alpha = 0.75)+theme_classic(base_size = 20, base_family = "Times")+ylab("GxE EMM P-value")+xlab("GxE Anova P-value")+
    geom_vline(xintercept = 0.05,colour = "red")+labs(shape = "Number of Populations",colour = "Sample Size")+
    geom_hline(yintercept = 0.05,colour = "red"))

(pval2 = ggplot(filter(dat_csv, replicate ==1),aes(x = GxE_Anova, y = GxE_means_pvalue,group = factor(total_samples),colour = factor(sample_size), shape = factor(n_pop)))+ # means data
    geom_point(size = 2, alpha = 0.75)+theme_classic(base_size = 20, base_family = "Times")+ylab("GxE Means P-value")+xlab("GxE Anova P-value")+
    geom_vline(xintercept = 0.05,colour = "red")+labs(shape = "Number of Populations",colour = "Sample Size")+
    geom_hline(yintercept = 0.05,colour = "red"))

(pval3 = ggplot(filter(dat_dub, replicate ==1),aes(x = GxE_Anova, y = GxE_emm_pvalue,colour = GxE_emm,shape = factor(n_pop)))+ # means data
    geom_point()+theme_classic(base_size = 20, base_family = "Times")+ylab("GxE EMM P-value")+xlab("GxE Anova P-value")+
    geom_vline(xintercept = 0.05,colour = "red")+labs(shape = "Number of Populations",colour = "Delta GxE")+
    geom_hline(yintercept = 0.05,colour = "red"))

(pval4 = ggplot(filter(dat_dub, replicate ==1),aes(x = GxE_Anova, y = GxE_means_pvalue,colour = GxE_means, shape = factor(n_pop)))+ # means data
    geom_point()+theme_classic(base_size = 20, base_family = "Times")+ylab("GxE Means P-value")+xlab("GxE Anova P-value")+
    geom_vline(xintercept = 0.05,colour = "red")+labs(shape = "Number of Populations",colour = "Delta GxE -- Group Means")+
    geom_hline(yintercept = 0.05,colour = "red"))

grid.arrange(pval1,pval2)
grid.arrange(pval1,pval2,pval3,pval4,ncol = 2)

# For those that don't match, is there any pattern ?

suspect.pvals = dat_csv %>% # raw 
  filter(GxE_emm_pvalue >= 0.05) %>%
  filter(GxE_Anova <= 0.05)

(ggplot(suspect.pvals,aes(x = GxE_Anova, y = GxE_emm_pvalue,colour = GxE_emm))+
    geom_point()+theme_classic()+ylab("GxE EMM Pvalue")+xlab("GxE Anova Pvalue")+
    geom_hline(yintercept = 0.05,colour = "red"))

suspect.pvals.mean = dat_csv %>% # mean
  filter(GxE_means_pvalue >= 0.05) %>%
  filter(GxE_Anova <= 0.05)

(ggplot(suspect.pvals.mean, aes(x = GxE_means, y = GxE_emm_pvalue)) + 
    geom_point()+theme_classic())

#############################
##      Phenotype plots    ##
#############################

# To find examples  
(rowpicker = dat_dub %>%
    filter(sample_size == 4)%>%
    filter(n_pop == 4)%>%
    #filter(GxE_emm < 0.25)%>% #No GxE
    filter(GxE_emm > 0.70)%>% # GxE
    filter(covariance < -.75)) # CnGV
    #filter(covariance > .75)) # CoGV

chosen_1 = c(16434, # CoGV, GxE, 2 pop #
             11021, # CnGV, GxE, 2 pop #
             8243, # CnGV, No GxE, 2 pop #
             1018, # CoGV, No GxE, 2 pop #
             
            # 12949, # CoGV, No GxE, 8 pop #
            # 16546, # CnGV, No GxE, 8 pop #
            # 13709, # CnGV, GxE, 8 pop #
            # 17324, # CoGV, GxE, 8 pops #
             
            12867, # CoGV, No GxE, 4 pop #
            299, # CnGV, No GxE, 4 pop #
            6484, # CnGV, GxE, 4 pop #
            13662) # CoGV, GxE, 4 pops #

chosen_2 = c(25869, # CoGV, GxE, 2 pop per env #
             30830, # CnGV, GxE, 2 pop per env #
             24478, # CnGV, No GxE, 2 pop per env #
             25827) # CoGV, No GxE, 2 pop per env #
             
            # 38536, # CoGV, No GxE, 8 pop per env #
           #  42175, # CnGV, No GxE, 8 pop #
            # 19648, # CnGV, GxE, 8 pop #
           #  32187) # CoGV, GxE, 8 pops #

chosen = c(chosen_1, chosen_2)

# Plotting Specs
short_env = c("E_1" = "Env. 1", "E_2" = "Env. 2")
short_gen = c("E_1" = "Native to Env. 1", "E_2" = "Native to Env. 2")

long_env = c("E_1" = "Env. 1", "E_2" = "Env. 2","E_3" = "Env. 3", "E_4" = "Env. 4","E_5" = "Env. 5", "E_6" = "Env. 6","E_7" = "Env. 7", "E_8" = "Env. 8")
mid_env = c("E_1" = "Env. 1", "E_2" = "Env. 2","E_3" = "Env. 3", "E_4" = "Env. 4")
mid_gen = c("E_1" = "Native to Env. 1", "E_2" = "Native to Env. 2","E_3" = "Native to Env. 3", "E_4" = "Native to Env. 4")

long_gen = c("E_1" = "Gen. 1\n(Native to Env. 1)", "E_2" = "Gen. 2\n(Native to Env. 2)","E_3" = "Gen. 3", "E_4" = "Gen. 4","E_5" = "Gen. 5", "E_6" = "Gen. 6","E_7" = "Gen. 7", "E_8" = "Gen. 8")
long_genlab = c("E_1" = "Native to Env. 1", "E_2" = "Native to Env. 2","E_3" = "Native to Env. 3", "E_4" = "Native to Env. 4",
                "E_5" = "Native to Env. 5", "E_6" = "Native to Env. 6","E_7" = "Native to Env. 7", "E_8" = "Native to Env. 8")


#for(i in 1:length(chosen)){
#phenRow = filter(start_df,row == chosen[i])
#plotdat = filter(phen_data, row == chosen[i])
(chosen = filter(dat_dub, row == 31969))
(phenRow = chosen)
plotdat = filter(phen_data, row == chosen$row)
label1 = paste0(phenRow$covariance,"; P = ",phenRow$covariance_pvalue)
label2 = paste0(phenRow$GxE_emm,"; P = ",phenRow$GxE_emm_pvalue)
colorpal = c("E_1" = "#3CBB75FF", "E_2" = "#453781FF")
shape_2 = c("E_1" = 15, "E_2" =  17)
shape_4 = c("E_1" = 15, "E_2" = 17,"E_3" = 19, "E_4"= 18)
shape_8 = c("E_1" = 15, "E_2" = 17,"E_3" = 19, "E_4"= 18,
            "E_5" = 15, "E_6" = 17,"E_7" = 19, "E_8"= 18)

p = ggplot(plotdat,aes(x = exp_env_factor, y = phen_corrected, group = gen_factor,colour = nat_env_factor))+
  geom_jitter(aes(shape = nat_env_factor),size = 5,width = 0.07)+ #position=position_dodge(width = 0.15))+
  geom_smooth(size = 4, se=FALSE)+
  theme_classic(base_size = 40, base_family = "Times")+
  ylim(-4,6)+
  labs(colour = "",shape = "")+
  guides(colour = guide_legend(ncol = 2))+
  ylab("Phenotype")+xlab("Environment")+
  #annotate("text", x = 0.5, y = 5.5, label = deparse(bquote("Cov"["GE"] ==~.(label1))),size=16, color = "black",hjust = 0,parse = T)+
  #annotate("text", x = 0.5, y = 4, label = deparse(bquote(bar(Delta)*""["GxE"] ==~.(label2))),size=16 , color = "black", hjust = 0,parse = T)+
  theme(axis.text=element_text(colour="black"))+
  theme(legend.position="bottom")

p1 = p + if(phenRow$n_pop == 2){scale_x_discrete(labels = short_env)
  }else if(phenRow$env_scenario == 2){scale_x_discrete(labels = short_env)
  }else{scale_x_discrete(labels = long_env)} 

p2 = p1 + if(phenRow$n_pop == 2){
  scale_colour_manual(labels = short_gen, values = colorpal)
}else if(phenRow$env_scenario == 2){scale_colour_manual(values = colorpal,labels = short_gen)
}else{scale_colour_viridis(discrete = TRUE,labels = long_genlab)}

p3 = p2 + if(phenRow$n_pop == 2){
  scale_shape_manual(values = shape_2,labels = short_gen)
}else if(phenRow$env_scenario == 2){scale_shape_manual(values = shape_2,labels = short_gen)
}else if(phenRow$n_pop == 4){scale_shape_manual(values = shape_4,labels = mid_gen)
}else{scale_shape_manual(values = shape_8,labels = long_genlab)}
p3
pdf(paste("Row_", phenRow$row, ".pdf", sep = ""), width=11, height=8.5) # start export
print(p3) 
dev.off() # finish export
#}
phenRow

# Check with Variance partition data frame
vp = read.csv("~/Desktop/Variance_output_results.csv")

filter(vp, row == chosen[2])
filter(dat_csv, row == chosen[2])

#model_df$phen_corrected = rep(c(-1,1,-.95,1.05),each = 5) For Vg, Ve, Vgxe
#Otherwise use simulated data from Cov_GxE_clusterFun.R

ggplot(model_df, aes(x = exp_env_factor, y = phen_corrected, group = gen_factor, shape = gen_factor, fill = gen_factor, colour = gen_factor)) + 
  geom_point(size = 5) + geom_line(size=1,lty = "dashed") + 
  scale_shape_manual(values = shape1,
                     labels = c("G_1" = "Genotype 1", "G_2" = "Genotype 2"))+
  scale_fill_manual(values = ColorFill,
                    labels = c("G_1" = "Genotype 1", "G_2" = "Genotype 2"))+
  scale_colour_manual(values = ColorFill,
                      labels = c("G_1" = "Genotype 1", "G_2" = "Genotype 2"))+
  ylab("Phenotype") + xlab(" ")+
  
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.text.y = element_blank(),
        #axis.text.x = element_blank(),
        axis.ticks.y = element_blank())+
  scale_x_discrete(labels=c("E_1" = "Environment 1","E_2" = "Environment 2"))+
  #theme(axis.title.x = element_text(face = "bold", size=26))+
  #theme(axis.title.y = element_text(face = "bold", size=26))+
  theme(axis.text= element_text(colour = "black"))+
  labs(colour = "Genotype",fill = "Genotype",shape = "Genotype")+
  theme(legend.position="none")





################################
##      Applied Data Plots    ##
################################
# Load Functions for meta-analysis 
setwd("~/Documents/GitHub/CnGV/CnGV/src/")
source("CovarianceDataFunctions.R")

# Set bootstraps
n_boot = 999

# Molly's age at metamorphosis data
mm = read.csv("~/Desktop/Work/DataSets/Tadpole Plasticity_2017/mortality2017.csv")
mm = mm[-which(is.na(mm$Jul_metamorph)),] # Include only those that metamorphosed
mm$age = mm$Jul_metamorph-mm$Jul_hatch # Calculate age at MM

# Wrangle data for analysis
mm1 = mm %>%
  filter(Pop != "BELL") %>% # Exclude bellamy bc missing data at 6ppt
  filter(tad %in% c(0,6)) %>% # Only use 0, 6ppt
  droplevels()
mm1$gen_factor = paste0("G_",as.numeric(factor(mm1$Pop)))
mm1$exp_env_factor = paste0("E_",as.numeric(as.factor(mm1$tad))) 
mm1$nat_env_factor = NULL
for(i in 1:nrow(mm1)){
  if(mm1$Pop[i] == "BOD" | mm1$Pop[i] == "CSI" |mm1$Pop[i] == "LH" |mm1$Pop[i] == "DQ"){ mm1$nat_env_factor[i] = "E_2"
  }else{ mm1$nat_env_factor[i] = "E_1"}
}

# Format final dataset
ma = data.frame("data_type" = rep("raw", nrow(mm1)),
                "gen_factor" = factor(mm1$gen_factor),
                "exp_env_factor" = factor(mm1$exp_env_factor), 
                "nat_env_factor" = factor(mm1$nat_env_factor), # E_2 = coastal; E_1 = inland
                "phen_data"= mm1$age)
dim(ma)
totalss = ma %>%
  group_by(gen_factor,nat_env_factor)%>%
  summarize("ss" = n())
length(unique(ma$gen_factor))*length(unique(ma$exp_env_factor))*mean(totalss$ss)

# Standardize by standard deviation of group means
ma$group = paste(ma$gen_factor,ma$exp_env_factor,sep = "-")
ma$phen_corrected = (ma$phen_data - mean(ma$phen_data))/sd(tapply(ma$phen_data, ma$group, mean))
MollyTest = amarillo_armadillo(ma, n_boot, data_type)

# Plot Molly's study
label1 = MollyTest$Covariance.Estimate
label2 = paste0(":  C.I. = ", MollyTest$Covariance.Lower.CI," - ", MollyTest$Covariance.Upper.CI)
label3 = round(MollyTest$GxE.Estimate,3)
label4 = paste0(":  P-value = ",round(MollyTest$GxE.p.value,3))

(mollyplot = ggplot(ma, aes(x = exp_env_factor, y = phen_corrected, group = gen_factor)) + 
    geom_smooth(aes(colour = nat_env_factor),size = 2, method = "glm",se=FALSE)+
    geom_point(aes(shape = gen_factor, fill = gen_factor),size = 4, alpha = 0.7,position = position_dodge(width = 0.25))+   
    xlab("")+ylab("Normalized Phenotype \n(age at metamorphosis in days)")+
    scale_x_discrete(breaks=c("E_1","E_2"),labels=c("Freshwater", "Saltwater"))+
    scale_colour_manual(values =  c("E_1" = "#3CBB75FF", "E_2" = "#453781FF"),
                        labels = c("E_1" = "Inland populations", "E_2" = "Coastal populations"))+
    scale_fill_manual(values =  c("G_1" = "#453781FF","G_2" = "#453781FF","G_3" = "#453781FF","G_4" = "#453781FF","G_5" = "#3CBB75FF", "G_6" ="#3CBB75FF","G_7" =  "#3CBB75FF"),
                      labels = c("G_1" = "Genotype1","G_2" = "Genotype2","G_3" = "Genotype3","G_4" = "Genotype4",
                                 "G_5" = "Genotype5", "G_6" = "Genotype6","G_7" = "Genotype7"))+
    scale_shape_manual(values = c("G_1" = 21,"G_2" = 22,"G_3" = 24,"G_4" = 23,"G_5" = 21, "G_6" = 22,"G_7" = 24),
                       labels = c("G_1" = "Genotype1","G_2" = "Genotype2","G_3" = "Genotype3","G_4" = "Genotype4",
                                  "G_5" = "Genotype5", "G_6" = "Genotype6","G_7" = "Genotype7"))+
    labs(col=" ",fill=" ",shape = " ")+
    theme_classic(base_size = 20, base_family = "Times")+
    theme(legend.position="none")+
    ggtitle("A")+
    theme(axis.text.x=element_text(colour="black"))+
    theme(axis.text.y=element_text(colour="black")))
    #annotate("text", x = 1, y = 3.5, label = deparse(bquote("Cov"["GE"] ==~.(label1))),size=6, color = "black",hjust = 0,parse = T)+
   # annotate("text", x = 1.35, y = 3.5, label = label2 ,size=6, color = "black",hjust = 0,parse = F)+
    #annotate("text", x = 1, y = 3, label = deparse(bquote(bar(Delta)*""["GxE"] ==~.(label3))),size=6 , color = "black", hjust = 0,parse = T)+
    #annotate("text", x = 1.3, y = 3, label = label4 ,size=6, color = "black",hjust = 0,parse = F))
#Save as PDF A4

## Wing Length Study
setwd("~/Documents/GitHub/CnGV/CnGV/data/")
data.df1 <- read.csv("meta_df.csv")
data.df1$Data_file_name <- as.character(data.df1$Data_file_name)
malewings = filter(data.df1, data.df1$Data_file_name == "630_Cenzer_male_wing_length")
femalewings = filter(data.df1, data.df1$Data_file_name == "630_Cenzer_female_wing_length")

# Format data for each
malewings$gen_factor = factor(malewings$gen_factor)
malewings$exp_env_factor = factor(malewings$exp_env_factor)
malewings$nat_env_factor = factor(malewings$nat_env_factor)

femalewings$gen_factor = factor(femalewings$gen_factor)
femalewings$exp_env_factor = factor(femalewings$exp_env_factor)
femalewings$nat_env_factor = factor(femalewings$nat_env_factor)

# Standardize by standard deviation of group means
malewings$group = paste(malewings$gen_factor,malewings$exp_env_factor,sep = "-")
malewings$phen_corrected = (malewings$phen_data - mean(malewings$phen_data))/sd(tapply(malewings$phen_data, malewings$group, mean))
femalewings$group = paste(femalewings$gen_factor,femalewings$exp_env_factor,sep = "-")
femalewings$phen_corrected = (femalewings$phen_data - mean(femalewings$phen_data))/sd(tapply(femalewings$phen_data, femalewings$group, mean))

# Run tests
MaleWingsTest = amarillo_armadillo(malewings, n_boot, data_type)
FemaleWingsTest = amarillo_armadillo(femalewings, n_boot, data_type)

totalss = femalewings %>%
  group_by(gen_factor,nat_env_factor)%>%
  summarize("ss" = n())
length(unique(femalewings$gen_factor))*length(unique(femalewings$exp_env_factor))*mean(totalss$ss)



# Plot MaleWings
(malewingplot = ggplot(malewings, aes(x = exp_env_factor, y = phen_corrected, group = gen_factor)) + 
    geom_smooth(aes(colour = nat_env_factor),size = 2, method = "glm",se=FALSE)+
    geom_point(aes(shape = gen_factor, fill = gen_factor),size = 4, alpha = 0.7,position = position_dodge(width = 0.25))+
    xlab("")+ylab("Normalized Phenotype \n(forewing length in mm)")+
    scale_shape_manual(values = c("G_1" = 21,"G_2" = 22,"G_3" = 24,"G_4" = 21,"G_5" = 22, "G_6" = 23,"G_7" = 24,"G_8" = 25),
                       labels = c("G_1" = "Genotype1","G_2" = "Genotype2","G_3" = "Genotype3","G_4" = "Genotype4",
                                  "G_5" = "Genotype5", "G_6" = "Genotype6","G_7" = "Genotype7","G_8" = "Genotype8"))+
    scale_x_discrete(breaks=c("E_1","E_2"),labels=c("Native \n host plant", "Introduced \n host plant"))+
    scale_colour_manual(values =  c("E_1" = "#3CBB75FF", "E_2" = "#453781FF"),labels = male_labels)+
    scale_fill_manual(values = c("G_1" = "#3CBB75FF","G_2" = "#3CBB75FF","G_3" = "#3CBB75FF","G_4" = "#453781FF",
                                 "G_5" = "#453781FF", "G_6" = "#453781FF","G_7" = "#453781FF","G_8" = "#453781FF"),
                      labels = c("G_1" = "Genotype1","G_2" = "Genotype2","G_3" = "Genotype3","G_4" = "Genotype4",
                                 "G_5" = "Genotype5", "G_6" = "Genotype6","G_7" = "Genotype7","G_8" = "Genotype8"))+
    labs(col=" ",fill=" ",shape = " ")+
    theme_classic(base_size = 20, base_family = "Times")+
    theme(axis.text.x=element_text(colour="black"))+
    ggtitle("B  Male")+
    theme(legend.position="none")+    
    theme(axis.text.y=element_text(colour="black")))

# Plot femaleWings
(femalewingplot = ggplot(femalewings, aes(x = exp_env_factor, y = phen_corrected, group = gen_factor)) + 
    geom_smooth(aes(colour = nat_env_factor),size = 2, method = "glm",se=FALSE)+
    geom_point(aes(shape = gen_factor, fill = gen_factor),size = 4, alpha = 0.7,position = position_dodge(width = 0.25))+
    xlab("")+ylab("Normalized Phenotype \n(forewing length in mm)")+
    scale_shape_manual(values = c("G_1" = 21,"G_2" = 22,"G_3" = 24,"G_4" = 21,"G_5" = 22, "G_6" = 23,"G_7" = 24,"G_8" = 25),
                       labels = c("G_1" = "Genotype1","G_2" = "Genotype2","G_3" = "Genotype3","G_4" = "Genotype4",
                                  "G_5" = "Genotype5", "G_6" = "Genotype6","G_7" = "Genotype7","G_8" = "Genotype8"))+
    scale_x_discrete(breaks=c("E_1","E_2"),labels=c("Native \n host plant", "Introduced \n host plant"))+
    scale_colour_manual(values =  c("E_1" = "#3CBB75FF", "E_2" = "#453781FF"),labels = male_labels)+
    scale_fill_manual(values = c("G_1" = "#3CBB75FF","G_2" = "#3CBB75FF","G_3" = "#3CBB75FF","G_4" = "#453781FF",
                                 "G_5" = "#453781FF", "G_6" = "#453781FF","G_7" = "#453781FF","G_8" = "#453781FF"),
                      labels = c("G_1" = "Genotype1","G_2" = "Genotype2","G_3" = "Genotype3","G_4" = "Genotype4",
                                 "G_5" = "Genotype5", "G_6" = "Genotype6","G_7" = "Genotype7","G_8" = "Genotype8"))+
    labs(col=" ",fill=" ",shape = " ")+
    theme_classic(base_size = 20, base_family = "Times")+
    theme(axis.text.x=element_text(colour="black"))+
    ggtitle("C  Female")+
    theme(legend.position="none")+    
    theme(axis.text.y=element_text(colour="black")))
    
grid.arrange(mollyplot,malewingplot,femalewingplot,
             layout_matrix = cbind(c(1,2),c(1,3))) #Device size 11x11
#########################################
##      Supplemental Materials Plots   ##
#########################################

## False Negative Heatmap - GxE 
allgxe = rbind(gxeperm1,gxeboot1,gxeanova1)

allgxe[is.nan(allgxe)] <- 0
gxeFNR1 = allgxe %>%
  filter(between(bin,0.3,0.6))%>%  
  
  group_by(sample_size,n_pop,ID) %>%
  summarize("fnr" = mean(fnr))
gxeFNR1$n_pop<-as.factor(gxeFNR1$n_pop)
gxeFNR2 <- 
  gxeFNR1 %>% 
  mutate(#xadj = ifelse(ID == "GxE_Anova", -.45/3,ifelse(ID == "GxE_Perm", 0.0, .45/3)),
    yadj = ifelse(ID == "GxE_Anova", -0.325,
                  ifelse(ID == "GxE_Perm", 0, 0.325)),
    # ifelse(ID == "GxE_Perm", 0.0, 0.9/3)),
    xpos = as.numeric(factor(sample_size)),
    ypos = as.numeric(n_pop)+yadj,
    col = ifelse(fnr >= 0.5, "black","white"),
    ID2 = ifelse(ID == "GxE_Anova", "Anova",  
                 ifelse(ID == "GxE_Perm", "Perm.", "Boot.")))

(gxeFNFRT = gxeFNR2 %>% 
    ggplot(aes(as.factor(sample_size), ypos, fill = fnr)) + 
    geom_tile(height = 0.33, width = 0.95, color= "white") +
    geom_text(aes(label = paste0(ID2,"\n",round(fnr, 2)), colour = col), 
              size = 4, family = "Times", show.legend = F) +
    scale_colour_identity()+
    scale_fill_viridis(breaks=seq(0,1,0.25), #breaks in the scale bar
                       limits=c(0,1))+
    scale_x_discrete(name = "Sample Size",
                     labels = unique(gxeFNR2$sample_size)) +
    scale_y_continuous(breaks = 1:3, name = "Number of Genotypes",
                       labels = c(2,4,8)) +
    # ggtitle(expression("Reciprocal Transplant: "*bar(Delta)*""["GxE"]*" False Positive Rates"))+  
    labs(fill = "False Negative Rate")+
    theme_classic(base_size = 18, base_family = "Times")+
    theme(axis.text = element_text(colour = "black"))+
    theme(legend.position = "none")+
    ggtitle("FRT: False Negative Rates")+
    theme(plot.title = element_text(size = 18, face = "bold")))


## CG False Negative Heatmap - GxE 
allgxe_dub = rbind(gxeperm2,gxeboot2,gxeanova2)


allgxe_dub[is.nan(allgxe)] <- 0
gxeFNRdub = allgxe_dub %>%
  filter(between(bin,0.3,0.6))%>%  
  
  group_by(sample_size,n_pop,ID) %>%
  summarize("fnr" = mean(fnr))
gxeFNRdub$n_pop<-as.factor(gxeFNRdub$n_pop)
gxeFNRdub2 <- 
  gxeFNRdub %>% 
  mutate(#xadj = ifelse(ID == "GxE_Anova", -.45/3,ifelse(ID == "GxE_Perm", 0.0, .45/3)),
    yadj = ifelse(ID == "GxE_Anova", -0.325,
                  ifelse(ID == "GxE_Perm", 0, 0.325)),
    # ifelse(ID == "GxE_Perm", 0.0, 0.9/3)),
    xpos = as.numeric(factor(sample_size)),
    ypos = as.numeric(n_pop)+yadj,
    col = ifelse(fnr >= 0.5, "black","white"),
    ID2 = ifelse(ID == "GxE_Anova", "Anova",  
                 ifelse(ID == "GxE_Perm", "Perm.", "Boot.")))

(gxeFNRdub2plot = gxeFNRdub2 %>% 
    ggplot(aes(as.factor(sample_size), ypos, fill = fnr)) + 
    geom_tile(height = 0.33, width = 0.95, color= "white") +
    geom_text(aes(label = paste0(ID2,"\n",round(fnr, 2)), colour = col), 
              size = 4, family = "Times", show.legend = F) +
    scale_colour_identity()+
    scale_fill_viridis(breaks=seq(0,1,0.25), #breaks in the scale bar
                       limits=c(0,1))+
    scale_x_discrete(name = "Sample Size",
                     labels = unique(gxeFNRdub2$sample_size)) +
    scale_y_continuous(breaks = 1:3, name = "Number of Genotypes",
                       labels = c(4,8,16)) +
    # ggtitle(expression("Reciprocal Transplant: "*bar(Delta)*""["GxE"]*" False Positive Rates"))+  
    labs(fill = "False Negative Rate")+
    theme_classic(base_size = 18, base_family = "Times")+
    theme(axis.text = element_text(colour = "black"))+
    theme(legend.position = "none")+
    ggtitle("CG: False Negative Rates")+
    theme(plot.title = element_text(size = 18, face = "bold")))

## FRT - False Positive Heatmap - GxE 
gxe_fp = rbind(raw_conf3,raw_conf4,raw_conf5)

gxe_fp[is.nan(gxe_fp)] <- 0
gxe_fp$n_pop<-as.factor(gxe_fp$n_pop)

gxe_fp1 <- 
  gxe_fp %>% 
  filter(name == "False Positive")%>%
  mutate(#xadj = ifelse(ID == "GxE_Anova", -.45/3,ifelse(ID == "GxE_Perm", 0.0, .45/3)),
    yadj = ifelse(ID == "GxE_Anova", -0.325,
                  ifelse(ID == "GxE_Perm", 0, 0.325)),
    # ifelse(ID == "GxE_Perm", 0.0, 0.9/3)),
    xpos = as.numeric(factor(sample_size)),
    ypos = as.numeric(n_pop)+yadj,
    col = ifelse(rate >= 0.5, "black","white"),
    ID2 = ifelse(ID == "GxE_Anova", "Anova",  
                 ifelse(ID == "GxE_Perm", "Perm.", "Boot.")))

(gxeFPFRT = gxe_fp1 %>% 
    ggplot(aes(as.factor(sample_size), ypos, fill = rate)) + 
    geom_tile(height = 0.33, width = 0.95, color= "white") +
    geom_text(aes(label = paste0(ID2,"\n",round(rate, 2)), colour = col), 
              size = 4, family = "Times", show.legend = F) +
    scale_colour_identity()+
    scale_fill_viridis(breaks=seq(0,1,0.25), #breaks in the scale bar
                       limits=c(0,1))+
    scale_x_discrete(name = "Sample Size",
                     labels = unique(gxe_fp1$sample_size)) +
    scale_y_continuous(breaks = 1:3, name = "Number of Genotypes",
                       labels = c(2,4,8)) +
    # ggtitle(expression("Reciprocal Transplant: "*bar(Delta)*""["GxE"]*" False Positive Rates"))+  
    labs(fill = "False Positive Rate")+
    theme_classic(base_size = 18, base_family = "Times")+
    theme(axis.text = element_text(colour = "black"))+
    theme(legend.position = "none")+
    ggtitle("FRT: False Positive Rates")+
    theme(plot.title = element_text(size = 18, face = "bold")))


## CG - False Positive Heatmap - GxE 
dub_gxefp = rbind(dub_raw_conf3,dub_raw_conf4,dub_raw_conf5)
dub_gxefp = dub_gxefp %>% filter(name == "False Positive") 
dub_gxefp[is.nan(dub_gxefp)] <- 0
dub_gxefp1 <- 
  dub_gxefp %>% 
  mutate(#xadj = ifelse(ID == "GxE_Anova", -.45/3,ifelse(ID == "GxE_Perm", 0.0, .45/3)),
    yadj = ifelse(ID == "GxE_Anova", -0.325,
                  ifelse(ID == "GxE_Perm", 0, 0.325)),
    #ifelse(ID == "GxE_Perm", 0.0, 0.9/3)),
    xpos = as.numeric(factor(sample_size)),
    ypos = as.numeric(factor(n_pop))+yadj,
    col = ifelse(rate >= 0.5, "black","white"),
    ID2 = ifelse(ID == "GxE_Anova", "Anova",  
                 ifelse(ID == "GxE_Perm", "Perm.", "Boot.")))

(gxeFPCG = dub_gxefp1 %>% 
    ggplot(aes(as.factor(sample_size), ypos, fill = rate)) + 
    geom_tile(height = 0.33, width = 0.95, color= "white") +
    geom_text(aes(label = paste0(ID2,"\n",round(rate, 2)), colour = col), 
              size = 4, family = "Times", show.legend = F) +
    scale_colour_identity()+
    scale_fill_viridis(breaks=seq(0,1,0.25), #breaks in the scale bar
                       limits=c(0,1))+
    scale_x_discrete(name = "Sample Size",
                     labels = unique(gxeFPR1$sample_size)) +
    scale_y_continuous(breaks = 1:3, name = "Number of Genotypes",
                       labels = c(4,8,16)) +
    #ggtitle(expression("Common Garden: "*bar(Delta)*""["GxE"]*" False Positive Rates"))+
    labs(fill = "False Positive Rate")+
    theme_classic(base_size = 18, base_family = "Times")+
    theme(axis.text = element_text(colour = "black"))+
    theme(legend.position = "none")+
    ggtitle("CG: False Positive Rates")+
    theme(plot.title = element_text(size = 18, face = "bold"))) 

grid.arrange(gxeFNFRT,gxeFNRdub2plot,gxeFPFRT, gxeFPCG)

#######################################
#########     Extra Code      #########
#######################################


## Covariance x GxE ##


# Assign colors for Cov (CI) x GxE (pval) plot
dat_csv$col = NULL
for(i in 1:nrow(dat_csv)){
  if(dat_csv$GxE_Anova[i] <= 0.05 & dat_csv$covariance_lwrCI[i] > 0 & dat_csv$covariance_uprCI[i] > 0){dat_csv$col[i] = "red" # Both significant
  }else if(dat_csv$GxE_Anova[i] <= 0.05 & dat_csv$covariance_lwrCI[i] < 0 & dat_csv$covariance_uprCI[i] < 0){dat_csv$col[i] = "red"
  }else if(dat_csv$GxE_Anova[i] > 0.05 & dat_csv$covariance_lwrCI[i] > 0 & dat_csv$covariance_uprCI[i] > 0){dat_csv$col[i] = "darkgreen" # Cov significant
  }else if(dat_csv$GxE_Anova[i] > 0.05 & dat_csv$covariance_lwrCI[i] < 0 & dat_csv$covariance_uprCI[i] < 0){dat_csv$col[i] = "darkgreen" 
  }else if(dat_csv$GxE_Anova[i] <= 0.05 & dat_csv$covariance_lwrCI[i] < 0 & dat_csv$covariance_uprCI[i] > 0){dat_csv$col[i] = "dodgerblue4" # GxE significant
  }else if(dat_csv$GxE_Anova[i] <= 0.05 & dat_csv$covariance_lwrCI[i] > 0 & dat_csv$covariance_uprCI[i] < 0){dat_csv$col[i] = "dodgerblue4"
  }else{dat_csv$col[i] = "grey"} # None significant
}

dat_dub$col = NULL
for(i in 1:nrow(dat_dub)){
  if(dat_dub$GxE_Anova[i] <= 0.05 & dat_dub$covariance_lwrCI[i] > 0 & dat_dub$covariance_uprCI[i] > 0){dat_dub$col[i] = "red" # Both significant
  }else if(dat_dub$GxE_Anova[i] <= 0.05 & dat_dub$covariance_lwrCI[i] < 0 & dat_dub$covariance_uprCI[i] < 0){dat_dub$col[i] = "red"
  }else if(dat_dub$GxE_Anova[i] > 0.05 & dat_dub$covariance_lwrCI[i] > 0 & dat_dub$covariance_uprCI[i] > 0){dat_dub$col[i] = "darkgreen" # Cov significant
  }else if(dat_dub$GxE_Anova[i] > 0.05 & dat_dub$covariance_lwrCI[i] < 0 & dat_dub$covariance_uprCI[i] < 0){dat_dub$col[i] = "darkgreen" 
  }else if(dat_dub$GxE_Anova[i] <= 0.05 & dat_dub$covariance_lwrCI[i] < 0 & dat_dub$covariance_uprCI[i] > 0){dat_dub$col[i] = "dodgerblue4" # GxE significant
  }else if(dat_dub$GxE_Anova[i] <= 0.05 & dat_dub$covariance_lwrCI[i] > 0 & dat_dub$covariance_uprCI[i] < 0){dat_dub$col[i] = "dodgerblue4"
  }else{dat_dub$col[i] = "grey"} # None significant
}


# Cov x GxE Plot
cge1 = filter(dat_csv, total_samples == 128)
cge1$n_pop_factor = NULL
for(i in 1:nrow(cge1)){
  if(cge1$n_pop[i] == 4){cge1$n_pop_factor[i] = "4 Genotypes"}
  else{cge1$n_pop_factor[i] = "8 Genotypes"}
}

(cge1_plot = ggplot(cge1, aes(x = covariance, y = GxE_emm, group = n_pop_factor, alpha = 0.1,colour = col)) + 
    geom_point() + ylim(0,1) + xlim(-1,1)+
    theme_classic(base_size = 24, base_family = "Times")+ 
    ylab(expression(bar(Delta)*""["GxE"]))+xlab(expression("Cov"["GE"]))+
    theme(legend.position = "none")+
    scale_colour_identity()+
    theme(axis.text.x=element_text(colour="black"))+
    theme(axis.text.y=element_text(colour="black"))+
    facet_wrap(~n_pop_factor)+
    ggtitle("Full Reciprocal Transplant Design\nTotal Samples = 128")+
    theme(plot.title = element_text(size = 24, face = "bold")))

cge2 = filter(dat_dub, total_samples == 128)
cge2$n_pop_factor = NULL
for(i in 1:nrow(cge2)){
  if(cge2$n_pop[i] == 4){cge2$n_pop_factor[i] = "2 Gen. per Env."
  }else if(cge2$n_pop[i] == 8){cge2$n_pop_factor[i] = "4 Gen. per Env."
  }else{cge2$n_pop_factor[i] = "8 Gen. per Env."}
}

(cge2_plot = ggplot(cge2, aes(x = covariance, y = GxE_emm, group = n_pop_factor, alpha = 0.1,colour = col)) + 
    geom_point() + ylim(0,1) + xlim(-1,1)+
    theme_classic(base_size = 24, base_family = "Times")+ 
    ylab(expression(bar(Delta)*""["GxE"]))+xlab(expression("Cov"["GE"]))+  theme(legend.position = "none")+
    scale_colour_identity()+
    theme(axis.text.x=element_text(colour="black"))+
    theme(axis.text.y=element_text(colour="black"))+
    facet_wrap(~n_pop_factor)+
    ggtitle("Paired Common Garden Design\nTotal Samples = 128")+
    theme(plot.title = element_text(size = 24, face = "bold")))

grid.arrange(cge1_plot, cge2_plot)

# Cov x Omega2 Plot

ggplot(dat_csv, aes(x = covariance, y = GxE_omega, group = factor(n_pop), alpha = 0.1,colour = col)) + 
  geom_point() + theme_classic() + 
  xlab("Covariance Estimate") + ylab("GxE Estimate (Omega^2)") +
  theme(legend.position = "none")+
  scale_colour_identity()+
  facet_grid(sample_size~n_pop)

ggplot((filter(dat_dub, sample_size != 16)), aes(x = covariance, y = GxE_omega, group = factor(n_pop), alpha = 0.1,colour = col)) + 
  geom_point() + theme_classic() + 
  xlab("Covariance Estimate") + ylab("GxE Estimate (Omega^2)") +
  theme(legend.position = "none")+
  scale_colour_identity()+
  facet_grid(sample_size~n_pop)


# GxE - FRT
gxe_hm = fnr.effsize(dat_csv1, metric = "gxe",  data.type = "raw",analysis = "perm", resolution = "fine")
gxe_hm1 = gxe_hm %>%
  filter(between(bin,0.3,0.6))%>%
  group_by(sample_size,n_pop) %>%
  summarize("avgpower" = mean(power))
gxe_hm1$totals = gxe_hm1$n_pop *gxe_hm1$n_pop * gxe_hm1$sample_size

# GxE - CG
gxe_hm2 = fnr.effsize(dat_dub1, metric = "gxe",  data.type = "raw",analysis = "perm", resolution = "fine")
gxe_hm3 = gxe_hm2 %>%
  filter(between(bin,0.3,0.6))%>%
  group_by(sample_size,n_pop) %>%
  summarize("avgpower" = mean(power))
gxe_hm3$totals = gxe_hm3$n_pop * 2 * gxe_hm3$sample_size

# HeatMaps - Covariance based on Permutation - FRT
(GxEPower1 = ggplot(gxe_hm1,aes(x = factor(sample_size), y = factor(n_pop), fill = avgpower)) + 
    geom_tile() + 
    geom_text(aes(label= paste(gxe_hm1$totals, round(gxe_hm1$avgpower,3),sep = '\n'),family = "Times"), size = 4) +
    theme_classic(base_size = 18, base_family = "Times")+ 
    scale_fill_viridis(
      breaks=seq(0,1,0.25), #breaks in the scale bar
      limits=c(0,1))+
    xlab("Sample Size") + ylab("Number of Genotypes")+
    labs(fill = "Power")+
    theme(axis.text = element_text(colour = "black"))+  
    theme(legend.position = "none")+
    ggtitle("FRT: Power")+
    theme(plot.title = element_text(size = 18, face = "bold")))

# HeatMaps - Covariance based on Permutation - FRT
(GxEPower2 = ggplot(gxe_hm3,aes(x = factor(sample_size), y = factor(n_pop), fill = avgpower)) + 
    geom_tile() + 
    geom_text(aes(label= paste(gxe_hm3$totals, round(gxe_hm3$avgpower,3),sep = '\n'),family = "Times"), size = 4) +
    theme_classic(base_size = 18, base_family = "Times")+ 
    scale_fill_viridis(
      breaks=seq(0,1,0.25), #breaks in the scale bar
      limits=c(0,1))+
    xlab("Sample Size") + ylab("Number of Genotypes")+
    labs(fill = "Power")+
    theme(axis.text = element_text(colour = "black"))+  
    theme(legend.position = "none")+
    ggtitle("CG: Power")+
    theme(plot.title = element_text(size = 18, face = "bold")))
# Full block (non-divided tiles) power analysis
cov_hm2 = fnr.effsize(dat_dub1, metric = "cov", data.type = "raw", analysis = "perm", resolution = "fine")
cov_hm3 = cov_hm2 %>%
  filter(between(bin,0.2,0.6))%>%
  group_by(sample_size,n_pop) %>%
  summarize("avgpower" = mean(power))
cov_hm3$totals = cov_hm3$n_pop * 2 * cov_hm3$sample_size

# HeatMaps - Covariance based on Permutation - FRT
(CovPower1 = ggplot(cov_hm1,aes(x = factor(sample_size), y = factor(n_pop), fill = avgpower)) + 
    geom_tile() + 
    geom_text(aes(label= paste(cov_hm1$totals, round(cov_hm1$avgpower,3),sep = '\n')), family = "Times", size = 5) +
    theme_classic(base_size = 18, base_family = "Times")+ 
    scale_fill_viridis(
      breaks=seq(0,1,0.25), #breaks in the scale bar
      limits=c(0,1))+
    xlab("Sample Size") + ylab("Number of Genotypes")+
    labs(fill = "Power")+
    theme(axis.text = element_text(colour = "black"))+
    theme(legend.position = "none")+
    ggtitle("FRT: Power")+
    theme(plot.title = element_text(size = 18, face = "bold")))

cov_hm = fnr.effsize(dat_csv1, metric = "cov", data.type = "raw", analysis = "perm", resolution = "fine")
cov_hm1 = cov_hm %>%
  filter(between(bin,0.2,0.6))%>%
  group_by(sample_size,n_pop) %>%
  summarize("avgpower" = mean(power))
cov_hm1$totals = cov_hm1$n_pop *cov_hm1$n_pop * cov_hm1$sample_size


# HeatMaps - Covariance based on Permutation - FRT
(CovPower2 = ggplot(cov_hm3,aes(x = factor(sample_size), y = factor(n_pop), fill = avgpower)) + 
    geom_tile() + 
    geom_text(aes(label= paste(cov_hm3$totals, round(cov_hm3$avgpower,3),sep = '\n')),family = "Times",size = 5) +
    theme_classic(base_size = 18, base_family = "Times")+ 
    scale_fil
  l_viridis(
      breaks=seq(0,1,0.25), #breaks in the scale bar
      limits=c(0,1))+
    xlab("Sample Size") + ylab("Number of Genotypes")+
    labs(fill = "Power")+
    theme(axis.text = element_text(colour = "black"))+
    theme(legend.position = "none")+
    ggtitle("CG: Power")+
    theme(plot.title = element_text(size = 18, face = "bold")))
#theme(legend.position = c(0.85,0.85)))
## False Positive Rates 
(falsePos2 = ggplot(dub_fpdf, aes(x = ID, y = rate, group = sample_size, fill = factor(sample_size)))+ 
   geom_bar(position = "dodge", stat = "identity") + 
   geom_hline(aes(yintercept = 0.05),linetype = "dashed")+
   ylab("False Positive Rate") + xlab("")+
   ggtitle("False Positive Rates: Paired Common Garden")+
   facet_wrap(~npop_plot) + 
   labs(fill = "Sample Size")+
   #ylim(0,0.3)+
   scale_fill_viridis(discrete = TRUE)+
   scale_x_discrete(labels=c("Cov_Perm" = "Perm. \n CovGE", 
                             "Cov_Boot" = "Boot. \n CovGE",
                             "GxE_Boot" = "Boot. \n GxE",
                             "GxE_Perm" = "Perm. \n GxE",
                             "GxE_Anova" = "ANOVA \n GxE"))+
   theme_classic(base_family = "Times",base_size = 16) + 
   theme(axis.text = element_text(colour = "black")))

grid.arrange(falsePos, falsePos2)

## False Negative Rates 
gxePow2$totals = gxePow2$n_pop * 2 * gxePow2$sample_size
(gxePow3 = ggplot(gxePow2,aes(x = factor(sample_size), y = power, group = bin,fill = bin)) + 
    geom_bar(position = "dodge", stat = "identity")+
    #geom_text(aes(label= paste(gxePow$totals, round(gxePow$power,3),sep = '\n')), size = 5) +
    theme_classic(base_size = 20, base_family = "Times")+ 
    scale_fill_viridis(discrete = TRUE)+
    #breaks=seq(0,1,0.25), #breaks in the scale bar
    #limits=c(0,1))+
    xlab("Sample Size") + ylab(expression("Power (1-"*beta*")"))+
    labs(fill = "Effect Size")+
    theme(axis.text = element_text(colour = "black"))+
    #theme(legend.position = "none")+
    ggtitle(expression("Common Garden: "*bar(Delta)*""["GxE"]*" Power"))+
    theme(plot.title = element_text(size = 24, face = "bold"))+
    facet_grid(~label*n_pop))

covPow2$totals = covPow2$n_pop * 2 * covPow2$sample_size
(covPow3 = ggplot(covPow2,aes(x = factor(sample_size), y = power, group = bin, fill = bin)) + 
    geom_bar(position = "dodge", stat = "identity")+
    theme_classic(base_size = 20, base_family = "Times")+ 
    scale_fill_viridis(discrete = TRUE)+
    xlab("Sample Size") + ylab(expression("Power (1-"*beta*")"))+
    labs(fill = "Effect Size")+
    theme(axis.text = element_text(colour = "black"))+
    ggtitle(expression("Common Garden: Cov"["GE"]*" Power"))+
    theme(plot.title = element_text(size = 24, face = "bold"))+
    facet_grid(~label*n_pop))
grid.arrange(covPow1, gxePow1, covPow3,gxePow3,ncol = 1)

## False Negative Rates 
(falseNeg = ggplot(fndf, aes(x = label, y = fnr, group = sample_size, fill = factor(sample_size)))+ 
   geom_bar(position = "dodge", stat = "identity") + 
   #geom_hline(aes(yintercept = 0.05),linetype = "dashed")+
   ylab("False Negative Rate") + xlab("")+
   ggtitle("Full Reciprocal Transplant")+
   facet_wrap(~n_pop) + 
   labs(fill = "Sample Size")+
   scale_fill_viridis(discrete = TRUE)+
   scale_x_discrete(labels=c("Cov_Perm" = "Perm. \n CovGE", 
                             "Cov_Boot" = "Boot. \n CovGE",
                             "GxE_Boot" = "Boot. \n GxE",
                             "GxE_Perm" = "Perm. \n GxE",
                             "GxE_Anova" = "ANOVA \n GxE"))+
   theme_classic(base_family = "Times",base_size = 16) + 
   theme(axis.text = element_text(colour = "black")))
bigcov = c("4_8_FRT","16_4_FRT","8_16_CG","16_8_CG")

## False Positive BarPlot for Covariance
label128 = c("4_16_CG" = "Common Garden \n 4 Samples \n 16 Genotypes \n 2 Environments", 
             "8_8_CG" = "Common Garden \n 8 Samples \n 8 Genotypes \n 2 Environments",
             "16_4_CG" = "Common Garden \n 16 Samples \n 4 Genotypes \n 2 Environments",
             "2_8_FRT" = "Reciprocal Transplant \n 2 Samples \n 8 Genotypes \n 8 Environments",
             "8_4_FRT" = "Reciprocal Transplant \n 8 Samples \n 4 Genotypes \n 4 Environments")
label128 = c("4_8_FRT" = "Reciprocal Transplant \n 4 Samples \n 8 Genotypes \n 8 Environments",
             "16_4_FRT"= "Reciprocal Transplant \n 16 Samples \n 4 Genotypes \n 4 Environments",
             "8_16_CG"= "Common Garden \n 8 Samples \n 16 Genotypes \n 2 Environments",
             "16_8_CG"="Common Garden \n 16 Samples \n 8 Genotypes \n 2 Environments")
dub_cov$scenario = rep("CG",nrow(dub_cov))
covFPR$scenario = rep("FRT",nrow(covFPR))
covbo = rbind(dub_cov,covFPR)
covbo$grp = paste(covbo$sample_size,covbo$n_pop,covbo$scenario,sep="_")
covbo = filter(covbo,totsamp == 128)
#covbo = filter(covbo, grp %in% bigcov)
covbo[is.nan(covbo)] <- 0

(falsePosCov = ggplot(covbo, aes(x = reorder(factor(ID),-totsamp), y = rate, group = grp,colour = factor(grp),fill = factor(grp)))+ 
    geom_bar(position = "dodge", stat = "identity") + 
    geom_hline(aes(yintercept = 0.05),linetype = "dashed")+
    ylab("False Positive Rate") + xlab("")+
    ggtitle("False Positive Rates")+
    labs(colour = "Experimental Design",fill = "Experimental Design")+
    scale_fill_viridis(discrete = TRUE, labels  = label128)+
    #labels = c("Cov_Boot" = "Bootstrap", "Cov_Perm" = "Permutation"))+
    scale_colour_viridis(discrete = TRUE, labels  = label128)+
    #labels = c("Cov_Boot" = "Bootstrap", "Cov_Perm" = "Permutation"))+
    scale_x_discrete(labels= c("Cov_Boot" = "Bootstrap", "Cov_Perm" = "Permutation"))+ #label128)+#
    theme(axis.text.x = element_text(colour = "black",size = 22)) +
    theme_classic(base_family = "Times",base_size = 14) + 
    theme(axis.text = element_text(colour = "black")))

## False Positive BarPlot for GxE 

bigpow = c("8_8_CG","16_8_CG","16_4_CG","8_16_CG","4_16_CG","16_4_FRT","8_4_FRT","4_8_FRT","2_8_FRT")
label128 = c("4_16_CG" = "Common Garden \n 4 Samples \n 16 Genotypes \n 2 Environments", # Total = 128
             "8_8_CG" = "Common Garden \n 8 Samples \n 8 Genotypes \n 2 Environments",
             "16_4_CG" = "Common Garden \n 16 Samples \n 4 Genotypes \n 2 Environments",
             "2_8_FRT" = "Reciprocal Transplant \n 2 Samples \n 8 Genotypes \n 8 Environments",
             "8_4_FRT" = "Reciprocal Transplant \n 8 Samples \n 4 Genotypes \n 4 Environments")
label128 = c("8_8_CG" = "Common Garden \n 8 Samples \n 8 Genotypes \n 2 Environments", # Big Cov
             "16_8_CG" = "Common Garden \n 16 Samples \n 8 Genotypes \n 2 Environments",
             "16_4_CG" = "Common Garden \n 16 Samples \n 4 Genotypes \n 2 Environments",
             "8_16_CG" = "Common Garden \n 8 Samples \n 16 Genotypes \n 2 Environments",
             "4_16_CG" = "Common Garden \n 4 Samples \n 16 Genotypes \n 2 Environments",
             "16_4_FRT" = "Reciprocal Transplant \n 16 Samples \n 4 Genotypes \n 4 Environments",
             "8_4_FRT" = "Reciprocal Transplant \n 8 Samples \n 4 Genotypes \n 4 Environments",
             "4_8_FRT" = "Reciprocal Transplant \n 4 Samples \n 8 Genotypes \n 8 Environments",
             "2_8_FRT" = "Reciprocal Transplant \n 2 Samples \n 8 Genotypes \n 8 Environments")

dub_gxe$scenario = rep("CG",nrow(dub_gxe))
gxeFPR$scenario = rep("FRT",nrow(gxeFPR))
combo = rbind(dub_gxe,gxeFPR)
combo$grp = paste(combo$sample_size,combo$n_pop,combo$scenario,sep="_")
combo = filter(combo,totsamp == 128)
#combo = filter(combo, grp %in% bigpow)
combo[is.nan(combo)] <- 0


(falsePosGxE = ggplot(combo, aes(x = reorder(factor(ID),-totsamp), y = rate, group = grp,colour = factor(grp),fill = factor(grp)))+ 
    geom_bar(position = "dodge", stat = "identity") + 
    geom_hline(aes(yintercept = 0.05),linetype = "dashed")+
    ylab("False Positive Rate") + xlab("")+
    ggtitle("False Positive Rates")+
    labs(colour = "Experimental Design",fill = "Experimental Design")+
    scale_fill_viridis(discrete = TRUE, labels  = label128)+
    #labels = c("GxE_Anova" = "ANOVA", "GxE_Perm" = "Permutation"))+
    scale_colour_viridis(discrete = TRUE,labels  = label128)+
    #labels = c("GxE_Anova" = "ANOVA", "GxE_Perm" = "Permutation"))+
    scale_x_discrete(labels= c("GxE_Anova" = "ANOVA", "GxE_Perm" = "Permutation"))+#label128)+
    theme(axis.text.x = element_text(colour = "black",size = 22)) +
    theme_classic(base_family = "Times",base_size = 14) + 
    theme(axis.text = element_text(colour = "black")))

## Old Power Analysis approach:


# Estimate Power and wrangle datums
covpow1 = cov1pow %>%
  group_by(sample_size, n_pop, std_dev) %>%
  summarize("total_tick" = n(),
            "covtick" = sum(covtick))
covpow1$covpower = covpow1$covtick/covpow1$total_tick
rng.cov1 = range(covpow1$covpower)
covpow1$totals = covpow1$sample_size* covpow1$n_pop * covpow1$n_pop
covpow1low = covpow1 %>%
  filter(std_dev == min(covpow1$std_dev))%>%
  filter(totals > 17)
covpow1hi = covpow1 %>%
  filter(std_dev == max(covpow1$std_dev))%>%
  filter(totals > 17)

gxepow1 = gxe1pow %>%
  group_by(sample_size,n_pop,std_dev) %>%
  summarize("total_tick" = n(),
            "gxetick" = sum(gxetick))
gxepow1$gxepower = gxepow1$gxetick/gxepow1$total_tick
gxepow1$totals = gxepow1$sample_size* gxepow1$n_pop * gxepow1$n_pop
rng.gxe1 = range(gxepow1$gxepower)
gxepow1low = gxepow1 %>%
  filter(std_dev == min(gxepow1$std_dev))%>%
  filter(totals > 17)
gxepow1hi = gxepow1%>%
  filter(std_dev == max(gxepow1$std_dev))%>%
  filter(totals > 17)

covpow2 = cov2pow %>%
  group_by(sample_size, n_pop, std_dev) %>%
  summarize("total_tick" = n(),
            "covtick" = sum(covtick))
covpow2$covpower = covpow2$covtick/covpow2$total_tick
covpow2$totals = covpow2$sample_size* covpow2$n_pop * 2
rng.cov2 = range(covpow2$covpower)
covpow2low = covpow2 %>% 
  filter(std_dev == min(covpow2$std_dev))%>%
  filter(totals > 17)
covpow2hi = covpow2 %>%
  filter(std_dev == max(covpow2$std_dev))%>%
  filter(totals > 17)

gxepow2 = gxe2pow %>%
  group_by(sample_size,n_pop,std_dev) %>%
  summarize("total_tick" = n(),
            "gxetick" = sum(gxetick))
gxepow2$gxepower = gxepow2$gxetick/gxepow2$total_tick
gxepow2$totals = gxepow2$sample_size* gxepow2$n_pop * 2
rng.gxe2 = range(gxepow2$gxepower)
gxepow2low = gxepow2 %>%
  filter(std_dev == min(gxepow2$std_dev))%>%
  filter(totals > 17)
gxepow2hi = gxepow2 %>%
  filter(std_dev == max(gxepow2$std_dev))%>%
  filter(totals > 17)

(covpower_high2 = ggplot(covpow2hi,aes(x = factor(sample_size), y = factor(n_pop), fill = covpower)) + 
    geom_tile() + 
    geom_text(aes(label= paste(covpow2hi$totals, round(covpow2hi$covpower,4),sep = '\n')), size = 5) +
    theme_classic(base_size = 24, base_family = "Times")+ 
    scale_fill_viridis(breaks=seq(0,1,0.25), #breaks in the scale bar
                       limits=c(0,1))+
    xlab("Sample Size") + ylab("Number of Genotypes")+
    labs(fill = "Power")+
    theme(axis.text = element_text(colour = "black"))+
    theme(legend.position = "none")+
    ggtitle(expression("Cov"["GE"]*": Paired Common Garden"))+
    theme(plot.title = element_text(size = 24, face = "bold")))

(gxepower_high1 = ggplot(gxepow1hi,aes(x = factor(sample_size), y = factor(n_pop), fill = gxepower)) + 
    geom_tile() + 
    geom_text(aes(label= paste(gxepow1hi$totals, round(gxepow1hi$gxepower,4),sep = '\n')), size = 5) +
    theme_classic(base_size = 24, base_family = "Times")+ 
    scale_fill_viridis(breaks=seq(0,1,0.25), #breaks in the scale bar
                       limits=c(0,1))+
    xlab("Sample Size") + ylab("Number of Genotypes")+
    labs(fill = "Power")+
    theme(axis.text = element_text(colour = "black"))+
    theme(legend.position = "none")+
    ggtitle(expression(bar(Delta)*""["GxE"]*": Full Reciprocal Transplant"))+
    theme(plot.title = element_text(size = 24, face = "bold")))

(gxepower_high2 = ggplot(gxepow2hi,aes(x = factor(sample_size), y = factor(n_pop), fill = gxepower)) + 
    geom_tile() + 
    geom_text(aes(label= paste(gxepow2hi$totals, round(gxepow2hi$gxepower,4),sep = '\n')), size = 5) +
    theme_classic(base_size = 24, base_family = "Times")+ 
    scale_fill_viridis(
      #scale_fill_gradient2(low="#DDDDDD", mid="#99CCEE", high="#000044", #colors in the scale
      #midpoint=mean(rng.gxe2),    #same midpoint for plots (mean of the range)
      breaks=seq(0,1,0.25), #breaks in the scale bar
      limits=c(0,1))+
    xlab("Sample Size") + ylab("Number of Genotypes")+
    labs(fill = "Power")+
    theme(axis.text = element_text(colour = "black"))+
    theme(legend.position = "none")+
    ggtitle(expression(bar(Delta)*""["GxE"]*": Paired Common Garden"))+
    theme(plot.title = element_text(size = 24, face = "bold")))

grid.arrange(gxepower_high1,gxepower_high2,covpower_high1,covpower_high2,ncol = 2)

### Low Std. Deviation Heatmaps -- Supplementary Material ###

(covpower_low1 = ggplot(covpow1low,aes(x = factor(sample_size), y = factor(n_pop), fill = covpower)) + 
    geom_tile() + 
    geom_text(aes(label= paste(covpow1low$totals, round(covpow1low$covpower,4),sep = '\n')), size = 5) +
    theme_classic(base_size = 24, base_family = "Times")+ 
    scale_fill_viridis(breaks=seq(0,1,0.25), #breaks in the scale bar
                       limits=c(0,1))+
    xlab("Sample Size") + ylab("Number of Genotypes")+
    labs(fill = "Power")+
    theme(axis.text = element_text(colour = "black"))+  
    theme(legend.position = "none")+
    ggtitle(expression("Cov"["GE"]*": Full Reciprocal Transplant"))+
    theme(plot.title = element_text(size = 24, face = "bold")))

(covpower_low2 = ggplot(covpow2low,aes(x = factor(sample_size), y = factor(n_pop), fill = covpower)) + 
    geom_tile() + 
    geom_text(aes(label= paste(covpow2low$totals, round(covpow2low$covpower,4),sep = '\n')), size = 5) +
    theme_classic(base_size = 24, base_family = "Times")+ 
    scale_fill_viridis(breaks=seq(0,1,0.25), #breaks in the scale bar
                       limits=c(0,1))+
    xlab("Sample Size") + ylab("Number of Genotypes")+
    labs(fill = "Power")+
    theme(axis.text = element_text(colour = "black"))+
    theme(legend.position = "none")+
    ggtitle(expression("Cov"["GE"]*": Paired Common Garden"))+
    theme(plot.title = element_text(size = 24, face = "bold")))

(gxepower_low1 = ggplot(gxepow1low,aes(x = factor(sample_size), y = factor(n_pop), fill = gxepower)) + 
    geom_tile() + 
    geom_text(aes(label= paste(gxepow1low$totals, round(gxepow1low$gxepower,4),sep = '\n')), size = 5) +
    theme_classic(base_size = 24, base_family = "Times")+ 
    scale_fill_viridis(breaks=seq(0,1,0.25), #breaks in the scale bar
                       limits=c(0,1))+
    xlab("Sample Size") + ylab("Number of Genotypes")+
    labs(fill = "Power")+
    theme(axis.text = element_text(colour = "black"))+
    theme(legend.position = "none")+
    ggtitle(expression(bar(Delta)*""["GxE"]*": Full Reciprocal Transplant"))+
    theme(plot.title = element_text(size = 24, face = "bold")))

(gxepower_low2 = ggplot(gxepow2low,aes(x = factor(sample_size), y = factor(n_pop), fill = gxepower)) + 
    geom_tile() + 
    geom_text(aes(label= paste(gxepow2low$totals, round(gxepow2low$gxepower,4),sep = '\n')), size = 5) +
    theme_classic(base_size = 24, base_family = "Times")+ 
    scale_fill_viridis(breaks=seq(0,1,0.25), #breaks in the scale bar
                       limits=c(0,1))+
    xlab("Sample Size") + ylab("Number of Genotypes")+
    labs(fill = "Power")+
    theme(axis.text = element_text(colour = "black"))+
    theme(legend.position = "none")+
    ggtitle(expression(bar(Delta)*""["GxE"]*": Paired Common Garden"))+
    theme(plot.title = element_text(size = 24, face = "bold")))

grid.arrange(gxepower_low1,gxepower_low2,covpower_low1,covpower_low2,ncol = 2)



# Bin Covariance Bootstrap and See whats driving false/true pos's and neg's
dat_csv1$binCov = "NA"
for(i in 1:nrow(dat_csv1)){
  if(dat_csv1$true_cov[i] == 0){dat_csv1$binCov[i] = 0
  }else if(abs(dat_csv1$true_cov[i]) > 0 & abs(dat_csv1$true_cov[i]) <= 0.15){dat_csv1$binCov[i] = 0.1
  }else if(abs(dat_csv1$true_cov[i]) > 0.15 & abs(dat_csv1$true_cov[i]) <= 0.25){dat_csv1$binCov[i] = 0.2
  }else if(abs(dat_csv1$true_cov[i]) > 0.25 & abs(dat_csv1$true_cov[i]) <= 0.35){dat_csv1$binCov[i] = 0.3
  }else if(abs(dat_csv1$true_cov[i]) > 0.35 & abs(dat_csv1$true_cov[i]) <= 0.45){dat_csv1$binCov[i] = 0.4
  }else if(abs(dat_csv1$true_cov[i]) > 0.45 & abs(dat_csv1$true_cov[i]) <= 0.55){dat_csv1$binCov[i] = 0.5
  }else if(abs(dat_csv1$true_cov[i]) > 0.55 & abs(dat_csv1$true_cov[i]) <= 0.65){dat_csv1$binCov[i] = 0.6
  }else if(abs(dat_csv1$true_cov[i]) > 0.65 & abs(dat_csv1$true_cov[i]) <= 0.75){dat_csv1$binCov[i] = 0.7
  }else if(abs(dat_csv1$true_cov[i]) > 0.75 & abs(dat_csv1$true_cov[i]) <= 0.85){dat_csv1$binCov[i] = 0.8
  }else if(abs(dat_csv1$true_cov[i]) > 0.85 & abs(dat_csv1$true_cov[i]) <= 0.95){dat_csv1$binCov[i] = 0.9
  }else{dat_csv1$binCov[i] = 1}
}

bindf = dat_csv1 %>%
  group_by(binCov)%>%
  summarize("N_total" = n())
bindf$binCov[11]<-"0"

bindf2 = dat_csv1 %>%
  group_by(binCov,Covconfintboot)%>%
  summarize("N_category" = n())
bindf2$binCov[17] <- "0"
bindf3 = merge(bindf, bindf2)

(PropCovBoot = ggplot(bindf3, aes(x = binCov, y = (N_category/N_total), group = Covconfintboot, fill = Covconfintboot))+
    geom_bar(stat="identity") + theme_classic()+ labs(fill = "")+
    scale_fill_manual(values = c("True Positive" = "coral2", "False Positive" = "darkolivegreen4","True Negative" = "darkolivegreen4" ,"False Negative" = "deepskyblue4"))+
    labs(colour = "") +xlab("Binned Covariance")+ylab("Category")+ggtitle("Covariance Bootstrap"))

(CovBoot = ggplot(bindf3, aes(x = abs(binCov), y = N_category, group = Covconfintboot, colour = Covconfintboot))+
    geom_point() + theme_hc()+ geom_line()+
    scale_colour_manual(values = c("True Positive" = "#992266", "False Positive" = "#009988","True Negative" = "#CC8800" ,"False Negative" = "#0000AA"))+
    labs(colour = "") +xlab("Binned Covariance")+ylab("Category")+ggtitle("Covariance Bootstrap"))

## CovGE -- Permutation -- Scen2
(cov_perm_dub = ggplot(transform(dat_dub, Covconfintperm = factor(Covconfintperm, levels = c("True Positive", "True Negative", "False Positive", "False Negative")))) +
    geom_errorbar(aes(x = reorder(row,true_cov), ymin = covariance_lwrCI, ymax = covariance_uprCI),alpha = 0.5)+ 
    geom_point(aes(x = reorder(row,true_cov), y = covariance), shape = 24,colour = "black", fill = "#7570B3",size = 1.2)+
    geom_point(aes(x = reorder(row,true_cov), y = true_cov), shape = 21, colour = "black", fill = "#E6AB02",size = 1.2)+
    theme_classic() + geom_hline(aes(yintercept = 0))+ ggtitle(expression("Cov"["GE"]*": Permutation"))+
    ylab(expression("Cov"["GE"]))+xlab(expression("Actual Cov"["GE"]))+
    facet_wrap(~Covconfintperm,ncol = 2)+
    theme_classic(base_family = "Times")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()))

(cov_confusion_perm_dub = ggplot(filter(dat_dub, replicate == 1), aes(x = true_cov, y = covariance)) + 
    geom_abline(slope = 1, intercept = 0,colour = "black")+
    geom_point(aes(fill = factor(Covconfintperm), shape = factor(Covconfintperm)),size = 3) +
    theme_classic(base_size = 18,base_family = "Times")+
    scale_shape_manual(values = c("False Negative" = 24, "True Positive"=21,"True Negative"= 22,"False Positive"=23)) +   
    scale_fill_manual(values = c("False Negative"="#238A8DFF", "True Positive" ="#FDE725FF", "True Negative"="#481567FF","False Positive"="#55C667FF")) +
    ylab(expression("Cov"["GE"]))+  
    xlab(expression("Actual Cov"["GE"]))+
    labs(fill = "", shape = "")+
    ggtitle(expression("CG Design: Cov"["GE"]*" Permutation"))+
    theme(legend.position = "none")+
    theme(plot.title = element_text(size = 18,face = "bold"))+
    theme(axis.text.y=element_text(colour = "black"))+
    theme(axis.text.x=element_text(colour = "black")))

## CovGE -- Bootstrap -- Scen2

(cov_confusion_boot_dub = ggplot(filter(dat_dub, replicate == 1), aes(x = true_cov, y = covariance)) + 
    geom_abline(slope = 1, intercept = 0,colour = "black")+
    geom_point(aes(fill = factor(Covconfintboot), shape = factor(Covconfintboot)),size = 3) +
    theme_classic(base_size = 18,base_family = "Times")+
    scale_shape_manual(values = c("False Negative" = 24, "True Positive"=21,"True Negative"= 22,"False Positive"=23)) +   
    scale_fill_manual(values = c("False Negative"="#238A8DFF", "True Positive" ="#FDE725FF", "True Negative"="#481567FF","False Positive"="#55C667FF")) +
    ylab(expression("Cov"["GE"]))+  
    xlab(expression("Actual Cov"["GE"]))+
    labs(fill = "", shape = "")+
    ggtitle(expression("CG Design: Cov"["GE"]*" Bootstrap"))+
    theme(plot.title = element_text(size = 18,face = "bold"))+
    theme(legend.position = "none")+
    theme(axis.text.y=element_text(colour = "black"))+
    theme(axis.text.x=element_text(colour = "black")))

(cov_boot_dub = ggplot(transform(dat_dub, Covconfintboot = factor(Covconfintboot, levels = c("True Positive", "True Negative", "False Positive", "False Negative")))) +
    geom_errorbar(aes(x = reorder(row,true_cov), ymin = covariance_lwrCI, ymax = covariance_uprCI),alpha = 0.5)+ 
    geom_point(aes(x = reorder(row,true_cov), y = covariance), shape = 24,colour = "black", fill = "#7570B3",size = 1.2)+
    geom_point(aes(x = reorder(row,true_cov), y = true_cov), shape = 21, colour = "black", fill = "#E6AB02",size = 1.2)+
    theme_classic() + geom_hline(aes(yintercept = 0))+ ggtitle(expression("Cov"["GE"]*": Bootstrap"))+
    ylab(expression("Cov"["GE"]))+xlab(expression("Actual Cov"["GE"]))+
    facet_wrap(~Covconfintboot,ncol = 2)+
    theme_classic(base_family = "Times")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()))

## GxE -- Permutation -- Scen2

(gxe_perm_dub = ggplot(transform(dat_dub, GxEconfintperm = factor(GxEconfintperm, levels = c("True Positive", "True Negative", "False Positive", "False Negative")))) +
    geom_errorbar(aes(x = reorder(row,true_GxE_emm), ymin = GxE_emm_lwrCI, ymax = GxE_emm_uprCI),alpha = 0.5)+ 
    geom_point(aes(x = reorder(row,true_GxE_emm), y = GxE_emm), shape = 24,colour = "black", fill = "#7570B3",size = 1.2)+
    geom_point(aes(x = reorder(row,true_GxE_emm), y = true_GxE_emm), shape = 21, colour = "black", fill = "#E6AB02",size = 1.2)+
    theme_classic() + geom_hline(aes(yintercept = 0))+ ggtitle(expression(bar(Delta)*""["GxE"]*": Permutation"))+
    ylab(expression(bar(Delta)*""["GxE"]))+
    facet_wrap(~GxEconfintperm,ncol = 2)+
    theme_classic(base_family = "Times")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()))

(gxe_confusion_perm_dub = ggplot(filter(dat_dub, replicate == 1), aes(x = true_GxE_emm, y = GxE_emm)) + 
    geom_abline(slope = 1, intercept = 0,colour = "black")+
    geom_point(aes(fill = factor(GxEconfintperm), shape = factor(GxEconfintperm)),size = 3) +
    theme_classic(base_size = 18,base_family = "Times")+
    scale_shape_manual(values = c("False Negative" = 24, "True Positive"=21,"True Negative"= 22,"False Positive"=23)) +   
    scale_fill_manual(values = c("False Negative"="#238A8DFF", "True Positive" ="#FDE725FF", "True Negative"="#481567FF","False Positive"="#55C667FF")) +
    ylab(expression(bar(Delta)*""["GxE"]))+
    xlab(expression("Actual "*bar(Delta)*""["GxE"]))+
    ggtitle(expression("CG Design: "*bar(Delta)*""["GxE"]*" Permutation"))+
    theme(legend.position = "none")+
    theme(plot.title = element_text(size = 18,face = "bold"))+
    theme(axis.text.y=element_text(colour = "black"))+
    theme(axis.text.x=element_text(colour = "black")))

## GxE -- Bootstrap -- Scen2

(gxe_boot_dub = ggplot(transform(dat_dub, GxEconfintboot = factor(GxEconfintboot, levels = c("True Positive", "True Negative", "False Positive", "False Negative")))) +
    geom_errorbar(aes(x = reorder(row,true_GxE_emm), ymin = GxE_emm_lwrCI, ymax = GxE_emm_uprCI),alpha = 0.5)+ 
    geom_point(aes(x = reorder(row,true_GxE_emm), y = GxE_emm), shape = 24,colour = "black", fill = "#7570B3",size = 1.2)+
    geom_point(aes(x = reorder(row,true_GxE_emm), y = true_GxE_emm), shape = 21, colour = "black", fill = "#E6AB02",size = 1.2)+
    theme_classic() + geom_hline(aes(yintercept = 0))+ ggtitle(expression(bar(Delta)*""["GxE"]*": Bootstrap"))+
    ylab(expression(bar(Delta)*""["GxE"]))+
    facet_wrap(~GxEconfintboot,ncol = 2)+
    theme_classic(base_family = "Times")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()))

(gxe_confusion_boot_dub = ggplot(filter(dat_dub, replicate == 1), aes(x = true_GxE_emm, y = GxE_emm)) + 
    geom_abline(slope = 1, intercept = 0,colour = "black")+
    geom_point(aes(fill = factor(GxEconfintboot), shape = factor(GxEconfintboot)),size = 3) +
    scale_shape_manual(values = c("False Negative" = 24, "True Positive"=21,"True Negative"= 22,"False Positive"=23)) +   
    scale_fill_manual(values = c("False Negative"="#238A8DFF", "True Positive" ="#FDE725FF", "True Negative"="#481567FF","False Positive"="#55C667FF")) +
    ylab(expression(bar(Delta)*""["GxE"]))+
    xlab(expression("Actual "*bar(Delta)*""["GxE"]))+
    ggtitle(expression("CG Design: "*bar(Delta)*""["GxE"]*" Bootstrap"))+
    theme_classic(base_size = 18, base_family = "Times")+
    theme(legend.position = "none")+
    theme(plot.title = element_text(size = 18,face = "bold"))+
    theme(axis.text.y=element_text(colour = "black"))+
    theme(axis.text.x=element_text(colour = "black")))

## MEANS -- CovGE -- Permutation -- Scen 2
(cov_perm_means_dub = ggplot(transform(dat_dub, meansCovconfintperm = factor(meansCovconfintperm, levels = c("True Positive", "True Negative", "False Positive", "False Negative")))) +
    geom_errorbar(aes(x = reorder(row,true_cov), ymin = covariance_lwrCI, ymax = covariance_uprCI),alpha = 0.5)+ 
    geom_point(aes(x = reorder(row,true_cov), y = covariance), shape = 24,colour = "black", fill = "#7570B3",size = 1.2)+
    geom_point(aes(x = reorder(row,true_cov), y = true_cov), shape = 21, colour = "black", fill = "#E6AB02",size = 1.2)+
    theme_classic() + geom_hline(aes(yintercept = 0))+ ggtitle(expression("Cov"["GE"]*": Permutation - Means"))+
    ylab(expression("Cov"["GE"]))+xlab(expression("Actual Cov"["GE"]))+
    facet_wrap(~meansCovconfintperm,ncol = 2)+
    theme_classic(base_family = "Times")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()))

(cov_confusion_perm_means_dub = ggplot(filter(dat_dub, replicate == 1), aes(x = true_cov_means, y = cov_means)) + 
    geom_abline(slope = 1, intercept = 0,colour = "black")+
    geom_point(aes(fill = factor(meansCovconfintperm), shape = factor(meansCovconfintperm)),size = 3) +
    theme_classic(base_size = 18,base_family = "Times")+
    scale_shape_manual(values = c("False Negative" = 24, "True Positive"=21,"True Negative"= 22,"False Positive"=23)) +   
    scale_fill_manual(values = c("False Negative"="#238A8DFF", "True Positive" ="#FDE725FF", "True Negative"="#481567FF","False Positive"="#55C667FF")) +
    ylab(expression("Cov"["GE"]*" from group means"))+  
    xlab(expression("Actual Cov"["GE"]*" from group means"))+
    labs(fill = "", shape = "")+
    ggtitle(expression("CG Design: Cov"["GE"]*" Permutation"))+
    theme(plot.title = element_text(size = 18,face = "bold"))+
    theme(legend.position = "none")+
    theme(axis.text.y=element_text(colour = "black"))+
    theme(axis.text.x=element_text(colour = "black")))

## MEANS -- CovGE -- Bootstrap -- Scen 2

(cov_boot_means_dub = ggplot(transform(dat_dub, MeansCovconfintboot = factor(MeansCovconfintboot, levels = c("True Positive", "True Negative", "False Positive", "False Negative")))) +
    geom_errorbar(aes(x = reorder(row,true_cov), ymin = covariance_lwrCI, ymax = covariance_uprCI),alpha = 0.5)+ 
    geom_point(aes(x = reorder(row,true_cov), y = covariance), shape = 24,colour = "black", fill = "#7570B3",size = 1.2)+
    geom_point(aes(x = reorder(row,true_cov), y = true_cov), shape = 21, colour = "black", fill = "#E6AB02",size = 1.2)+
    theme_classic() + geom_hline(aes(yintercept = 0))+ ggtitle(expression("Cov"["GE"]*": Bootstrap - Means"))+
    ylab(expression("Cov"["GE"]))+xlab(expression("Actual Cov"["GE"]))+
    facet_wrap(~MeansCovconfintboot,ncol = 2)+
    theme_classic(base_family = "Times")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()))

(cov_confusion_boot_means_dub = ggplot(filter(dat_dub, replicate == 1), aes(x = true_cov_means, y = cov_means)) + 
    geom_abline(slope = 1, intercept = 0,colour = "black")+
    geom_point(aes(fill = factor(MeansCovconfintboot), shape = factor(MeansCovconfintboot)),size = 3) +
    theme_classic(base_size = 18,base_family = "Times")+
    scale_shape_manual(values = c("False Negative" = 24, "True Positive"=21,"True Negative"= 22,"False Positive"=23)) +   
    scale_fill_manual(values = c("False Negative"="#238A8DFF", "True Positive" ="#FDE725FF", "True Negative"="#481567FF","False Positive"="#55C667FF")) +
    ylab(expression("Cov"["GE"]*" from group means"))+  
    xlab(expression("Actual Cov"["GE"]*" from group means"))+
    labs(fill = "", shape = "")+
    ggtitle(expression("CG Design: Cov"["GE"]*" Bootstrap"))+
    theme(plot.title = element_text(size = 18,face = "bold"))+
    theme(legend.position = "none")+
    theme(axis.text.y=element_text(colour = "black"))+
    theme(axis.text.x=element_text(colour = "black")))


## MEANS -- GxE -- Bootstrap -- Scen 2

(gxe_boot_means_dub = ggplot(transform(dat_dub, MeanGxEconfintboot = factor(MeanGxEconfintboot, levels = c("True Positive", "True Negative", "False Positive", "False Negative")))) +
    geom_errorbar(aes(x = reorder(row,true_GxE_emm), ymin = GxE_emm_lwrCI, ymax = GxE_emm_uprCI),alpha = 0.5)+ 
    geom_point(aes(x = reorder(row,true_GxE_emm), y = GxE_emm), shape = 24,colour = "black", fill = "#7570B3",size = 1.2)+
    geom_point(aes(x = reorder(row,true_GxE_emm), y = true_GxE_emm), shape = 21, colour = "black", fill = "#E6AB02",size = 1.2)+
    theme_classic() + geom_hline(aes(yintercept = 0))+ ggtitle(expression(bar(Delta)*""["GxE"]*": Bootstrap - Means"))+
    ylab(expression(bar(Delta)*""["GxE"]))+
    facet_wrap(~GxEconfintperm,ncol = 2)+
    theme_classic(base_family = "Times")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()))

(gxe_confusion_boot_means_dub = ggplot(filter(dat_dub, replicate == 1), aes(x = true_GxE_means, y = GxE_means)) + 
    geom_abline(slope = 1, intercept = 0,colour = "black")+
    geom_point(aes(fill = factor(MeanGxEconfintboot), shape = factor(MeanGxEconfintboot)),size = 3) +
    theme_classic(base_size = 18,base_family = "Times")+
    scale_shape_manual(values = c("False Negative" = 24, "True Positive"=21,"True Negative"= 22,"False Positive"=23)) +   
    scale_fill_manual(values = c("False Negative"="#238A8DFF", "True Positive" ="#FDE725FF", "True Negative"="#481567FF","False Positive"="#55C667FF")) +
    ylab(expression(bar(Delta)*""["GxE"]*" from group means"))+  
    xlab(expression("Actual "*bar(Delta)*""["GxE"]*" from group means"))+  
    labs(fill = "", shape = "")+
    ggtitle(expression("CG Design: "*bar(Delta)*""["GxE"]*" Bootstrap"))+
    theme(plot.title = element_text(size = 18,face = "bold"))+
    theme(legend.position = "none")+
    theme(axis.text.y=element_text(colour = "black"))+
    theme(axis.text.x=element_text(colour = "black")))

## MEANS -- GxE -- Permutation -- Scen 2

(means_gxe_perm_dub = ggplot(transform(dat_dub, meansGxEconfintperm = factor(meansGxEconfintperm, levels = c("True Positive", "True Negative", "False Positive", "False Negative")))) +
    geom_errorbar(aes(x = reorder(row,true_GxE_emm), ymin = GxE_emm_lwrCI, ymax = GxE_emm_uprCI),alpha = 0.5)+ 
    geom_point(aes(x = reorder(row,true_GxE_emm), y = GxE_emm), shape = 24,colour = "black", fill = "#7570B3",size = 1.2)+
    geom_point(aes(x = reorder(row,true_GxE_emm), y = true_GxE_emm), shape = 21, colour = "black", fill = "#E6AB02",size = 1.2)+
    geom_hline(aes(yintercept = 0))+ ggtitle(expression(bar(Delta)*""["GxE"]*": Permutation - Means"))+
    ylab(expression(bar(Delta)*""["GxE"]))+
    facet_wrap(~meansGxEconfintperm,ncol = 2)+
    theme_classic(base_family = "Times")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())) 

(gxe_confusion_perm_means_dub = ggplot(filter(dat_dub, replicate == 1), aes(x = true_GxE_means, y = GxE_means)) + 
    geom_abline(slope = 1, intercept = 0,colour = "black")+
    geom_point(aes(fill = factor(meansGxEconfintperm), shape = factor(meansGxEconfintperm)),size = 3) +
    theme_classic(base_size = 18,base_family = "Times")+
    scale_shape_manual(values = c("False Negative" = 24, "True Positive"=21,"True Negative"= 22,"False Positive"=23)) +   
    scale_fill_manual(values = c("False Negative"="#238A8DFF", "True Positive" ="#FDE725FF", "True Negative"="#481567FF","False Positive"="#55C667FF")) +
    ylab(expression(bar(Delta)*""["GxE"]*" from group means"))+  
    xlab(expression("Actual "*bar(Delta)*""["GxE"]*" from group means"))+  
    labs(fill = "", shape = "")+
    ggtitle(expression("CG Design: "*bar(Delta)*""["GxE"]*" Permutation"))+
    theme(plot.title = element_text(size = 18,face = "bold"))+
    theme(legend.position = "none")+
    theme(axis.text.y=element_text(colour = "black"))+
    theme(axis.text.x=element_text(colour = "black")))

## MEANS - CovGE -- Permutation
(cov_perm_means = ggplot(transform(dat_csv, meansCovconfintperm = factor(meansCovconfintperm, levels = c("True Positive", "True Negative", "False Positive", "False Negative")))) +
    geom_errorbar(aes(x = reorder(row,true_cov), ymin = covariance_lwrCI, ymax = covariance_uprCI),alpha = 0.5)+ 
    geom_point(aes(x = reorder(row,true_cov), y = covariance), shape = 24,colour = "black", fill = "#7570B3",size = 1.2)+
    geom_point(aes(x = reorder(row,true_cov), y = true_cov), shape = 21, colour = "black", fill = "#E6AB02",size = 1.2)+
    theme_classic() + geom_hline(aes(yintercept = 0))+ ggtitle(expression("Cov"["GE"]*": Permutation - Means"))+
    ylab(expression("Cov"["GE"]))+xlab(expression("Actual Cov"["GE"]))+
    facet_wrap(~meansCovconfintperm,ncol = 2)+
    theme_classic(base_family = "Times")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()))

(cov_confusion_perm_means = ggplot(filter(dat_csv, replicate == 1), aes(x = true_cov_means, y = cov_means)) + 
    geom_abline(slope = 1, intercept = 0,colour = "black")+
    geom_point(aes(fill = factor(meansCovconfintperm), shape = factor(meansCovconfintperm)),size = 3) +
    theme_classic(base_size = 18,base_family = "Times")+
    scale_shape_manual(values = c("False Negative" = 24, "True Positive"=21,"True Negative"= 22,"False Positive"=23)) +   
    scale_fill_manual(values = c("False Negative"="#238A8DFF", "True Positive" ="#FDE725FF", "True Negative"="#481567FF","False Positive"="#55C667FF")) +
    ylab(expression("Cov"["GE"]*" from group means"))+  
    xlab(expression("Actual Cov"["GE"]*" from group means"))+
    labs(fill = "", shape = "")+
    ggtitle(expression("RT Design: Cov"["GE"]*" Permutation"))+
    theme(plot.title = element_text(size = 18,face = "bold"))+
    theme(legend.position = "none")+
    theme(axis.text.y=element_text(colour = "black"))+
    theme(axis.text.x=element_text(colour = "black")))

## MEANS - CovGE -- Bootstrap

(cov_boot_means = ggplot(transform(dat_csv, MeansCovconfintboot = factor(MeansCovconfintboot, levels = c("True Positive", "True Negative", "False Positive", "False Negative")))) +
    geom_errorbar(aes(x = reorder(row,true_cov), ymin = covariance_lwrCI, ymax = covariance_uprCI),alpha = 0.5)+ 
    geom_point(aes(x = reorder(row,true_cov), y = covariance), shape = 24,colour = "black", fill = "#7570B3",size = 1.2)+
    geom_point(aes(x = reorder(row,true_cov), y = true_cov), shape = 21, colour = "black", fill = "#E6AB02",size = 1.2)+
    theme_classic() + geom_hline(aes(yintercept = 0))+ ggtitle(expression("Cov"["GE"]*": Bootstrap - Means"))+
    ylab(expression("Cov"["GE"]))+xlab(expression("Actual Cov"["GE"]))+
    facet_wrap(~MeansCovconfintboot,ncol = 2)+
    theme_classic(base_family = "Times")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()))

(cov_confusion_boot_means = ggplot(filter(dat_csv, replicate == 1), aes(x = true_cov_means, y = cov_means)) + 
    geom_abline(slope = 1, intercept = 0,colour = "black")+
    geom_point(aes(fill = factor(MeansCovconfintboot), shape = factor(MeansCovconfintboot)),size = 3) +
    theme_classic(base_size = 18,base_family = "Times")+
    scale_shape_manual(values = c("False Negative" = 24, "True Positive"=21,"True Negative"= 22,"False Positive"=23)) +   
    scale_fill_manual(values = c("False Negative"="#238A8DFF", "True Positive" ="#FDE725FF", "True Negative"="#481567FF","False Positive"="#55C667FF")) +
    ylab(expression("Cov"["GE"]*" from group means"))+  
    xlab(expression("Actual Cov"["GE"]*" from group means"))+
    labs(fill = "", shape = "")+
    ggtitle(expression("RT Design: Cov"["GE"]*" Bootstrap"))+
    theme(plot.title = element_text(size = 18,face = "bold"))+
    theme(legend.position = "none")+
    theme(axis.text.y=element_text(colour = "black"))+
    theme(axis.text.x=element_text(colour = "black")))

## MEANS - GxE -- Bootstrap

(gxe_boot_means = ggplot(transform(dat_csv, MeanGxEconfintboot = factor(MeanGxEconfintboot, levels = c("True Positive", "True Negative", "False Positive", "False Negative")))) +
    geom_errorbar(aes(x = reorder(row,true_GxE_emm), ymin = GxE_emm_lwrCI, ymax = GxE_emm_uprCI),alpha = 0.5)+ 
    geom_point(aes(x = reorder(row,true_GxE_emm), y = GxE_emm), shape = 24,colour = "black", fill = "#7570B3",size = 1.2)+
    geom_point(aes(x = reorder(row,true_GxE_emm), y = true_GxE_emm), shape = 21, colour = "black", fill = "#E6AB02",size = 1.2)+
    theme_classic() + geom_hline(aes(yintercept = 0))+ ggtitle(expression(bar(Delta)*""["GxE"]*": Bootstrap - Means"))+
    ylab(expression(bar(Delta)*""["GxE"]))+
    facet_wrap(~MeanGxEconfintboot,ncol = 2)+
    theme_classic(base_family = "Times")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()))

(gxe_confusion_boot_means = ggplot(filter(dat_csv, replicate == 1), aes(x = true_GxE_means, y = GxE_means)) + 
    geom_abline(slope = 1, intercept = 0,colour = "black")+
    geom_point(aes(fill = factor(MeanGxEconfintboot), shape = factor(MeanGxEconfintboot)),size = 3) +
    theme_classic(base_size = 18,base_family = "Times")+
    scale_shape_manual(values = c("False Negative" = 24, "True Positive"=21,"True Negative"= 22,"False Positive"=23)) +   
    scale_fill_manual(values = c("False Negative"="#238A8DFF", "True Positive" ="#FDE725FF", "True Negative"="#481567FF","False Positive"="#55C667FF")) +
    ylab(expression(bar(Delta)*""["GxE"]*" from group means"))+  
    xlab(expression("Actual "*bar(Delta)*""["GxE"]*" from group means"))+  
    labs(fill = "", shape = "")+
    ggtitle(expression("RT Design: "*bar(Delta)*""["GxE"]*" Bootstrap"))+
    theme(plot.title = element_text(size = 18,face = "bold"))+
    theme(legend.position = "none")+
    theme(axis.text.y=element_text(colour = "black"))+
    theme(axis.text.x=element_text(colour = "black")))

## MEANS - GxE -- Permutation

(means_gxe_perm = ggplot(transform(dat_csv, meansGxEconfintperm = factor(meansGxEconfintperm, levels = c("True Positive", "True Negative", "False Positive", "False Negative")))) +
    geom_errorbar(aes(x = reorder(row,true_GxE_emm), ymin = GxE_emm_lwrCI, ymax = GxE_emm_uprCI),alpha = 0.5)+ 
    geom_point(aes(x = reorder(row,true_GxE_emm), y = GxE_emm), shape = 24,colour = "black", fill = "#7570B3",size = 1.2)+
    geom_point(aes(x = reorder(row,true_GxE_emm), y = true_GxE_emm), shape = 21, colour = "black", fill = "#E6AB02",size = 1.2)+
    theme_classic() + geom_hline(aes(yintercept = 0))+ ggtitle(expression(bar(Delta)*""["GxE"]*": Permutation - Means"))+
    ylab(expression(bar(Delta)*""["GxE"]))+
    facet_wrap(~meansGxEconfintperm,ncol = 2)+
    theme_classic(base_family = "Times")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()))  

(gxe_confusion_perm_means = ggplot(filter(dat_csv, replicate == 1), aes(x = true_GxE_means, y = GxE_means)) + 
    geom_abline(slope = 1, intercept = 0,colour = "black")+
    geom_point(aes(fill = factor(meansGxEconfintperm), shape = factor(meansGxEconfintperm)),size = 3) +
    theme_classic(base_size = 18,base_family = "Times")+
    scale_shape_manual(values = c("False Negative" = 24, "True Positive"=21,"True Negative"= 22,"False Positive"=23)) +   
    scale_fill_manual(values = c("False Negative"="#238A8DFF", "True Positive" ="#FDE725FF", "True Negative"="#481567FF","False Positive"="#55C667FF")) +
    ylab(expression(bar(Delta)*""["GxE"]*" from group means"))+  
    xlab(expression("Actual "*bar(Delta)*""["GxE"]*" from group means"))+  
    labs(fill = "", shape = "")+
    ggtitle(expression("RT Design: "*bar(Delta)*""["GxE"]*" Permutation"))+
    theme(plot.title = element_text(size = 18,face = "bold"))+
    theme(legend.position = "none")+
    theme(axis.text.y=element_text(colour = "black"))+
    theme(axis.text.x=element_text(colour = "black")))

## CovGE -- Permutation
(cov_perm = ggplot(transform(dat_csv1, Covconfintperm = factor(Covconfintperm, levels = c("True Positive", "True Negative", "False Positive", "False Negative")))) +
    geom_errorbar(aes(x = reorder(row,true_cov), ymin = covariance_lwrCI, ymax = covariance_uprCI),alpha = 0.5)+ 
    geom_point(aes(x = reorder(row,true_cov), y = covariance), shape = 24,colour = "black", fill = "#7570B3",size = 1.2)+
    geom_point(aes(x = reorder(row,true_cov), y = true_cov), shape = 21, colour = "black", fill = "#E6AB02",size = 1.2)+
    theme_classic() + geom_hline(aes(yintercept = 0))+ ggtitle(expression("Cov"["GE"]*": Permutation"))+
    ylab(expression("Cov"["GE"]))+xlab(expression("Actual Cov"["GE"]))+
    facet_wrap(~Covconfintperm,ncol = 2)+
    theme_classic(base_family = "Times")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_text(colour = "black"))+
    theme(plot.title = element_text(size = 14, face = "bold")))

(cov_confusion_perm = ggplot(filter(dat_csv1, replicate == 1), aes(x = true_cov, y = covariance)) + 
    geom_abline(slope = 1, intercept = 0,colour = "black")+
    geom_point(aes(fill = factor(Covconfintperm), shape = factor(Covconfintperm)),size = 3) +
    theme_classic(base_size = 18,base_family = "Times")+
    scale_shape_manual(values = c("False Negative" = 24, "True Positive"=21,"True Negative"= 22,"False Positive"=23)) +   
    scale_fill_manual(values = c("False Negative"="#238A8DFF", "True Positive" ="#FDE725FF", "True Negative"="#481567FF","False Positive"="#55C667FF")) +
    ylab(expression("Cov"["GE"]))+  
    xlab(expression("Actual Cov"["GE"]))+
    labs(fill = "", shape = "")+
    ggtitle(expression("RT Design: Cov"["GE"]*" Permutation"))+
    theme(plot.title = element_text(size = 18,face = "bold"))+
    theme(legend.position = "none")+
    theme(axis.text.y=element_text(colour = "black"))+
    theme(axis.text.x=element_text(colour = "black")))


## CovGE  -- Bootstrap
(cov_confusion_boot = ggplot(filter(dat_csv1, replicate == 1), aes(x = true_cov, y = covariance)) + 
    geom_abline(slope = 1, intercept = 0,colour = "black")+
    geom_point(aes(fill = factor(Covconfintboot), shape = factor(Covconfintboot)),size = 3) +
    theme_classic(base_size = 18,base_family = "Times")+
    scale_shape_manual(values = c("False Negative" = 24, "True Positive"=21,"True Negative"= 22,"False Positive"=23)) +   
    scale_fill_manual(values = c("False Negative"="#238A8DFF", "True Positive" ="#FDE725FF", "True Negative"="#481567FF","False Positive"="#55C667FF")) +
    ylab(expression("Cov"["GE"]))+  
    xlab(expression("Actual Cov"["GE"]))+
    labs(fill = "", shape = "")+
    theme(plot.title = element_text(size = 18,face = "bold"))+
    ggtitle(expression("RT Design: Cov"["GE"]*" Bootstrap"))+
    theme(legend.position = "none")+
    theme(axis.text.y=element_text(colour = "black"))+
    theme(axis.text.x=element_text(colour = "black")))

(cov_boot = ggplot(transform(dat_csv1, Covconfintboot = factor(Covconfintboot, levels = c("True Positive", "True Negative", "False Positive", "False Negative")))) +
    geom_errorbar(aes(x = reorder(row,true_cov), ymin = covariance_lwrCI, ymax = covariance_uprCI),alpha = 0.5)+ 
    geom_point(aes(x = reorder(row,true_cov), y = covariance), shape = 24,colour = "black", fill = "#7570B3",size = 1.2)+
    geom_point(aes(x = reorder(row,true_cov), y = true_cov), shape = 21, colour = "black", fill = "#E6AB02",size = 1.2)+
    theme_classic() + geom_hline(aes(yintercept = 0))+ ggtitle(expression("Cov"["GE"]*": Bootstrap"))+
    ylab(expression("Cov"["GE"]))+xlab(expression("Actual Cov"["GE"]))+
    facet_wrap(~Covconfintboot,ncol = 2)+
    theme_classic(base_family = "Times")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()))

## GxE -- Permutation
(gxe_perm = ggplot(transform(dat_csv1, GxEconfintperm = factor(GxEconfintperm, levels = c("True Positive", "True Negative", "False Positive", "False Negative")))) +
    geom_errorbar(aes(x = reorder(row,true_GxE_emm), ymin = GxE_emm_lwrCI, ymax = GxE_emm_uprCI),alpha = 0.5)+ 
    geom_point(aes(x = reorder(row,true_GxE_emm), y = GxE_emm), shape = 24,colour = "black", fill = "#7570B3",size = 1.2)+
    geom_point(aes(x = reorder(row,true_GxE_emm), y = true_GxE_emm), shape = 21, colour = "black", fill = "#E6AB02",size = 1.2)+
    theme_classic() + geom_hline(aes(yintercept = 0))+ ggtitle(expression(bar(Delta)*""["GxE"]*": Permutation"))+
    ylab(expression(bar(Delta)*""["GxE"]))+
    facet_wrap(~GxEconfintperm,ncol = 2)+
    theme_classic(base_family = "Times")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()))

(gxe_confusion_perm = ggplot(filter(dat_csv1, replicate == 1), aes(x = true_GxE_emm, y = GxE_emm)) + 
    geom_abline(slope = 1, intercept = 0,colour = "black")+
    geom_point(aes(fill = factor(GxEconfintperm), shape = factor(GxEconfintperm)),size = 3) +
    theme_classic(base_size = 18,base_family = "Times")+
    scale_shape_manual(values = c("False Negative" = 24, "True Positive"=21,"True Negative"= 22,"False Positive"=23)) +   
    scale_fill_manual(values = c("False Negative"="#238A8DFF", "True Positive" ="#FDE725FF", "True Negative"="#481567FF","False Positive"="#55C667FF")) +
    ylab(expression(bar(Delta)*""["GxE"]))+
    xlab(expression("Actual "*bar(Delta)*""["GxE"]))+
    ggtitle(expression("RT Design: "*bar(Delta)*""["GxE"]*" Permutation"))+
    theme(legend.position = "none")+
    theme(plot.title = element_text(size = 18,face = "bold"))+
    theme(axis.text.y=element_text(colour = "black"))+
    theme(axis.text.x=element_text(colour = "black")))

## GxE -- Bootstrap
(gxe_boot = ggplot(transform(dat_csv1, GxEconfintboot = factor(GxEconfintboot, levels = c("True Positive", "True Negative", "False Positive", "False Negative")))) +
    geom_errorbar(aes(x = reorder(row,true_GxE_emm), ymin = GxE_emm_lwrCI, ymax = GxE_emm_uprCI),alpha = 0.5)+ 
    geom_point(aes(x = reorder(row,true_GxE_emm), y = GxE_emm), shape = 24,colour = "black", fill = "#7570B3",size = 1.2)+
    geom_point(aes(x = reorder(row,true_GxE_emm), y = true_GxE_emm), shape = 21, colour = "black", fill = "#E6AB02",size = 1.2)+
    theme_classic() + geom_hline(aes(yintercept = 0))+ ggtitle(expression(bar(Delta)*""["GxE"]*": Bootstrap"))+
    ylab(expression(bar(Delta)*""["GxE"]))+
    facet_wrap(~GxEconfintboot,ncol = 2)+
    theme_classic(base_family = "Times")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()))

(gxe_confusion_boot = ggplot(filter(dat_csv1, replicate == 1), aes(x = true_GxE_emm, y = GxE_emm)) + 
    geom_abline(slope = 1, intercept = 0,colour = "black")+
    geom_point(aes(fill = factor(GxEconfintboot), shape = factor(GxEconfintboot)),size = 3) +
    scale_shape_manual(values = c("False Negative" = 24, "True Positive"=21,"True Negative"= 22,"False Positive"=23)) +   
    scale_fill_manual(values = c("False Negative"="#238A8DFF", "True Positive" ="#FDE725FF", "True Negative"="#481567FF","False Positive"="#55C667FF")) +
    ylab(expression(bar(Delta)*""["GxE"]))+
    xlab(expression("Actual "*bar(Delta)*""["GxE"]))+
    ggtitle(expression("RT Design: "*bar(Delta)*""["GxE"]*" Bootstrap"))+
    theme_classic(base_size = 18, base_family = "Times")+
    theme(plot.title = element_text(size = 18,face = "bold"))+
    theme(legend.position = "none")+
    theme(axis.text.y=element_text(colour = "black"))+
    theme(axis.text.x=element_text(colour = "black")))




# Panels for Confusion Matrices for both scenarios 

# Confusion Plots for Raw (both scenarios)
grid.arrange(gxe_confusion_boot,gxe_confusion_boot_dub,
             gxe_confusion_perm,gxe_confusion_perm_dub,
             cov_confusion_boot,cov_confusion_boot_dub,
             cov_confusion_perm,cov_confusion_perm_dub,
             ncol = 2)

# Confusion Plots for Means (both scenarios)
grid.arrange(gxe_confusion_boot_means,gxe_confusion_boot_means_dub,
             gxe_confusion_perm_means,gxe_confusion_perm_means_dub,
             cov_confusion_boot_means,cov_confusion_boot_means_dub,
             cov_confusion_perm_means,cov_confusion_perm_means_dub,
             ncol = 2)

# 16 panel grid for Env Scenario 1
title1=textGrob("Full Reciprocal Transplant Design", gp=gpar(font= 2))
grid.arrange(gxe_perm,gxe_boot,cov_perm,cov_boot,top=title1)

# 16 panel grid for Env Scenario 1 - Means
title2=textGrob("Full Reciprocal Transplant Design - Means", gp=gpar(font= 2))
grid.arrange(means_gxe_perm, gxe_boot_means,cov_perm_means, cov_boot_means,top=title2)

# 16 panel grid for Env Scenario 2
title3=textGrob("Paired Common Garden Design", gp=gpar(font= 2))
grid.arrange(gxe_perm_dub,gxe_boot_dub,cov_perm_dub,cov_boot_dub,top = title3)

# 16 panel grid for Env Scenario 2 - Means
title4=textGrob("Paired Common Garden Design - Means", gp=gpar(font= 2))
grid.arrange(means_gxe_perm_dub, gxe_boot_means_dub,cov_perm_means_dub, cov_boot_means_dub, top=title4)



