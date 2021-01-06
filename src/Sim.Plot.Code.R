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
setwd("~/Documents/GitHub/CnGV/CnGV/results/Sim_12.15.20/")
start_df = read.csv("Power_output_results.csv") 
phen_data = read.csv("~/Desktop/phenotype_output_results.csv")

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
start_params = read.csv("~/Desktop/df.csv")
start_params2 = read.csv("~/Desktop/df2pop.csv")
(missing_rows = anti_join(start_df, start_params, by = "row"))
start_df = filter(start_df, row %notin% missing_rows$row)
#write.csv(start_df1[,-1], "~/Desktop/rerun.csv")
range(start_df$row)
args = filter(start_params, row == 1824)
args = args[,-1]

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
   ylab(expression("Actual "*bar(Delta)*""["GxE"]))+xlab(expression("Actual Cov"["GE"]))+
   ggtitle("Full Reciprocal Transplant") + facet_grid(ss_f~np_f) +
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
    ylab(expression("Actual "*bar(Delta)*""["GxE"]))+xlab(expression("Actual Cov"["GE"]))+
    ggtitle("Paired Common Garden") + facet_grid(ss_f~np_f) +
    theme_classic(base_family = "Times"))

grid.arrange(hexy,hexy2) 

######################
## Covariance x GxE ##
######################

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


########################################################
##          Confusion Matrices  -- Env Scenario 1     ##
########################################################

# Covariance Permutation check
dat_csv$Covconfintperm = rep("NA",nrow(dat_csv))
for(i in 1:nrow(dat_csv)){
    if(dat_csv$true_cov[i] != 0 && dat_csv$covariance_pvalue[i] <= 0.025){dat_csv$Covconfintperm[i] = "True Positive"
    }else if(dat_csv$true_cov[i] == 0 & dat_csv$covariance_pvalue[i] <= 0.025){dat_csv$Covconfintperm[i] = "False Positive"
    }else if(dat_csv$true_cov[i]!= 0 & dat_csv$covariance_pvalue[i] > 0.025){dat_csv$Covconfintperm[i] = "False Negative"
    }else if(dat_csv$true_cov[i] == 0 & dat_csv$covariance_pvalue[i] > 0.025){dat_csv$Covconfintperm[i] = "True Negative"
    }else{dat_csv$Covconfintperm[i] = "None"}
}

# Cov Boot check
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

# GxE Perm check
dat_csv$GxEconfintperm = rep("NA", nrow(dat_csv))
for(i in 1:nrow(dat_csv)){
  if(dat_csv$true_GxE_emm[i] != 0 & round(dat_csv$GxE_emm_pvalue[i],2) <= 0.05){dat_csv$GxEconfintperm[i] = "True Positive"
  }else if(dat_csv$true_GxE_emm[i] == 0 & round(dat_csv$GxE_emm_pvalue[i],2) <= 0.05){dat_csv$GxEconfintperm[i] = "False Positive"
  }else if(dat_csv$true_GxE_emm[i] != 0 & round(dat_csv$GxE_emm_pvalue[i],2) > 0.05){dat_csv$GxEconfintperm[i] = "False Negative"
  }else if(dat_csv$true_GxE_emm[i] == 0 & round(dat_csv$GxE_emm_pvalue[i],2) > 0.05){dat_csv$GxEconfintperm[i] = "True Negative"
  }else{dat_csv$GxEconfintperm == "None"}
}

# GxE Boot check
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

# GxE Anova check
dat_csv$GxEanova_conf = rep("NA", nrow(dat_csv))
for(i in 1:nrow(dat_csv)){
  if(dat_csv$true_GxE_emm[i] != 0 & round(dat_csv$GxE_Anova[i],2) <= 0.05){dat_csv$GxEanova_conf[i] = "True Positive"
  }else if(dat_csv$true_GxE_emm[i] == 0 & round(dat_csv$GxE_Anova[i],2) <= 0.05){dat_csv$GxEanova_conf[i] = "False Positive"
  }else if(dat_csv$true_GxE_emm[i] != 0 & round(dat_csv$GxE_Anova[i],2) > 0.05){dat_csv$GxEanova_conf[i] = "False Negative"
  }else if(dat_csv$true_GxE_emm[i] == 0 & round(dat_csv$GxE_Anova[i],2) > 0.05){dat_csv$GxEanova_conf[i] = "True Negative"
  }else{dat_csv$GxEanova_conf == "None"}
}

## Counts for table

dat_csv1 <- dat_csv %>% 
  filter(std_dev == 1) %>%
  filter(total_samples > 17)
  #filter(total_samples == 128)

dat_csv2 <- dat_csv %>% 
  filter(std_dev == 1)%>%
  filter(total_samples > 17)
  #filter(total_samples == 128)

## Tables
cov_perm_table = dat_csv1 %>%
    group_by("name" = Covconfintperm) %>%
    summarize("n" = n())
fpr.fnr(cov_perm_table, divided = FALSE, scenario = 1)

cov_boot_table = dat_csv1 %>%
    group_by("name" = Covconfintboot) %>%
  summarize("n" = n())
fpr.fnr(cov_boot_table, divided = FALSE, scenario = 1)

gxe_anova_table = dat_csv2 %>%
  group_by("name" =GxEanova_conf) %>%
  summarize("n" = n())
fpr.fnr(gxe_anova_table, divided = FALSE, scenario = 1)

gxe_perm_table = dat_csv2 %>%
    group_by("name" =GxEconfintperm) %>%
    summarize("n" = n())
fpr.fnr(gxe_perm_table, divided = FALSE, scenario = 1)

gxe_boot_table = dat_csv2 %>%
  group_by("name" =GxEconfintboot) %>%
  summarize("n" = n())
fpr.fnr(gxe_boot_table, divided = FALSE, scenario = 1)

## Counts for heatmaps
(raw_confusion_hmap1 = dat_csv1 %>%
    group_by(sample_size, n_pop, "name" =Covconfintperm) %>%
    summarize("n" = n()))
raw_conf1 = fpr.fnr(raw_confusion_hmap1, divided = TRUE, scenario = 1)
raw_conf_plot1 <- heatmap_fun(raw_conf1,"rate") # Can also do "percent"

(raw_confusion_hmap2 = dat_csv1 %>%
    group_by(sample_size, n_pop, "name" =Covconfintboot) %>%
    summarize("n" = n()))
raw_conf2 = fpr.fnr(raw_confusion_hmap2, divided = TRUE, scenario = 1)
raw_conf_plot2 <- heatmap_fun(raw_conf2,"rate")

(raw_confusion_hmap3 = dat_csv2 %>%
    group_by(sample_size, n_pop,"name" = GxEconfintboot) %>%
    summarize("n" = n()))
raw_conf3 = fpr.fnr(raw_confusion_hmap3, divided = TRUE, scenario = 1)
raw_conf_plot3 <- heatmap_fun(raw_conf3,"rate")

raw_confusion_hmap4 = dat_csv2 %>%
    group_by(sample_size, n_pop,"name" = GxEconfintperm) %>%
    summarize("n" = n())
raw_conf4 = fpr.fnr(raw_confusion_hmap4, divided = TRUE, scenario = 1)
raw_conf_plot4 <- heatmap_fun(raw_conf4,"rate")

raw_confusion_hmap5 = dat_csv2 %>%
  group_by(sample_size, n_pop,"name" = GxEanova_conf) %>%
  summarize("n" = n())
raw_conf5 = fpr.fnr(raw_confusion_hmap5, divided = TRUE, scenario = 1)
raw_conf_plot5 <- heatmap_fun(raw_conf5,"rate")

# Compile FPs for plot
raw_conf1$ID = rep("Cov_Perm", nrow(raw_conf1))
raw_conf2$ID = rep("Cov_Boot", nrow(raw_conf2))
raw_conf3$ID = rep("GxE_Boot", nrow(raw_conf3))
raw_conf4$ID = rep("GxE_Perm", nrow(raw_conf4))
raw_conf5$ID = rep("GxE_Anova", nrow(raw_conf5))
fpdf = rbind(raw_conf1,raw_conf2,raw_conf3,raw_conf4,raw_conf5)
gxeFPR = rbind(raw_conf4,raw_conf5)
gxeFPR = gxeFPR[gxeFPR$name == "False Positive",]

## False Positive Rates
fpdf$npop_plot = NA
for(i in 1:nrow(fpdf)){
  if(fpdf$n_pop[i] == 2){fpdf$npop_plot[i] = "2 Genotypes"
  }else if(fpdf$n_pop[i] == 4){fpdf$npop_plot[i] = "4 Genotypes"
  }else{fpdf$npop_plot[i] = "8 Genotypes"}
}
fpdf1 = fpdf[fpdf$name == "False Positive",]

## False Negative Rates
covperm1 = fnr.effsize(dat_csv1, metric = "cov", analysis ="perm",resolution = "fine")
covperm1$label = rep("Cov_Perm",nrow(covperm1))
covboot1 = fnr.effsize(dat_csv1, metric = "cov", analysis ="boot",resolution = "fine")
covboot1$label = rep("Cov_Boot",nrow(covboot1))
gxeperm1 = fnr.effsize(dat_csv1, metric = "gxe", analysis ="perm",resolution = "fine")
gxeperm1$label = rep("GxE_Perm",nrow(gxeperm1))
gxeboot1 = fnr.effsize(dat_csv1, metric = "gxe", analysis = "boot",resolution = "fine")
gxeboot1$label = rep("GxE_Boot",nrow(gxeboot1))
gxeanova1 = fnr.effsize(dat_csv1, metric = "gxe", analysis = "anova",resolution = "fine")
gxeanova1$label = rep("GxE_Anova",nrow(gxeanova1))

fndf = rbind(covperm1,covboot1,gxeperm1,gxeboot1,gxeanova1)
gxePow = rbind(gxeperm1,gxeanova1)
covPow = rbind(covperm1,covboot1)

############ Confusion Plots  -- Env Scenario 1 ###############

## False Positive BarPlot
dub_gxe$scenario = rep("CG",nrow(dub_gxe))
gxeFPR$scenario = rep("FRT",nrow(gxeFPR))
combo = rbind(dub_gxe,gxeFPR)

combo = filter(combo,totsamp == 128)
combo$grp = paste(combo$sample_size,combo$n_pop,combo$scenario,sep="_")
combo[is.nan(combo)] <- 0

(falsePos = ggplot(combo, aes(x = reorder(factor(grp),-rate), y = rate, group = ID,colour = factor(ID),fill = factor(ID)))+ 
  geom_bar(position = "dodge", stat = "identity") + 
  geom_hline(aes(yintercept = 0.05),linetype = "dashed")+
  ylab("False Positive Rate") + xlab("")+
  ggtitle("False Positive Rates")+
 # facet_wrap(~npop_plot) + 
  labs(colour = "Experimental Design",fill = "Experimental Design")+
  scale_fill_viridis(discrete = TRUE)+
    scale_colour_viridis(discrete = TRUE)+
  scale_x_discrete(labels=c("4_16_CG" = "Common Garden \n 4 Samples \n 16 Genotypes", 
                            "8_8_CG" = "Common Garden \n 8 Samples \n 8 Genotypes",
                            "16_4_CG" = "Common Garden \n 16 Samples \n 4 Genotypes",
                            "2_8_FRT" = "Reciprocal Transplant \n 2 Samples \n 8 Genotypes",
                            "8_4_FRT" = "Reciprocal Transplant \n 8 Samples \n 4 Genotypes"))+
  theme_classic(base_family = "Times",base_size = 16) + 
  theme(axis.text = element_text(colour = "black")))

## False Positive Heatmap 
(falsePos = ggplot(filter(raw_conf1, name == "False Positive"), aes(x = sample_size, y = n_pop, group = sample_size, fill = rate))+ 
  geom_tile() + 
  geom_text(aes(label= paste(total, round(rate,2),sep = '\n')), size = 5) +
  theme_classic(base_size = 24, base_family = "Times")+ 
  scale_fill_viridis(
    breaks=seq(0,1,0.25), #breaks in the scale bar
    limits=c(0,1))+
  xlab("Sample Size") + ylab("Number of Genotypes")+
  labs(fill = "Power")+
  theme(axis.text = element_text(colour = "black"))+  
  theme(legend.position = "none")+
  ggtitle(expression("Cov"["GE"]*": Full Reciprocal Transplant"))+
  theme(plot.title = element_text(size = 24, face = "bold")))

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


## GxE Power Bar plot
gxePow$totals = gxePow$n_pop * gxePow$n_pop * gxePow$sample_size
(gxePow1 = ggplot(gxePow,aes(x = factor(sample_size), y = power, group = bin,fill = bin)) + 
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
    ggtitle(expression("Reciprocal Transplant: "*bar(Delta)*""["GxE"]*" Power"))+    
    theme(plot.title = element_text(size = 24, face = "bold"))+
    facet_grid(~label*n_pop))

## Cov Power Bar plot
covPow$totals = covPow$n_pop * covPow$n_pop * covPow$sample_size
(covPow1 = ggplot(covPow,aes(x = factor(sample_size), y = power, group = bin, fill = bin)) + 
    geom_bar(position = "dodge", stat = "identity")+
    theme_classic(base_size = 20, base_family = "Times")+ 
    scale_fill_viridis(discrete = TRUE)+
    xlab("Sample Size") + ylab(expression("Power (1-"*beta*")"))+
    labs(fill = "Effect Size")+
    theme(axis.text = element_text(colour = "black"))+
    ggtitle(expression("Reciprocal Transplant: Cov"["GE"]*" Power"))+
    theme(plot.title = element_text(size = 24, face = "bold"))+
    facet_grid(~label*n_pop))

###################################
###  Confusion Matrix on Means  ###
###################################

# MEANS -- CovGE -- Permutation
dat_csv1$meansCovconfintperm = rep("NA",nrow(dat_csv1))
for(i in 1:nrow(dat_csv1)){
  if(dat_csv1$true_cov_means[i] != 0 && dat_csv1$cov_means_pvalue[i] <= 0.025){dat_csv1$meansCovconfintperm[i] = "True Positive"
  }else if(dat_csv1$true_cov_means[i] == 0 & dat_csv1$cov_means_pvalue[i] <= 0.025){dat_csv1$meansCovconfintperm[i] = "False Positive"
  }else if(dat_csv1$true_cov_means[i]!= 0 & dat_csv1$cov_means_pvalue[i] > 0.025){dat_csv1$meansCovconfintperm[i] = "False Negative"
  }else if(dat_csv1$true_cov_means[i] == 0 & dat_csv1$cov_means_pvalue[i] > 0.025){dat_csv1$meansCovconfintperm[i] = "True Negative"
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

## Counts for table

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

## Counts for heatmaps
(means_confusion_hmap1 = dat_csv1 %>%
    group_by(sample_size, n_pop, "name" =meansCovconfintperm) %>%
    summarize("n" = n()))
means_conf1 = fpr.fnr(means_confusion_hmap1, divided = TRUE, scenario = 1)
means_conf_plot1 <- heatmap_fun(means_conf1,"rate") # Can also do "percent"

(means_confusion_hmap2 = dat_csv1 %>%
    group_by(sample_size, n_pop, "name" =MeansCovconfintboot) %>%
    summarize("n" = n()))
means_conf2 = fpr.fnr(means_confusion_hmap2, divided = TRUE, scenario = 1)
means_conf_plot2 <- heatmap_fun(means_conf2,"rate")

(means_confusion_hmap3 = dat_csv1 %>%
    group_by(sample_size, n_pop,"name" = MeanGxEconfintboot) %>%
    summarize("n" = n()))
means_conf3 = fpr.fnr(means_confusion_hmap3, divided = TRUE, scenario = 1)
means_conf_plot3 <- heatmap_fun(means_conf3,"rate")

means_confusion_hmap4 = dat_csv1 %>%
  group_by(sample_size, n_pop,"name" = meansGxEconfintperm) %>%
  summarize("n" = n())
means_conf4 = fpr.fnr(raw_confusion_hmap4, divided = TRUE, scenario = 1)
means_conf_plot4 <- heatmap_fun(means_conf4,"rate")

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

############ Confusion Plots  -- Means -- Env Scenario 1 ###############

## False Positive Rates 

#(falsePos = ggplot(filter(fpdf, ID %in% c("Cov_Perm","GxE_Perm","GxE_Anova")), aes(x = ID, y = rate, group = sample_size, fill = factor(sample_size)))+ 
(means_falsePos = ggplot(fpdf_means, aes(x = ID, y = rate, group = sample_size, fill = factor(sample_size)))+ 
   geom_bar(position = "dodge", stat = "identity") + 
   geom_hline(aes(yintercept = 0.05),linetype = "dashed")+
   ylab("False Positive Rate") + xlab("")+
   facet_wrap(~npop_plot) + 
   labs(fill = "Sample Size")+
   #ylim(0,0.3)+
   scale_fill_viridis(discrete = TRUE)+
   scale_x_discrete(labels=c("Cov_Perm" = "Perm. \n CovGE", 
                             "Cov_Boot" = "Boot. \n CovGE",
                             "GxE_Boot" = "Boot. \n GxE",
                             "GxE_Perm" = "Perm. \n GxE"))+
   theme_classic2(base_family = "Times",base_size = 16) + 
   theme(axis.text = element_text(colour = "black")))



########################################################
##          Confusion Matrices  -- Env Scenario 2     ##
########################################################

# CovGE -- Permutation - Scenario 2
dat_dub$Covconfintperm = rep("NA",nrow(dat_dub))

for(i in 1:nrow(dat_dub)){
  if(dat_dub$true_cov[i] != 0 && dat_dub$covariance_pvalue[i] <= 0.025){dat_dub$Covconfintperm[i] = "True Positive"
  }else if(dat_dub$true_cov[i] == 0 & dat_dub$covariance_pvalue[i] <= 0.025){dat_dub$Covconfintperm[i] = "False Positive"
  }else if(dat_dub$true_cov[i]!= 0 & dat_dub$covariance_pvalue[i] > 0.025){dat_dub$Covconfintperm[i] = "False Negative"
  }else if(dat_dub$true_cov[i] == 0 & dat_dub$covariance_pvalue[i] > 0.025){dat_dub$Covconfintperm[i] = "True Negative"
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

cov_perm_table = dat_dub1 %>%
  group_by("name" = Covconfintperm) %>%
  summarize("n" = n())
fpr.fnr(cov_perm_table, divided = FALSE, scenario = 2)

cov_boot_table = dat_dub1 %>%
  group_by("name" =Covconfintboot) %>%
  summarize("n" = n())
fpr.fnr(cov_boot_table, divided = FALSE, scenario = 2)

gxe_anova_table = dat_dub2 %>%
  group_by("name" = GxEanova_conf) %>%
  summarize("n" = n())
fpr.fnr(gxe_anova_table, divided = FALSE, scenario = 2)

gxe_perm_table = dat_dub2 %>%
  group_by("name" =GxEconfintperm) %>%
  summarize("n" = n())
fpr.fnr(gxe_perm_table, divided = FALSE, scenario = 2)

gxe_boot_table = dat_dub2 %>%
  group_by("name" =GxEconfintboot) %>%
  summarize("n" = n())
fpr.fnr(gxe_boot_table, divided = FALSE, scenario = 2)

## Counts for heatmaps
(raw_confusion_hmap1 = dat_dub1 %>%
    group_by(sample_size, n_pop, "name" =Covconfintperm) %>%
    summarize("n" = n()))
dub_raw_conf1 = fpr.fnr(raw_confusion_hmap1, divided = TRUE, scenario = 2)
raw_conf_plot1 <- heatmap_fun(dub_raw_conf1,"rate") # Can also do "percent"

(raw_confusion_hmap2 = dat_dub1 %>%
    group_by(sample_size, n_pop, "name" =Covconfintboot) %>%
    summarize("n" = n()))
dub_raw_conf2 = fpr.fnr(raw_confusion_hmap2, divided = TRUE, scenario = 2)
raw_conf_plot2 <- heatmap_fun(dub_raw_conf2,"rate")

(raw_confusion_hmap3 = dat_dub2 %>%
    group_by(sample_size, n_pop,"name" = GxEconfintboot) %>%
    summarize("n" = n()))
dub_raw_conf3 = fpr.fnr(raw_confusion_hmap3, divided = TRUE, scenario = 2)
raw_conf_plot3 <- heatmap_fun(dub_raw_conf3,"rate")

raw_confusion_hmap4 = dat_dub2 %>%
  group_by(sample_size, n_pop,"name" = GxEconfintperm) %>%
  summarize("n" = n())
dub_raw_conf4 = fpr.fnr(raw_confusion_hmap4, divided = TRUE, scenario = 2)
raw_conf_plot4 <- heatmap_fun(dub_raw_conf4,"rate")

raw_confusion_hmap5 = dat_dub2 %>%
  group_by(sample_size, n_pop,"name" = GxEanova_conf) %>%
  summarize("n" = n())
dub_raw_conf5 = fpr.fnr(raw_confusion_hmap5, divided = TRUE, scenario = 2)
raw_conf_plot5 <- heatmap_fun(dub_raw_conf5,"rate")

# Compile FPs for plot
dub_raw_conf1$ID = rep("Cov_Perm", nrow(dub_raw_conf1))
dub_raw_conf2$ID = rep("Cov_Boot", nrow(dub_raw_conf2))
dub_raw_conf3$ID = rep("GxE_Boot", nrow(dub_raw_conf3))
dub_raw_conf4$ID = rep("GxE_Perm", nrow(dub_raw_conf4))
dub_raw_conf5$ID = rep("GxE_Anova", nrow(dub_raw_conf5))
dub_fpdf = rbind(dub_raw_conf1,dub_raw_conf2,dub_raw_conf3,dub_raw_conf4,dub_raw_conf5)
dub_fpdf = dub_fpdf[dub_fpdf$name == "False Positive",]

dub_gxe = rbind(dub_raw_conf4,dub_raw_conf5)
dub_gxe = dub_gxe[dub_gxe$name == "False Positive",]


## False Positive Rates
dub_fpdf$npop_plot = NA
for(i in 1:nrow(dub_fpdf)){
  if(dub_fpdf$n_pop[i] == 4){dub_fpdf$npop_plot[i] = "4 Genotypes"
  }else if(dub_fpdf$n_pop[i] == 8){dub_fpdf$npop_plot[i] = "8 Genotypes"
  }else{dub_fpdf$npop_plot[i] = "16 Genotypes"}
}
dub_fpdf$npop_plot = factor(dub_fpdf$npop_plot, levels=c("4 Genotypes","8 Genotypes","16 Genotypes")) 


## False Negative Rates
covperm2 = fnr.effsize(dat_dub1, metric = "cov", analysis ="perm",scenario = 2)
covperm2$label = rep("Cov_Perm",nrow(covperm2))
covboot2 = fnr.effsize(dat_dub1, metric = "cov", analysis ="boot",scenario = 2)
covboot2$label = rep("Cov_Boot",nrow(covboot2))
gxeperm2 = fnr.effsize(dat_dub1, metric = "gxe", analysis ="perm",scenario = 2)
gxeperm2$label = rep("GxE_Perm",nrow(gxeperm2))
gxeboot2 = fnr.effsize(dat_dub1, metric = "gxe", analysis = "boot",scenario = 2)
gxeboot2$label = rep("GxE_Boot",nrow(gxeboot2))
gxeanova2 = fnr.effsize(dat_dub1, metric = "gxe", analysis = "anova",scenario = 2)
gxeanova2$label = rep("GxE_Anova",nrow(gxeanova2))

fndf2 = rbind(covperm2,covboot2,gxeperm1,gxeboot1,gxeanova1)
gxePow2 = rbind(gxeperm2,gxeboot2,gxeanova2)
covPow2 = rbind(covperm2,covboot2)

############ Confusion Plots -- Env Scenario 2 ###############

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




###############################################
###  Confusion Matrix on Means - Env Scen 2 ###
###############################################

dat_dub$meansCovconfintperm = rep("NA",nrow(dat_dub))

for(i in 1:nrow(dat_dub)){
  if(dat_dub$true_cov_means[i] != 0 && dat_dub$cov_means_pvalue[i] <= 0.025){dat_dub$meansCovconfintperm[i] = "True Positive"
  }else if(dat_dub$true_cov_means[i] == 0 & dat_dub$cov_means_pvalue[i] <= 0.025){dat_dub$meansCovconfintperm[i] = "False Positive"
  }else if(dat_dub$true_cov_means[i]!= 0 & dat_dub$cov_means_pvalue[i] > 0.025){dat_dub$meansCovconfintperm[i] = "False Negative"
  }else if(dat_dub$true_cov_means[i] == 0 & dat_dub$cov_means_pvalue[i] > 0.025){dat_dub$meansCovconfintperm[i] = "True Negative"
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

means_cov_perm_table = dat_dub %>%
  group_by("name" = meansCovconfintperm) %>%
  summarize("n" = n())
fpr.fnr(means_cov_perm_table, divided = FALSE, scenario = 1)

means_cov_boot_table = dat_dub %>%
  group_by("name" =MeansCovconfintboot) %>%
  summarize("n" = n())
fpr.fnr(means_cov_boot_table, divided = FALSE, scenario = 1)

means_gxe_boot_table = dat_dub %>%
  group_by("name" =MeanGxEconfintboot) %>%
  summarize("n" = n())
fpr.fnr(means_gxe_boot_table, divided = FALSE, scenario = 1)

means_gxe_perm_table = dat_dub %>%
  group_by("name" =meansGxEconfintperm) %>%
  summarize("n" = n())
fpr.fnr(means_gxe_perm_table, divided = FALSE, scenario = 1)

## Counts for heatmaps
(means_confusion_hmap1 = dat_dub %>%
    group_by(sample_size, n_pop, "name" =meansCovconfintperm) %>%
    summarize("n" = n()))
dub_means_conf1 = fpr.fnr(means_confusion_hmap1, divided = TRUE, scenario = 1)
means_conf_plot1 <- heatmap_fun(dub_means_conf1,"rate") # Can also do "percent"

(means_confusion_hmap2 = dat_dub %>%
    group_by(sample_size, n_pop, "name" =MeansCovconfintboot) %>%
    summarize("n" = n()))
dub_means_conf2 = fpr.fnr(means_confusion_hmap2, divided = TRUE, scenario = 1)
means_conf_plot2 <- heatmap_fun(dub_means_conf2,"rate")

(means_confusion_hmap3 = dat_dub %>%
    group_by(sample_size, n_pop,"name" = MeanGxEconfintboot) %>%
    summarize("n" = n()))
dub_means_conf3 = fpr.fnr(means_confusion_hmap3, divided = TRUE, scenario = 1)
means_conf_plot3 <- heatmap_fun(dub_means_conf3,"rate")

means_confusion_hmap4 = dat_dub %>%
  group_by(sample_size, n_pop,"name" = meansGxEconfintperm) %>%
  summarize("n" = n())
dub_means_conf4 = fpr.fnr(raw_confusion_hmap4, divided = TRUE, scenario = 1)
means_conf_plot4 <- heatmap_fun(dub_means_conf4,"rate")

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

############ Confusion Plots -- Env Scenario 2 -- MEANS ###############

# False Positive Rates 
(dub_means_falsePos = ggplot(fpdf_means, aes(x = ID, y = rate, group = sample_size, fill = factor(sample_size)))+ 
   geom_bar(position = "dodge", stat = "identity") + 
   geom_hline(aes(yintercept = 0.05),linetype = "dashed")+
   ylab("False Positive Rate") + xlab("")+
   facet_wrap(~npop_plot) + 
   labs(fill = "Sample Size")+
   #ylim(0,0.3)+
   scale_fill_viridis(discrete = TRUE)+
   scale_x_discrete(labels=c("Cov_Perm" = "Perm. \n CovGE", 
                             "Cov_Boot" = "Boot. \n CovGE",
                             "GxE_Boot" = "Boot. \n GxE",
                             "GxE_Perm" = "Perm. \n GxE"))+
   theme_classic2(base_family = "Times",base_size = 16) + 
   theme(axis.text = element_text(colour = "black")))



#############################
## Power Analysis Heatmaps ##
#############################

# Covariance - FRT
cov_hm = fnr.effsize(dat_csv1, metric = "cov", analysis = "perm", resolution = "fine")
cov_hm1 = cov_hm %>%
  filter(between(bin,0.2,0.6))%>%
  group_by(sample_size,n_pop) %>%
  summarize("avgpower" = mean(power))
cov_hm1$totals = cov_hm1$n_pop *cov_hm1$n_pop * cov_hm1$sample_size

# Covariance - CG
cov_hm2 = fnr.effsize(dat_dub1, metric = "cov", analysis = "perm", resolution = "fine")
cov_hm3 = cov_hm2 %>%
  filter(between(bin,0.2,0.6))%>%
  group_by(sample_size,n_pop) %>%
  summarize("avgpower" = mean(power))
cov_hm3$totals = cov_hm3$n_pop * 2 * cov_hm3$sample_size

# HeatMaps - Covariance based on Permutation - FRT
(CovPower1 = ggplot(cov_hm1,aes(x = factor(sample_size), y = factor(n_pop), fill = avgpower)) + 
    geom_tile() + 
    geom_text(aes(label= paste(cov_hm1$totals, round(cov_hm1$avgpower,3),sep = '\n')), size = 5) +
    theme_classic(base_size = 24, base_family = "Times")+ 
    scale_fill_viridis(
      breaks=seq(0,1,0.25), #breaks in the scale bar
      limits=c(0,1))+
    xlab("Sample Size") + ylab("Number of Genotypes")+
    labs(fill = "Power")+
    theme(axis.text = element_text(colour = "black"))+  
    theme(legend.position = "none")+
    ggtitle(expression("Cov"["GE"]*": Full Reciprocal Transplant"))+
    theme(plot.title = element_text(size = 24, face = "bold")))

# HeatMaps - Covariance based on Permutation - FRT
(CovPower2 = ggplot(cov_hm3,aes(x = factor(sample_size), y = factor(n_pop), fill = avgpower)) + 
    geom_tile() + 
    geom_text(aes(label= paste(cov_hm3$totals, round(cov_hm3$avgpower,3),sep = '\n')), size = 5) +
    theme_classic(base_size = 24, base_family = "Times")+ 
    scale_fill_viridis(
      breaks=seq(0,1,0.25), #breaks in the scale bar
      limits=c(0,1))+
    xlab("Sample Size") + ylab("Number of Genotypes")+
    labs(fill = "Power")+
    theme(axis.text = element_text(colour = "black"))+  
    theme(legend.position = "none")+
    ggtitle(expression("Cov"["GE"]*": Paired Common Garden"))+
    theme(plot.title = element_text(size = 24, face = "bold")))

# Covariance - FRT
gxe_hm = fnr.effsize(dat_csv1, metric = "gxe", analysis = "perm", resolution = "fine")
gxe_hm1 = gxe_hm %>%
  filter(between(bin,0.3,0.6))%>%
  group_by(sample_size,n_pop) %>%
  summarize("avgpower" = mean(power))
gxe_hm1$totals = gxe_hm1$n_pop *gxe_hm1$n_pop * gxe_hm1$sample_size

# Covariance - CG
gxe_hm2 = fnr.effsize(dat_dub1, metric = "gxe", analysis = "perm", resolution = "fine")
gxe_hm3 = gxe_hm2 %>%
  filter(between(bin,0.3,0.6))%>%
  group_by(sample_size,n_pop) %>%
  summarize("avgpower" = mean(power))
gxe_hm3$totals = gxe_hm3$n_pop * 2 * gxe_hm3$sample_size

# HeatMaps - Covariance based on Permutation - FRT
(GxEPower1 = ggplot(gxe_hm1,aes(x = factor(sample_size), y = factor(n_pop), fill = avgpower)) + 
    geom_tile() + 
    geom_text(aes(label= paste(gxe_hm1$totals, round(gxe_hm1$avgpower,3),sep = '\n')), size = 5) +
    theme_classic(base_size = 24, base_family = "Times")+ 
    scale_fill_viridis(
      breaks=seq(0,1,0.25), #breaks in the scale bar
      limits=c(0,1))+
    xlab("Sample Size") + ylab("Number of Genotypes")+
    labs(fill = "Power")+
    theme(axis.text = element_text(colour = "black"))+  
    theme(legend.position = "none")+
    ggtitle(expression(bar(Delta)*""["GxE"]*": Full Reciprocal Transplant"))+
    theme(plot.title = element_text(size = 24, face = "bold")))

# HeatMaps - Covariance based on Permutation - FRT
(GxEPower2 = ggplot(gxe_hm3,aes(x = factor(sample_size), y = factor(n_pop), fill = avgpower)) + 
    geom_tile() + 
    geom_text(aes(label= paste(gxe_hm3$totals, round(gxe_hm3$avgpower,3),sep = '\n')), size = 5) +
    theme_classic(base_size = 24, base_family = "Times")+ 
    scale_fill_viridis(
      breaks=seq(0,1,0.25), #breaks in the scale bar
      limits=c(0,1))+
    xlab("Sample Size") + ylab("Number of Genotypes")+
    labs(fill = "Power")+
    theme(axis.text = element_text(colour = "black"))+  
    #theme(legend.position = "none")+
    ggtitle(expression(bar(Delta)*""["GxE"]*": Paired Common Garden"))+
    theme(plot.title = element_text(size = 24, face = "bold")))

grid.arrange(CovPower1,GxEPower1,CovPower2,GxEPower2,ncol = 2)


######################################
## Tradeoff with GxE and Covariance ##
######################################

dat_csv$sig = NULL
for(i in 1:nrow(dat_csv)){ # Use only if one or the other is significant
  if(dat_csv$GxE_emm_pvalue[i] <= 0.05 | dat_csv$covariance_lwrCI[i] > 0 & dat_csv$covariance_uprCI[i] > 0){dat_csv$sig[i] = TRUE
  }else if(dat_csv$GxE_emm_pvalue[i] <= 0.05 | dat_csv$covariance_lwrCI[i] < 0 & dat_csv$covariance_uprCI[i] < 0){dat_csv$sig[i] = TRUE
  }else{dat_csv$sig[i]=FALSE} 
}

dat_dub$sig = NULL
for(i in 1:nrow(dat_dub)){ # Use only if one or the other is significant
  if(dat_dub$GxE_emm_pvalue[i] <= 0.05 | dat_dub$covariance_lwrCI[i] > 0 & dat_dub$covariance_uprCI[i] > 0){dat_dub$sig[i] = TRUE
  }else if(dat_dub$GxE_emm_pvalue[i] <= 0.05 | dat_dub$covariance_lwrCI[i] < 0 & dat_dub$covariance_uprCI[i] < 0){dat_dub$sig[i] = TRUE
  }else{dat_dub$sig[i]=FALSE} 
}

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

sigGxE = dat_csv %>% 
  filter(sig ==TRUE) %>% # filter out false positives as potential solution to weed out messiness.
  filter(true_GxE_emm != 0) %>%
  filter(true_cov != 0) %>%
  filter(Covconfintboot != "false positive") %>%
  filter(GxEconfintperm != "false positive")

sigGxE2 = dat_dub %>% 
  filter(sig ==TRUE) %>% # filter out false positives as potential solution to weed out messiness.
  filter(true_GxE_emm != 0) %>%
  filter(true_cov != 0) %>%
  filter(Covconfintboot != "false positive") %>%
  filter(GxEconfintperm != "false positive")

(bin = ggplot(sigGxE, aes(x = true_GxE_emm, y = covtick))+
    geom_smooth(method = "glm",method.args = list(family = "binomial"),se = T,colour = "black") + 
    geom_point()+
    xlab(expression("Actual "*bar(Delta)*""["GxE"]))+ylab(expression("Proportion significant Cov"["GE"]))+
    theme_bw(base_size = 18, base_family = "Times")+
    theme(axis.text.x = element_text(colour = "black"))+
    theme(axis.text.y = element_text(colour = "black"))+
    ggtitle("")+
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 2)))


(lin = ggplot(filter(sigGxE,replicate ==1), aes(x = true_GxE_emm, y = abs(true_cov)))+
    geom_point(alpha = 0.15)+
    geom_smooth(method = "glm", colour = "black", size = 1.5)+
    xlab(expression("Actual "*bar(Delta)*""["GxE"]))+ylab(expression("| Actual Cov"["GE"]*" |"))+
    theme_bw(base_size = 18, base_family = "Times")+
    theme(axis.text.x = element_text(colour = "black"))+
    theme(axis.text.y = element_text(colour = "black"))+
    ggtitle("Full Reciprocal Transplant Design")+
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 2)))

(bin2 = ggplot(sigGxE2, aes(x = true_GxE_emm, y = covtick))+
    geom_smooth(method = "glm",method.args = list(family = "binomial"),se = T,colour = "black") + 
    geom_point()+
    xlab(expression("Actual "*bar(Delta)*""["GxE"]))+ylab(expression("Proportion significant Cov"["GE"]))+
    theme_bw(base_size = 18, base_family = "Times")+
    theme(axis.text.x = element_text(colour = "black"))+
    theme(axis.text.y = element_text(colour = "black"))+
    ggtitle("")+
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 2)))


(lin2 = ggplot(filter(sigGxE2, replicate == 1), aes(x = true_GxE_emm, y = abs(true_cov)))+
    geom_point(alpha = 0.15)+
    geom_smooth(method = "glm", colour = "black", size = 1.5)+
    xlab(expression("Actual "*bar(Delta)*""["GxE"]))+ylab(expression("| Actual Cov"["GE"]*" |"))+
    theme_bw(base_size = 18, base_family = "Times")+
    theme(axis.text.x = element_text(colour = "black"))+
    theme(axis.text.y = element_text(colour = "black"))+
    ggtitle("Paired Common Garden Design")+
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 2)))

grid.arrange(lin, bin, lin2, bin2, ncol = 2)

### Same but with Omega^2
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
oe = dat_csv %>% filter(replicate == 1) %>% filter(std_dev == 1)
                                        
(eff_var1 = ggplot(oe, aes(x = true_GxE_emm, y = GxE_omega)) + 
    geom_point(aes(shape = factor(sigsig), fill = factor(sigsig)), alpha = 0.65, size = 4)+
    scale_fill_manual(values = devcol)+
                      #breaks = c("0.5", "1"))+
    scale_shape_manual(values = devshape)+
    geom_smooth(method = "lm", formula = y ~ splines::bs(x, 3),se = F,colour = "black",size = 1.5) + 
    xlab(expression(""*bar(Delta)*""["GxE"]*" of Population"))+
    #scale_linetype_manual(values = c("0.5" = "dotdash", "1" = "solid"))+
    ylab(expression(""*omega^2*" of Population"))+
    #guides(fill=guide_legend(override.aes=list(shape=21,size =6)))+
    labs(fill = "Significance", shape =  "Significance", linetype = "Standard Deviation")+
    guides(shape = guide_legend(override.aes = list(size = 6)))+
    #ggtitle("Full Reciprocal Transplant Design")+
    theme_classic(base_size = 24, base_family = "Times")+
    theme(axis.text =element_text(colour = "black")))


dat_dub$sigsig = NA
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
(covpopcheck = ggplot(dat_csv_2,aes(x = abs(true_cov), y = abs(covariance)))+
    geom_point(aes(colour = factor(n_pop)))+ylab(expression("Cov"["GE"]*": Group Means"))+xlab(expression("Cov"["GE"]*": Raw data"))+
    geom_abline(slope = 1, intercept = 0,colour = "red")+
    theme_classic(base_size = 20, base_family = "Times")+
    ggtitle("Full reciprocal transplant design")+
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

grid.arrange(coverrorcheck_lwr, coverrorcheck_upr)

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

grid.arrange(gxeerrorcheck, gxeerrorcheck_lwr, gxeerrorcheck_upr, ncol = 1)

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
(rowpicker = dat_csv %>%
    filter(sample_size == 4)%>%
    filter(n_pop == 4)%>%
    filter(GxE_emm < 0.35)%>% #No GxE
    #filter(GxE_emm > 0.65)%>% # GxE
    #filter(covariance < -.35)) # CnGV
    filter(covariance > .75)) # CoGV

chosen_1 = c(45657, # CoGV, GxE, 2 pop #
             11021, # CnGV, GxE, 2 pop #
             44014, # CnGV, No GxE, 2 pop #
             16434, # CoGV, No GxE, 2 pop #
             
            # 12949, # CoGV, No GxE, 8 pop #
            # 16546, # CnGV, No GxE, 8 pop #
            # 13709, # CnGV, GxE, 8 pop #
            # 17324, # CoGV, GxE, 8 pops #
             
             9279, # CoGV, No GxE, 4 pop #
             299, # CnGV, No GxE, 4 pop #
             6484, # CnGV, GxE, 4 pop #
             13662) # CoGV, GxE, 4 pops #

chosen_2 = c(21997, # CoGV, GxE, 2 pop per env #
             30830, # CnGV, GxE, 2 pop per env #
             54780, # CnGV, No GxE, 2 pop per env #
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

for(i in 1:length(chosen)){

phenRow = filter(start_df,row == chosen[i])
plotdat = filter(phen_data, row == chosen[i])
label1 = paste0(phenRow$covariance,"; P = ",phenRow$covariance_pvalue)
label2 = paste0(phenRow$GxE_emm,"; P = ",phenRow$GxE_emm_pvalue)
colorpal = c("E_1" = "#404788FF", "E_2" = "#73D055FF")
shape_2 = c("E_1" = 15, "E_2" =  17)
shape_4 = c("E_1" = 15, "E_2" = 17,"E_3" = 19, "E_4"= 18)
shape_8 = c("E_1" = 15, "E_2" = 17,"E_3" = 19, "E_4"= 18,
            "E_5" = 15, "E_6" = 17,"E_7" = 19, "E_8"= 18)


p = ggplot(plotdat,aes(x = exp_env_factor, y = phen_corrected, group = gen_factor,colour = nat_env_factor))+
  geom_point(aes(shape = nat_env_factor),size = 5, position=position_dodge(width = 0.15))+
  geom_smooth(size = 4, se=FALSE)+
  theme_classic(base_size = 40, base_family = "Times")+
  ylim(-3,6)+
  labs(colour = "",shape = "")+
  guides(colour = guide_legend(ncol = 2))+
  ylab("Phenotype")+xlab("Environment")+
  annotate("text", x = 0.5, y = 5.5, label = deparse(bquote("Cov"["GE"] ==~.(label1))),size=16, color = "black",hjust = 0,parse = T)+
  annotate("text", x = 0.5, y = 4, label = deparse(bquote(bar(Delta)*""["GxE"] ==~.(label2))),size=16 , color = "black", hjust = 0,parse = T)+
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
pdf(paste("Row_", phenRow$row, ".pdf", sep = ""), width=11, height=8.5) # start export
print(p3) 
dev.off() # finish export
}
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




#######################################
#########     Extra Code      #########
#######################################

## Old Power Analysis approach:

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



