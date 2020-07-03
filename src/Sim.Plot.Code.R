##########################################################
###             Simulation Visualization Code          ###
##########################################################

## Load packages

library(ggplot2)
library(readr)
library(tidyverse)
library(gridExtra)

## Compile Data

setwd("~/Desktop/power_output/")
Param_temp <- list.files(path = "~/Desktop/power_output/",pattern= "*Parameter") 
Cov_temp <-list.files(path = "~/Desktop/power_output/",pattern= "*Covariance") 
GxE_temp <- list.files(path = "~/Desktop/power_output/",pattern= "*GxE") 
GE = plyr::ldply(GxE_temp, read_csv)
Cov = plyr::ldply(Cov_temp, read_csv)
Par = plyr::ldply(Param_temp, read_csv)

dat_csv. = full_join(Par, GE, by= "row")
dat_csv = full_join(dat_csv., Cov, by = "row")

# Once code is read in, write csv to skip lengthy loading
# write.csv(dat_csv,"~/Desktop/dat_csv.csv")
dat_csv <- read.csv("~/Desktop/dat_csv.csv")

######################
## Covariance x GxE ##
######################

# Assign colors for Cov x GxE plot

dat_csv$col = NULL
for(i in 1:nrow(dat_csv)){
  if(dat_csv$covariance_pvalue[i] <= 0.05 & dat_csv$GxE_emm_pvalue[i] <= 0.05){dat_csv$col[i] = "red" # Both significant
  }else if(dat_csv$covariance_pvalue[i] <= 0.05 & dat_csv$GxE_emm_pvalue[i] > 0.05){dat_csv$col[i] = "darkgreen" # Cov significant
  }else if(dat_csv$covariance_pvalue[i] > 0.05 & dat_csv$GxE_emm_pvalue[i] <= 0.05){dat_csv$col[i] = "dodgerblue4" # GxE significant
  }else{dat_csv$col[i] = "grey"} # None significant
}

# Cov x GxE Plot

ggplot(dat_csv, aes(x = covariance, y = GxE_emm, group = factor(n_pop), alpha = 0.1,colour = col)) + 
  geom_jitter() + theme_classic() + ylim(0,1) + xlim(-1,1)+
  xlab("Covariance Estimate") + ylab("GxE Estimate") +
  theme(legend.position = "none")+
  scale_colour_identity()+
  facet_grid(sample_size~n_pop)

# Cov x Omega2 Plot

ggplot(dat_csv, aes(x = covariance, y = GxE_omega, group = factor(n_pop), alpha = 0.1,colour = col)) + 
  geom_jitter() + theme_classic() + 
  xlab("Covariance Estimate") + ylab("GxE Estimate (Omega^2)") +
  theme(legend.position = "none")+
  scale_colour_identity()+
  facet_grid(sample_size~n_pop)


#############################
## Power Analysis Heatmaps ##
#############################

# Assign 1's and 0's for significant outcomes

dat_csv$covtick <- NULL
for(i in 1:nrow(dat_csv)){
  if(dat_csv$covariance_pvalue[i] > 0.05){dat_csv$covtick[i]=0}else{dat_csv$covtick[i]=1} 
}

dat_csv$gxetick <- NULL
for(i in 1:nrow(dat_csv)){
  if(dat_csv$GxE_emm_pvalue[i] > 0.05){dat_csv$gxetick[i]=0}else{dat_csv$gxetick[i]=1}
}

# Power at moderate levels of Cov and GxE

covpow1 = dat_csv %>%
  filter(between(covariance, 0.4,0.6))
gxepow1 = dat_csv %>%
  filter(between(GxE_emm,0.4,0.6))

covpow = covpow1 %>%
  group_by(delta_env, delta_gen, sample_size, n_pop, std_dev, interaction) %>%
  summarize("total_tick" = n(),
            "covtick" = sum(covtick))
covpow$covpower = covpow$covtick/covpow$total_tick

gxepow = gxepow1 %>%
  group_by(delta_env, delta_gen, sample_size,n_pop,std_dev,interaction) %>%
  summarize("total_tick" = n(),
            "gxetick" = sum(gxetick))
gxepow$gxepower = gxepow$gxetick/gxepow$total_tick

# Divide datasets according to standard deviation, take the mean power across different parameter sets

covpow_lowsd = covpow %>%
  filter(std_dev == min(covpow$std_dev)) %>%
  group_by(sample_size, n_pop) %>%
  summarize("meancovpower" = mean(covpower))
covpow_lowsd$pop_samp = paste(covpow_lowsd$n_pop,covpow_lowsd$sample_size)

covpow_highsd = covpow %>%
  filter(std_dev == max(covpow$std_dev)) %>%
  group_by(sample_size, n_pop) %>%
  summarize("meancovpower" = mean(covpower))
covpow_highsd$pop_samp = paste(covpow_highsd$n_pop,covpow_highsd$sample_size)

gxepow_lowsd = gxepow %>%
  filter(std_dev == min(gxepow$std_dev)) %>%
  group_by(sample_size, n_pop) %>%
  summarize("meangxepower" = mean(gxepower))
gxepow_lowsd$pop_samp = paste(gxepow_lowsd$n_pop,gxepow_lowsd$sample_size)

gxepow_highsd = gxepow %>%
  filter(std_dev == max(gxepow$std_dev)) %>%
  group_by(sample_size, n_pop) %>%
  summarize("meangxepower" = mean(gxepower))
gxepow_highsd$pop_samp = paste(gxepow_highsd$n_pop,gxepow_highsd$sample_size)

# Plot Moderate levels of each in grid

(csdlow = ggplot(covpow_lowsd,aes(x = factor(sample_size), y = factor(n_pop), fill = meancovpower)) + 
  geom_tile() + scale_fill_gradient(low="grey", high="darkgreen") +
  xlab("Sample Size") + ylab("Number of Populations")+
  labs(fill = "Power")+
  theme_classic() + ggtitle(paste0("Covariance: Standard Deviation = ",min(covpow$std_dev))))

(csdhigh = ggplot(covpow_highsd,aes(x = factor(sample_size), y = factor(n_pop), fill = meancovpower)) + 
  geom_tile() + scale_fill_gradient(low="grey", high="darkgreen") +
  xlab("Sample Size") + ylab("Number of Populations")+
  labs(fill = "Power")+
  theme_classic()+ ggtitle(paste0("Covariance: Standard Deviation = ",max(covpow$std_dev))))

(gsdlow = ggplot(gxepow_lowsd,aes(x = factor(sample_size), y = factor(n_pop), fill = meangxepower)) + 
  geom_tile() + scale_fill_gradient(low="grey", high="dodgerblue4") +
  xlab("Sample Size") + ylab("Number of Populations")+
  labs(fill = "Power")+
  theme_classic() + ggtitle(paste0("GxE: Standard Deviation = ",min(gxepow$std_dev))))

(gsdhigh = ggplot(gxepow_highsd,aes(x = factor(sample_size), y = factor(n_pop), fill = meangxepower)) + 
  geom_tile() + scale_fill_gradient(low="grey", high="dodgerblue4") +
  xlab("Sample Size") + ylab("Number of Populations")+
  labs(fill = "Power") +
  theme_classic() + ggtitle(paste0("GxE: Standard Deviation = ",max(gxepow$std_dev))))

grid.arrange(gsdhigh,gsdlow,csdhigh,csdlow,ncol = 2)

# Line-style Power plots

dat_csv$mod = NULL
for(i in 1:nrow(dat_csv)){
  if(dat_csv$covariance[i] >= 0.4 & dat_csv$covariance[i] <= 0.6){dat_csv$mod[i] = TRUE
  }else if(dat_csv$GxE_emm[i] >= 0.4 & dat_csv$GxE_emm[i] <= 0.6){dat_csv$mod[i] = TRUE
  }else{dat_csv$mod[i] = FALSE}
}

col_df = dat_csv %>%
  filter(mod == TRUE) %>%
  group_by(n_pop,sample_size,col) %>%
  summarize(frequency = n()) 

new_df = data.frame()
for(i in unique(col_df$sample_size)){
  for(j in unique(col_df$n_pop)){
    subset = col_df[col_df$sample_size==i,]
    subsetsubset = subset[subset$n_pop == j,]
    colproportion = subsetsubset$frequency/sum(subsetsubset$frequency)
    new_df. = data.frame(subsetsubset,"proportion"=colproportion)
    new_df = rbind(new_df,new_df.)
  }
}

(linepower =ggplot(new_df,aes(x = factor(n_pop), y = proportion, group = sample_size, fill = factor(col,levels=c("grey","dodgerblue4","darkgreen","red"))))+
  geom_bar(position = "stack", stat="identity") +
    scale_fill_identity() +
    facet_wrap(~sample_size) +theme_classic())

red_dat = filter(new_df,col=="red")
ggplot(red_dat,aes(x = factor(n_pop), y = proportion,group = factor(sample_size), colour = factor(sample_size)))+
  geom_point(size = 4)+geom_line()+#theme_classic()+
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  theme_bw(base_size = 20, base_family = "Times")+
  theme(axis.text.x = element_text(size=16,colour = "black"),
        axis.title.x = element_text(size=18,face="bold")) +
  theme(axis.text.y = element_text(size=16,colour = "black"),
        axis.title.y = element_text(size=18,face="bold")) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2))+
  xlab("Number of Populations")+ylab("Proportion of significant Cov and GxE")+
  labs(colour = "Sample Size")

blue_dat = filter(new_df,col=="dodgerblue4")
ggplot(blue_dat,aes(x = factor(n_pop), y = proportion,group = factor(sample_size), colour = factor(sample_size)))+
  geom_point(size = 4)+geom_line()+#theme_classic()+
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  theme_bw(base_size =20, base_family = "Times")+
  theme(axis.text.x = element_text(size=16,colour = "black"),
        axis.title.x = element_text(size=18,face="bold")) +
  theme(axis.text.y = element_text(size=16,colour = "black"),
        axis.title.y = element_text(size=18,face="bold")) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2))+
  xlab("Number of Populations")+ylab("Proportion of significant GxE")+
  labs(colour = "Sample Size")

green_dat = filter(new_df,col=="darkgreen")
ggplot(green_dat,aes(x = factor(n_pop), y = proportion,group = factor(sample_size), colour = factor(sample_size)))+
  geom_point(size = 4)+geom_line()+#theme_classic()+
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  theme_bw(base_size = 20, base_family = "Times")+
  theme(axis.text.x = element_text(size=16,colour = "black"),
        axis.title.x = element_text(size=18,face="bold")) +
  theme(axis.text.y = element_text(size=16,colour = "black"),
        axis.title.y = element_text(size=18,face="bold")) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2))+
  xlab("Number of Populations")+ylab("Proportion of significant CovGE)")+
  labs(colour = "Sample Size")

# Tradeoff between GxE and Covariance plot

sigGxE = filter(dat_csv, GxE_emm_pvalue <=0.05 | covariance_pvalue <= 0.05)
ggplot(sigGxE, aes(x = GxE_emm, y = covtick))+
  geom_smooth(method = "glm",method.args = list(family = "binomial"),se = T,colour = "black") + 
  xlab("Magnitude of GxE")+ylab("Proportion of significant CovGE values (p < 0.05)")+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(axis.text.x = element_text(size=16,colour = "black"),
        axis.title.x = element_text(size=18,face="bold")) +
  theme(axis.text.y = element_text(size=16,colour = "black"),
        axis.title.y = element_text(size=18,face="bold")) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2))

######################################
##          Sanity Checks           ##
######################################

# Do we have good parameter coverage? Hex plot

(hexy = ggplot(dat_csv, aes(x = covariance, y = GxE_emm)) + 
  geom_hex()+ theme_classic() + xlab("Covariance Estimate") + ylab("GxE Estimate") +
  ggtitle("HexPlot") + facet_grid(sample_size~n_pop))

# Do confidence intervals match with p-values? 

dat_csv$Covconfint = NULL
dat_csv$GxEconfint = NULL

for(i in 1:nrow(dat_csv)){
  
  if(dat_csv$covariance_pvalue[i] > 0.025 & dat_csv$covariance_lwrCI[i] < 0 & dat_csv$covariance_uprCI[i] > 0){dat_csv$Covconfint[i] = "Result: Not Significant - Correct"
  }else if(dat_csv$covariance_pvalue[i] > 0.025 & dat_csv$covariance_lwrCI[i] > 0 & dat_csv$covariance_uprCI[i] > 0){dat_csv$Covconfint[i] = "Result: Type II Error"
  }else if(dat_csv$covariance_pvalue[i] > 0.025 & dat_csv$covariance_lwrCI[i] < 0 & dat_csv$covariance_uprCI[i] < 0){dat_csv$Covconfint[i] = "Result: Type II Error"
  }else if(dat_csv$covariance_pvalue[i] <= 0.025 & dat_csv$covariance_lwrCI[i] > 0 & dat_csv$covariance_uprCI[i] > 0){dat_csv$Covconfint[i] = "Result: Significant - Correct"
  }else if(dat_csv$covariance_pvalue[i] <= 0.025 & dat_csv$covariance_lwrCI[i] < 0 & dat_csv$covariance_uprCI[i] < 0){dat_csv$Covconfint[i] = "Result: Significant - Correct"
  }else{dat_csv$Covconfint[i] = "Result: Type I Error"}
  
  if(dat_csv$GxE_emm_pvalue[i] > 0.05 & dat_csv$GxE_emm_lwrCI[i] <= 0){dat_csv$GxEconfint[i] = "Result: Not Significant - Correct"
  }else if(dat_csv$GxE_emm_pvalue[i] > 0.05 & dat_csv$GxE_emm_lwrCI[i] > 0){dat_csv$GxEconfint[i] = "Result: Type II Error"
  }else if(dat_csv$GxE_emm_pvalue[i] <= 0.05 & dat_csv$GxE_emm_lwrCI[i] > 0){dat_csv$GxEconfint[i] = "Result: Significant - Correct"
  }else{dat_csv$GxEconfint[i] = "Result: Type I Error"}
  
}

(cov_ci = ggplot(dat_csv[dat_csv$replicate==1,], aes(x = reorder(row,covariance))) +
  geom_point(aes(y = covariance), colour = "firebrick",alpha = 1)+
  geom_point(aes(y = true_cov), colour = "springgreen4",alpha = 1)+
  geom_errorbar(aes(ymin = covariance_lwrCI, ymax = covariance_uprCI))+
  theme_classic() + geom_hline(aes(yintercept = 0))+facet_wrap(~Covconfint,ncol = 2))

(gxe_ci = ggplot(dat_csv[dat_csv$replicate==1,], aes(x = reorder(row,GxE_emm))) +
    geom_point(aes(y = GxE_emm), colour = "firebrick",alpha = 1)+
    geom_point(aes(y = true_GxE_emm), colour = "springgreen4",alpha = 1)+
    geom_errorbar(aes(ymin = GxE_emm_lwrCI, ymax = GxE_emm_uprCI))+
    theme_classic() + geom_hline(aes(yintercept = 0))+facet_wrap(~GxEconfint,ncol = 2))


# Do confidence intervals for each overlap with their true values? Nope. 
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


# Do means estimates match raw estimates? If yes should fall along 1:1 line

dat_csv$meancoverror = abs(dat_csv$cov_means_correct_uprCI - dat_csv$cov_means_correct_lwrCI)
dat_csv$coverror = abs(dat_csv$covariance_uprCI - dat_csv$covariance_lwrCI)
dat_csv$meangxeerror = dat_csv$GxE_means_uprCI - dat_csv$GxE_means_lwrCI
dat_csv$gxeerror = dat_csv$GxE_emm_uprCI - dat_csv$GxE_emm_lwrCI

# Covariance
(covmeancheck = ggplot(dat_csv,aes(x = covariance, y = cov_means_correct))+
    geom_point()+theme_classic()+ylab("Covariance from Means")+xlab("Covariance from Raw")+
    geom_abline(slope = 1, intercept = 0,colour = "red"))

(covcheck = ggplot(dat_csv,aes(x = true_cov, y = covariance, colour = col))+ # Cov = 1 happens when delta_env = delta_gen, they are both positive, and there is no interaction
    geom_point()+theme_classic()+ylab("Covariance")+xlab("True Covariance")+scale_color_identity()+ # In the no error scenarios, this makes perfectly parallel lines that always have true_cov = 1
    geom_abline(slope = 1, intercept = 0,colour = "red"))

(coverrorcheck = ggplot(dat_csv, aes(x = coverror,y = meancoverror)) + 
    geom_point(alpha = 0.5) + theme_classic() + ylab("Length of 95% CI for CovGE means") + xlab("Length of 95% CI for CovGE raw")+
    geom_abline(slope = 1, intercept = 0,colour = "red")+scale_colour_identity())

suspect = dat_csv %>% filter(true_cov == 1) %>% filter(covariance_pvalue <0.05) # Check on those weirdos

# GxE
(gxemeancheck = ggplot(dat_csv,aes(x = true_GxE_emm, y = true_GxE_means))+
    geom_point()+theme_classic()+ylab("GxE from Means")+xlab("GxE from Raw")+
    geom_abline(slope = 1, intercept = 0,colour = "red"))

(gxecheck = ggplot(dat_csv,aes(x = true_GxE_emm, y = GxE_emm,colour = col))+
    geom_point()+theme_classic()+ylab("GxE with error")+xlab("True GxE")+scale_color_identity()+
    geom_abline(slope = 1, intercept = 0,colour = "red"))

(gxeerrorcheck = ggplot(dat_csv, aes(x = gxeerror,y = meangxeerror)) + 
    geom_point(alpha = 0.3) + theme_classic() + ylab("Length of 95% CI for GxE means") + xlab("Length of 95% CI for GxE raw")+
    geom_abline(slope = 1, intercept = 0,colour = "red")+scale_colour_identity())


# Do GxE pvalues from permutation match pvalues from anova? 

(ggplot(dat_csv,aes(x = GxE_Anova, y = GxE_emm_pvalue,colour = GxE_emm))+ #raw data 
    geom_jitter()+theme_classic()+ylab("GxE EMM Pvalue")+xlab("GxE Anova Pvalue")+
    geom_vline(xintercept = 0.05,colour = "red")+
    geom_hline(yintercept = 0.05,colour = "red"))

(ggplot(dat_csv,aes(x = GxE_Anova, y = GxE_means_pvalue))+ # means data
    geom_point()+theme_classic()+ylab("GxE EMM Pvalue")+xlab("GxE Anova Pvalue")+
    geom_vline(xintercept = 0.05,colour = "red")+
    geom_hline(yintercept = 0.05,colour = "red"))

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


