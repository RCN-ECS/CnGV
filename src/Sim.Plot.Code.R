##########################################################
###             Simulation Visualization Code          ###
##########################################################

## Load packages

library(ggplot2)
library(readr)
library(tidyverse)
library(gridExtra)
library(ggthemes)

## Compile Data

PL <- list.files(path = "~/Documents/GitHub/CnGV/results/SimResults_1-10/PL_output/",pattern= "*PL_") 
Sim <- list.files(path = "~/Documents/GitHub/CnGV/results/SimResults_1-10/power_output/",pattern= "*Results_") 

setwd("~/Documents/GitHub/CnGV/results/SimResults_1-10/PL_output/")
PL1 = plyr::ldply(PL, read_csv)

setwd("~/Documents/GitHub/CnGV/results/SimResults_1-10/power_output/")
dat_csv = plyr::ldply(Sim, read_csv) 

# Once code is read in, write csv to skip lengthy loading
write.csv(dat_csv,"~/Documents/GitHub/CnGV/results/SimResults_1-10/dat_csv_1to10.csv")
write.csv(PL1, "~/Documents/GitHub/CnGV/results/SimResults_1-10/PLreps_1to10")
#dat_csv <- read.csv("~/Desktop/dat_csv.csv")

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

# Cov Perm check
dat_csv1 = dat_csv[dat_csv$replicate==1,] # Only need to show one replicate
dat_csv1$Covconfintperm = rep("NA",nrow(dat_csv1))

for(i in 1:nrow(dat_csv1)){
    if(dat_csv1$true_cov[i] != 0 && dat_csv1$covariance_pvalue[i] <= 0.025){dat_csv1$Covconfintperm[i] = "True Positive"
    }else if(dat_csv1$true_cov[i] == 0 & dat_csv1$covariance_pvalue[i] <= 0.025){dat_csv1$Covconfintperm[i] = "False Positive"
    }else if(dat_csv1$true_cov[i]!= 0 & dat_csv1$covariance_pvalue[i] > 0.025){dat_csv1$Covconfintperm[i] = "False Negative"
    }else if(dat_csv1$true_cov[i] == 0 & dat_csv1$covariance_pvalue[i] > 0.025){dat_csv1$Covconfintperm[i] = "True Negative"
    }else{dat_csv1$Covconfintperm[i] = "None"}
}

(cov_perm = ggplot(transform(dat_csv1, Covconfintperm = factor(Covconfintperm, levels = c("True Positive", "True Negative", "False Positive", "False Negative")))) +
    geom_point(aes(x = reorder(row,true_cov), y = covariance), colour = "firebrick",alpha = 1)+
    geom_point(aes(x = reorder(row,true_cov), y = true_cov), colour = "springgreen4",alpha = 1)+
    geom_errorbar(aes(x = reorder(row,true_cov), ymin = covariance_lwrCI, ymax = covariance_uprCI))+ ggtitle("Covariance - Permutation Check")+
    theme_classic() + geom_hline(aes(yintercept = 0))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+ 
    facet_wrap(~Covconfintperm,ncol = 2))


# Cov Boot check
dat_csv1$Covconfintboot = rep(NA, nrow(dat_csv1))

for(i in 1:nrow(dat_csv1)){
  
  if(dat_csv1$true_cov[i] != 0 &&
     dat_csv1$covariance_lwrCI[i] < 0 &&
    dat_csv1$covariance_uprCI[i] < 0
    ){dat_csv1$Covconfintboot[i] = "True Positive"
    }else if(dat_csv1$true_cov[i] != 0 &&
             dat_csv1$covariance_lwrCI[i] > 0 &&
    dat_csv1$covariance_uprCI[i] > 0
    ){dat_csv1$Covconfintboot[i] = "True Positive"
    }else if(dat_csv1$true_cov[i] == 0 &&
             dat_csv1$covariance_lwrCI[i] < 0 &&
    dat_csv1$covariance_uprCI[i] < 0
    ){dat_csv1$Covconfintboot[i] = "False Positive"
    }else if(dat_csv1$true_cov[i] == 0 &&
             dat_csv1$covariance_lwrCI[i] > 0 &&
    dat_csv1$covariance_uprCI[i] > 0
    ){dat_csv1$Covconfintboot[i] = "False Positive"
    }else if(dat_csv1$true_cov[i] != 0 && 
             dat_csv1$covariance_lwrCI[i] <= 0 && 
    dat_csv1$covariance_uprCI[i] >= 0
    ){dat_csv1$Covconfintboot[i] = "False Negative"
    }else if(dat_csv1$true_cov[i]== 0 && 
             dat_csv1$covariance_lwrCI[i] <= 0 && 
    dat_csv1$covariance_uprCI[i] >= 0
   ){dat_csv1$Covconfintboot[i] = "True Negative"
   }else{dat_csv1$Covconfintboot[i] = "None"}

  }

(cov_boot = ggplot(transform(dat_csv1, Covconfintboot = factor(Covconfintboot, levels = c("True Positive", "True Negative", "False Positive", "False Negative")))) +
    geom_point(aes(x = reorder(row,true_cov),y = covariance), colour = "firebrick",alpha = 1)+
    geom_point(aes(x = reorder(row,true_cov),y = true_cov), colour = "springgreen4",alpha = 1)+
    geom_errorbar(aes(x = reorder(row,true_cov),ymin = covariance_lwrCI, ymax = covariance_uprCI))+ ggtitle("Covariance - Bootstrap Check")+
    theme_classic() + geom_hline(aes(yintercept = 0))+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())+ 
      facet_wrap(~Covconfintboot,ncol = 2))
 
# Bin Covariance Bootstrap and See whats driving false/true pos's and neg's
dat_csv1$binCov = "NA"#abs(round(dat_csv1$covariance,1))
for(i in 1:nrow(dat_csv1)){
  if(dat_csv1$true_cov[i] == 0){dat_csv$binCov[i] = 0
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


# GxE Perm check
dat_csv1$GxEconfintperm = rep("NA", nrow(dat_csv1))

for(i in 1:nrow(dat_csv1)){
  if(dat_csv1$true_GxE_emm[i] != 0 & dat_csv1$GxE_emm_pvalue[i] <= 0.05){dat_csv1$GxEconfintperm[i] = "True Positive"
  }else if(dat_csv1$true_GxE_emm[i] == 0 & dat_csv1$GxE_emm_pvalue[i] <= 0.05){dat_csv1$GxEconfintperm[i] = "False Positive"
  }else if(dat_csv1$true_GxE_emm[i] != 0 & dat_csv1$GxE_emm_pvalue[i] > 0.05){dat_csv1$GxEconfintperm[i] = "False Negative"
  }else if(dat_csv1$true_GxE_emm[i] == 0 & dat_csv1$GxE_emm_pvalue[i] > 0.05){dat_csv1$GxEconfintperm[i] = "True Negative"
  }else{dat_csv1$GxEconfintperm == "None"}
}

(gxe_perm = ggplot(transform(dat_csv1, GxEconfintperm = factor(GxEconfintperm, levels = c("True Positive", "True Negative", "False Positive", "False Negative")))) +
    geom_point(aes( x = reorder(row,GxE_emm),y = GxE_emm), colour = "firebrick",alpha = 1)+
    geom_point(aes( x = reorder(row,GxE_emm),y = true_GxE_emm), colour = "springgreen4",alpha = 1)+
    geom_errorbar(aes( x = reorder(row,GxE_emm),ymin = GxE_emm_lwrCI, ymax = GxE_emm_uprCI))+ggtitle("GxE - Permutation Check")+
    theme_classic() + geom_hline(aes(yintercept = 0))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+ 
    facet_wrap(~GxEconfintperm,ncol = 2))


# GxE Boot check
dat_csv1$GxEconfintboot = rep("NA", nrow(dat_csv1))

for(i in 1:nrow(dat_csv1)){
  if(dat_csv1$true_GxE_emm[i] != 0 & dat_csv1$GxE_emm_lwrCI[i] > 0)
  {dat_csv1$GxEconfintboot[i] = "True Positive"
  }else if(dat_csv1$true_GxE_emm[i] == 0 & dat_csv1$GxE_emm_lwrCI[i] > 0)
  {dat_csv1$GxEconfintboot[i] = "False Positive"
  }else if(dat_csv1$true_GxE_emm[i] != 0 & dat_csv1$GxE_emm_lwrCI[i] <= 0)
  {dat_csv1$GxEconfintboot[i] = "False Negative"
  }else if(dat_csv1$true_GxE_emm[i] == 0 & dat_csv1$GxE_emm_lwrCI[i] <= 0)
  {dat_csv1$GxEconfintboot[i] = "True Negative"
  }else{dat_csv1$GxEconfintboot[i] = "None"}
}
  
(gxe_boot = ggplot(transform(dat_csv1, GxEconfintboot = factor(GxEconfintboot, levels = c("True Positive", "True Negative", "False Positive", "False Negative")))) +
    geom_point(aes(x = reorder(row,true_GxE_emm),y = GxE_emm), colour = "firebrick",alpha = 1)+
    geom_point(aes(x = reorder(row,true_GxE_emm),y = true_GxE_emm), colour = "springgreen4",alpha = 1)+
    geom_errorbar(aes(x = reorder(row,true_GxE_emm),ymin = GxE_emm_lwrCI, ymax = GxE_emm_uprCI))+ggtitle("GxE - Bootstrap Check")+
    theme_classic() + geom_hline(aes(yintercept = 0))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+ 
    facet_wrap(~GxEconfintboot,ncol = 2))

# Cov Perm check MEANS
dat_csv1$meansCovconfintperm = rep("NA",nrow(dat_csv1))

for(i in 1:nrow(dat_csv1)){
  if(dat_csv1$true_cov_means_correct[i] != 0 && dat_csv1$covariance_pvalue[i] <= 0.025){dat_csv1$meansCovconfintperm[i] = "True Positive"
  }else if(dat_csv1$true_cov_means_correct[i] == 0 & dat_csv1$covariance_pvalue[i] <= 0.025){dat_csv1$meansCovconfintperm[i] = "False Positive"
  }else if(dat_csv1$true_cov_means_correct[i]!= 0 & dat_csv1$covariance_pvalue[i] > 0.025){dat_csv1$meansCovconfintperm[i] = "False Negative"
  }else if(dat_csv1$true_cov_means_correct[i] == 0 & dat_csv1$covariance_pvalue[i] > 0.025){dat_csv1$meansCovconfintperm[i] = "True Negative"
  }else{dat_csv1$meansCovconfintperm[i] = "None"}
}

(means_cov_perm = ggplot(transform(dat_csv1, meansCovconfintperm = factor(meansCovconfintperm, levels = c("True Positive", "True Negative", "False Positive", "False Negative")))) +
    geom_point(aes(x = reorder(row,true_cov_means_correct), y = cov_means_correct), colour = "firebrick",alpha = 1)+
    geom_point(aes(x = reorder(row,true_cov_means_correct), y = true_cov_means_correct), colour = "springgreen4",alpha = 1)+
    geom_errorbar(aes(x = reorder(row,true_cov_means_correct), ymin = cov_means_correct_lwrCI, ymax = cov_means_correct_uprCI))+ ggtitle("Means: Covariance - Permutation Check")+
    theme_classic() + geom_hline(aes(yintercept = 0))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+ 
    facet_wrap(~meansCovconfintperm,ncol = 2))

# MEANS Cov Boot check
dat_csv1$MeansCovconfintboot = rep(NA, nrow(dat_csv1))

for(i in 1:nrow(dat_csv1)){
  
  if(dat_csv1$true_cov[i] != 0 &&
     dat_csv1$covariance_lwrCI[i] < 0 &&
     dat_csv1$covariance_uprCI[i] < 0
  ){dat_csv1$MeansCovconfintboot[i] = "True Positive"
  }else if(dat_csv1$true_cov[i] != 0 &&
           dat_csv1$covariance_lwrCI[i] > 0 &&
           dat_csv1$covariance_uprCI[i] > 0
  ){dat_csv1$MeansCovconfintboot[i] = "True Positive"
  }else if(dat_csv1$true_cov[i] == 0 &&
           dat_csv1$covariance_lwrCI[i] < 0 &&
           dat_csv1$covariance_uprCI[i] < 0
  ){dat_csv1$MeansCovconfintboot[i] = "False Positive"
  }else if(dat_csv1$true_cov[i] == 0 &&
           dat_csv1$covariance_lwrCI[i] > 0 &&
           dat_csv1$covariance_uprCI[i] > 0
  ){dat_csv1$MeansCovconfintboot[i] = "False Positive"
  }else if(dat_csv1$true_cov[i] != 0 && 
           dat_csv1$covariance_lwrCI[i] <= 0 && 
           dat_csv1$covariance_uprCI[i] >= 0
  ){dat_csv1$MeansCovconfintboot[i] = "False Negative"
  }else if(dat_csv1$true_cov[i]== 0 && 
           dat_csv1$covariance_lwrCI[i] <= 0 && 
           dat_csv1$covariance_uprCI[i] >= 0
  ){dat_csv1$MeansCovconfintboot[i] = "True Negative"
  }else{dat_csv1$MeansCovconfintboot[i] = "None"}
  
}

(cov_boot = ggplot(transform(dat_csv1, MeansCovconfintboot = factor(MeansCovconfintboot, levels = c("True Positive", "True Negative", "False Positive", "False Negative")))) +
    geom_point(aes(x = reorder(row,true_cov_means_correct), y = cov_means_correct), colour = "firebrick",alpha = 1)+
    geom_point(aes(x = reorder(row,true_cov_means_correct), y = true_cov_means_correct), colour = "springgreen4",alpha = 1)+
    geom_errorbar(aes(x = reorder(row,true_cov_means_correct), ymin = cov_means_correct_lwrCI, ymax = cov_means_correct_uprCI))+ ggtitle("Means: Covariance - Bootstrap Check")+
    theme_classic() + geom_hline(aes(yintercept = 0))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+ 
    facet_wrap(~MeansCovconfintboot,ncol = 2))

# MEANS GxE Boot check
dat_csv1$MeanGxEconfintboot = rep("NA", nrow(dat_csv1))

for(i in 1:nrow(dat_csv1)){
  if(dat_csv1$true_GxE_emm[i] != 0 & dat_csv1$GxE_emm_lwrCI[i] > 0)
  {dat_csv1$MeanGxEconfintboot[i] = "True Positive"
  }else if(dat_csv1$true_GxE_emm[i] == 0 & dat_csv1$GxE_emm_lwrCI[i] > 0)
  {dat_csv1$MeanGxEconfintboot[i] = "False Positive"
  }else if(dat_csv1$true_GxE_emm[i] != 0 & dat_csv1$GxE_emm_lwrCI[i] <= 0)
  {dat_csv1$MeanGxEconfintboot[i] = "False Negative"
  }else if(dat_csv1$true_GxE_emm[i] == 0 & dat_csv1$GxE_emm_lwrCI[i] <= 0)
  {dat_csv1$MeanGxEconfintboot[i] = "True Negative"
  }else{dat_csv1$MeanGxEconfintboot[i] = "None"}
}

(mean_gxe_boot = ggplot(transform(dat_csv1, MeanGxEconfintboot = factor(MeanGxEconfintboot, levels = c("True Positive", "True Negative", "False Positive", "False Negative")))) +
    geom_point(aes(x = reorder(row,true_GxE_means),y = GxE_means), colour = "firebrick",alpha = 1)+
    geom_point(aes(x = reorder(row,true_GxE_means),y = true_GxE_means), colour = "springgreen4",alpha = 1)+
    geom_errorbar(aes(x = reorder(row,true_GxE_means),ymin = GxE_means_lwrCI, ymax = GxE_means_uprCI))+ggtitle("Means: GxE - Bootstrap Check")+
    theme_classic() + geom_hline(aes(yintercept = 0))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+ 
    facet_wrap(~MeanGxEconfintboot,ncol = 2))
  
# MEANS GxE Perm check
dat_csv1$meansGxEconfintperm = rep("NA", nrow(dat_csv1))

for(i in 1:nrow(dat_csv1)){
  if(dat_csv1$true_GxE_emm[i] != 0 & dat_csv1$GxE_emm_pvalue[i] <= 0.05){dat_csv1$meansGxEconfintperm[i] = "True Positive"
  }else if(dat_csv1$true_GxE_emm[i] == 0 & dat_csv1$GxE_emm_pvalue[i] <= 0.05){dat_csv1$meansGxEconfintperm[i] = "False Positive"
  }else if(dat_csv1$true_GxE_emm[i] != 0 & dat_csv1$GxE_emm_pvalue[i] > 0.05){dat_csv1$meansGxEconfintperm[i] = "False Negative"
  }else if(dat_csv1$true_GxE_emm[i] == 0 & dat_csv1$GxE_emm_pvalue[i] > 0.05){dat_csv1$meansGxEconfintperm[i] = "True Negative"
  }else{dat_csv1$meansGxEconfintperm == "None"}
}

(means_gxe_perm = ggplot(transform(dat_csv1, meansGxEconfintperm = factor(meansGxEconfintperm, levels = c("True Positive", "True Negative", "False Positive", "False Negative")))) +
    geom_point(aes( x = reorder(row,true_GxE_means),y = GxE_means), colour = "firebrick",alpha = 1)+
    geom_point(aes( x = reorder(row,true_GxE_means),y = true_GxE_means), colour = "springgreen4",alpha = 1)+
    geom_errorbar(aes( x = reorder(row,true_GxE_means),ymin = GxE_means_lwrCI, ymax = GxE_means_uprCI))+ggtitle("Means: GxE - Permutation Check")+
    theme_classic() + geom_hline(aes(yintercept = 0))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+ 
    facet_wrap(~meansGxEconfintperm,ncol = 2))              

# ANOVA vs. GxE_EMM comparison (to show that prop variance explained doesn't equal effect size of GxE)
(ggplot(dat_csv, aes(x = true_GxE_emm, y = GxE_omega, group = factor(std_dev), colour = factor(std_dev))) + geom_point()+
    geom_smooth(method = "lm", formula = y ~ splines::bs(x, 3),se = F) + 
    geom_abline(slope = 1, intercept = 0,colour = "black")+theme_classic()+
    scale_colour_manual(values = c("0.5" = "#6600CC", "1.5" = "#66BBAA"))+
    xlab("True GxE - Estimated Marginal Means")+ylab("GxE - Omega Squared")+
    labs(colour = "Standard Deviation"))


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

######### Test out some GxE
G = c(1,2)
phen1 = 2 + (-1*G)  aov.test <- lm(phen_corrected ~ exp_env_factor * gen_factor, data = input_df) 

# Estimated Marginal Means
emm_E = as.data.frame(emmeans(aov.test,"exp_env_factor"))
emm_G = as.data.frame(emmeans(aov.test, "gen_factor"))
emm_GxE = as.data.frame(emmeans(aov.test, ~ exp_env_factor*gen_factor))

# Gmeans
#E_means <- tapply(emm_GxE$emmean, emm_GxE$exp_env_factor, mean)
#G_means <- tapply(emm_GxE$emmean, emm_GxE$gen_factor, mean)
G_matrix <- data.frame("G_means" = emm_G$emmean, "gen_factor" = emm_G$gen_factor)
E_matrix <- data.frame("E_means" = emm_E$emmean, "exp_env_factor" = emm_E$exp_env_factor)

# Match Genotypes to Native Environment
Cov_matrix = G_matrix
Cov_matrix$exp_env_factor <- input_df$nat_env_factor[match(G_matrix$gen_factor,input_df$gen_factor)]
Cov_matrix$E_means <- E_matrix$E_means[match(Cov_matrix$exp_env_factor,E_matrix$exp_env_factor)]

# Magnitude of GxE -- EMMs
GxE_emm_original<- abs(emm_GxE$emmean[emm_GxE$gen_factor == "G_1" & emm_GxE$exp_env_factor == "E_1"] - # GxE (Phenotype of ith genotype in jth environment)
                         emm_G$emmean[emm_G$gen_factor == "G_1"] - # phenotype of ith Genotype
                         emm_E$emmean[emm_E$exp_env_factor == "E_1"] + # phenotype of jth Environment
                         mean(emm_GxE$emmean)) 
phen2 = -2 + (1*G)
phen3 = 0 + (0*G)
phen = rbind (phen1, phen2, phen3)


playdf = data.frame(gen_factor = rep(c("G_1","G_2","G_3"), each = 3),
                    exp_env_factor = rep(c("E_1","E_2","E_3"), times = 3), 
                    phen = c(-2,0,2,2,0,-2,0,0,0))
playdf = data.frame(gen_factor = rep(c("G_1","G_2"), each = 2),
                    exp_env_factor = rep(c("E_1","E_2"), times = 2), 
                    phen = c(1,-.5,-1,0.5))
  
ggplot(playdf, aes(x = exp_env_factor, y = phen, group = gen_factor))+theme_hc()+ geom_line()


aov.test <- lm(phen ~ exp_env_factor * gen_factor, data = playdf) 

# Estimated Marginal Means
emm_E = as.data.frame(emmeans(aov.test,"exp_env_factor"))
emm_G = as.data.frame(emmeans(aov.test, "gen_factor"))
emm_GxE = as.data.frame(emmeans(aov.test, ~ exp_env_factor*gen_factor))


# Magnitude of GxE -- Loop
allGE = c()
loopGxE = NULL
for (i in 1:nlevels(emm_GxE$gen_factor)){
  for (j in 1:nlevels(emm_GxE$exp_env_factor)){
    G_levels <- levels(emm_GxE$gen_factor)
    E_levels <- levels(emm_GxE$exp_env_factor)
    loopGxE <- abs(emm_GxE$emmean[emm_GxE$gen_factor == G_levels[i] & emm_GxE$exp_env_factor == E_levels[j]] - # GxE (Phenotype of ith genotype in jth environment)
                     emm_G$emmean[emm_G$gen_factor == G_levels[i]] - # mean phenotype of ith Genotype
                     emm_E$emmean[emm_E$exp_env_factor == E_levels[j]] + # mean phenotype of jth Environment
                     mean(emm_GxE$emmean)) # Overall mean
    allGE <- c(allGE, loopGxE)
  }
}
#hist(allGE)
(GxE_emm_loop = mean(allGE))

