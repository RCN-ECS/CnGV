##########################################################
###             Simulation Visualization Code          ###
##########################################################

## Load packages
library(ggplot2)
library(readr)
library(tidyverse)
library(gridExtra)
library(ggthemes)

# Load Data compiled on cluster
dat_csv1 = read.csv("~/Desktop/Power_output_results.csv")
#dat_csv1 = filter(dat_csv1, row %in% c(1:31900))
####### Check Average time for longest sims ##########
(timecheck = dat_csv1 %>%
  #filter(n_pop == max(n_pop)) %>%
  #filter(std_dev == max(std_dev)) %>%
  #filter(sample_size == max(sample_size)) %>%
  summarize(average_time = mean(Sim_time)))

######################
## Covariance x GxE ##
######################

dat_csv = filter(dat_csv1, env_scenario == 1)
dat_dub = filter(dat_csv1, env_scenario == 2)

####### Check how many results per N-Pop/sample size ##########
(sizecheck = dat_csv1 %>%
   filter(n_pop == min(n_pop)) %>%
   filter(sample_size == min(sample_size)) %>%
   #filter(n_pop == max(n_pop)) %>%
   #filter(sample_size == max(sample_size)) %>%
   summarize(size = n()))

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
ggplot(dat_csv, aes(x = covariance, y = GxE_emm, group = factor(n_pop), alpha = 0.1,colour = col)) + 
  geom_jitter() + theme_classic() + ylim(0,1) + xlim(-1,1)+
  xlab("Covariance Estimate") + ylab("GxE Estimate") +
  theme(legend.position = "none")+
  scale_colour_identity()+
  facet_grid(sample_size~n_pop)

 ggplot(dat_dub, aes(x = covariance, y = true_GxE_emm, group = factor(n_pop), alpha = 0.1,colour = col)) + 
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

ggplot((filter(dat_dub, sample_size != 16)), aes(x = covariance, y = GxE_omega, group = factor(n_pop), alpha = 0.1,colour = col)) + 
  geom_jitter() + theme_classic() + 
  xlab("Covariance Estimate") + ylab("GxE Estimate (Omega^2)") +
  theme(legend.position = "none")+
  scale_colour_identity()+
  facet_grid(sample_size~n_pop)


#############################
## Power Analysis Heatmaps ##
#############################

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

# Power at moderate levels of Cov and GxE -- Scenario 1
summary(abs(dat_csv$covariance))
summary(dat_csv$GxE_emm)
covpow1 = dat_csv %>% filter(between(abs(covariance),0.15,0.5)) #Between 1st and 3rd quartile ish
gxepow1 = dat_csv %>% filter(between(GxE_emm,0.2,0.7)) 

covpow = covpow1 %>%
  group_by(sample_size, n_pop, std_dev) %>%
  summarize("total_tick" = n(),
            "covtick" = sum(covtick))
covpow$covpower = covpow$covtick/covpow$total_tick

gxepow = gxepow1 %>%
  group_by(sample_size,n_pop,std_dev) %>%
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

# Plot 

# Set midrange for plotting
rng.gxe = range(gxepow_highsd$meangxepower, gxepow_lowsd$meangxepower)
rng.cov = range(covpow_highsd$meancovpower, covpow_lowsd$meancovpower)


(csdlow = ggplot(covpow_lowsd,aes(x = factor(sample_size), y = factor(n_pop), fill = meancovpower)) + 
  geom_tile() + scale_fill_gradient2(low="lemon chiffon", mid="#00AA00", high="#003300", #colors in the scale
                         midpoint=mean(rng.cov),    #same midpoint for plots (mean of the range)
                         breaks=seq(0,1,0.25), #breaks in the scale bar
                         limits=c(0,1))+
  xlab("Sample Size") + ylab("Number of Populations")+
  labs(fill = "Power")+
  theme_classic() + ggtitle(paste0("Covariance: Standard Deviation = ",min(covpow$std_dev))))

(csdhigh = ggplot(covpow_highsd,aes(x = factor(sample_size), y = factor(n_pop), fill = meancovpower)) + 
  geom_tile() + scale_fill_gradient2(low="lemon chiffon", mid="#00AA00", high="#003300", #colors in the scale
                                    midpoint=mean(rng.cov),    #same midpoint for S (mean of the range)
                                    breaks=seq(0,1,0.25), #breaks in the scale bar
                                    limits=c(0,1))+
  xlab("Sample Size") + ylab("Number of Populations")+
  labs(fill = "Power")+
  theme_classic()+ ggtitle(paste0("Covariance: Standard Deviation = ",max(covpow$std_dev))))

(gsdlow = ggplot(gxepow_lowsd,aes(x = factor(sample_size), y = factor(n_pop), fill = meangxepower)) + 
  geom_tile() + scale_fill_gradient2(low="#DDDDDD", mid="#99CCEE", high="#000044", #colors in the scale
                                     midpoint=mean(rng.gxe),    #same midpoint for plots (mean of the range)
                                     breaks=seq(0,1,0.25), #breaks in the scale bar
                                     limits=c(floor(rng.gxe[1]), ceiling(rng.gxe[2])))+
  xlab("Sample Size") + ylab("Number of Populations")+
  labs(fill = "Power")+
  theme_classic() + ggtitle(paste0("GxE: Standard Deviation = ",min(gxepow$std_dev))))

(gsdhigh = ggplot(gxepow_highsd,aes(x = factor(sample_size), y = factor(n_pop), fill = meangxepower)) + 
  geom_tile() + scale_fill_gradient2(low="#DDDDDD", mid="#99CCEE", high="#000044", #colors in the scale
                                     midpoint=mean(rng.gxe),    #same midpoint for plots (mean of the range)
                                     breaks=seq(0,1,0.25), #breaks in the scale bar
                                     limits=c(floor(rng.gxe[1]), ceiling(rng.gxe[2])))+
  xlab("Sample Size") + ylab("Number of Populations")+
  labs(fill = "Power") +
  theme_classic() + ggtitle(paste0("GxE: Standard Deviation = ",max(gxepow$std_dev))))

library(gridExtra)
grid.arrange(gsdhigh,gsdlow,csdhigh,csdlow,ncol = 2)

# Power at moderate levels of Cov and GxE for Env_scenario = 2

covpow1 = data.frame(dat_dub)
gxepow1 = dat_dub %>% filter(between(GxE_emm,0,1))

covpow = covpow1 %>%
  group_by(sample_size, n_pop, std_dev) %>%
  summarize("total_tick" = n(),
            "covtick" = sum(covtick))
covpow$covpower = covpow$covtick/covpow$total_tick

gxepow = gxepow1 %>%
  group_by(sample_size,n_pop,std_dev) %>%
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

# Plot Power

(csdlow = ggplot(covpow_lowsd,aes(x = factor(sample_size), y = factor(n_pop), fill = meancovpower)) + 
    geom_tile() + scale_fill_gradient2(low="lemon chiffon", mid="#00AA00", high="#003300", #colors in the scale
                                       midpoint=mean(rng.cov),    #same midpoint for S (mean of the range)
                                       breaks=seq(0,1,0.25), #breaks in the scale bar
                                       limits=c(0,1))+
    xlab("Sample Size") + ylab("Number of Populations")+
    labs(fill = "Power")+
    theme_classic() + ggtitle(paste0("Covariance: Standard Deviation = ",min(covpow$std_dev))))

(csdhigh = ggplot(covpow_highsd,aes(x = factor(sample_size), y = factor(n_pop), fill = meancovpower)) + 
    geom_tile() + scale_fill_gradient2(low="lemon chiffon", mid="#00AA00", high="#003300", #colors in the scale
                                       midpoint=mean(rng.cov),    #same midpoint for S (mean of the range)
                                       breaks=seq(0,1,0.25), #breaks in the scale bar
                                       limits=c(0,1))+
    xlab("Sample Size") + ylab("Number of Populations")+
    labs(fill = "Power")+
    theme_classic()+ ggtitle(paste0("Covariance: Standard Deviation = ",max(covpow$std_dev))))

(gsdlow = ggplot(gxepow_lowsd,aes(x = factor(sample_size), y = factor(n_pop), fill = meangxepower)) + 
    geom_tile() + scale_fill_gradient2(low="#DDDDDD", mid="#99CCEE", high="#000044", #colors in the scale
                                       midpoint=mean(rng.gxe),    #same midpoint for plots (mean of the range)
                                       breaks=seq(0,1,0.25), #breaks in the scale bar
                                       limits=c(floor(rng.gxe[1]), ceiling(rng.gxe[2])))+
    xlab("Sample Size") + ylab("Number of Populations")+
    labs(fill = "Power")+
    theme_classic() + ggtitle(paste0("GxE: Standard Deviation = ",min(gxepow$std_dev))))

(gsdhigh = ggplot(gxepow_highsd,aes(x = factor(sample_size), y = factor(n_pop), fill = meangxepower)) + 
    geom_tile() + scale_fill_gradient2(low="#DDDDDD", mid="#99CCEE", high="#000044", #colors in the scale
                                       midpoint=mean(rng.gxe),    #same midpoint for plots (mean of the range)
                                       breaks=seq(0,1,0.25), #breaks in the scale bar
                                       limits=c(floor(rng.gxe[1]), ceiling(rng.gxe[2])))+
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

(barpower = ggplot(new_df,aes(x = factor(n_pop), y = proportion, group = sample_size, fill = factor(col,levels=c("grey","dodgerblue4","darkgreen","red"))))+
  geom_bar(position = "stack", stat="identity", alpha = 0.9) + ylab("Proportion of simulation results")+ xlab("Number of Populations")+
    scale_fill_identity() +
    facet_wrap(~sample_size) + theme_classic())

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

# Line-style Power plots

dat_dub$mod = NULL
for(i in 1:nrow(dat_dub)){
  if(dat_dub$covariance[i] >= 0.4 & dat_dub$covariance[i] <= 0.6){dat_dub$mod[i] = TRUE
  }else if(dat_dub$GxE_emm[i] >= 0.4 & dat_dub$GxE_emm[i] <= 0.6){dat_dub$mod[i] = TRUE
  }else{dat_dub$mod[i] = FALSE}
}

col_df = dat_dub %>%
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

(barpower2 = ggplot(new_df,aes(x = factor(n_pop), y = proportion, group = sample_size, fill = factor(col,levels=c("grey","dodgerblue4","darkgreen","red"))))+
    geom_bar(position = "stack", stat="identity", alpha = 0.9) + ylab("Proportion of simulation results")+ xlab("Number of Populations")+
    scale_fill_identity() +
    facet_wrap(~sample_size) + theme_classic())

# Tradeoff between GxE and Covariance plot

dat_csv$sig = NULL
for(i in 1:nrow(dat_csv)){ # Use only if one or the other is significant
if(dat_csv$GxE_emm_pvalue[i] <= 0.05 | dat_csv$covariance_lwrCI[i] > 0 & dat_csv$covariance_uprCI[i] > 0){dat_csv$sig[i] = TRUE
}else if(dat_csv$GxE_emm_pvalue[i] <= 0.05 | dat_csv$covariance_lwrCI[i] < 0 & dat_csv$covariance_uprCI[i] < 0){dat_csv$sig[i] = TRUE
}else{dat_csv$sig[i]=FALSE} 
  }

sigGxE = filter(dat_csv, sig ==TRUE) # filter out false positives as potential solution to weed out messiness.


ggplot(sigGxE, aes(x = GxE_emm, y = covtick))+
  geom_smooth(method = "glm",method.args = list(family = "binomial"),se = T,colour = "black") + 
  geom_point()+
  xlab("Magnitude of GxE")+ylab("Proportion of significant CovGE values")+
  theme_bw(base_size = 18, base_family = "Times")+
  theme(axis.text.x = element_text(size=16,colour = "black"),
        axis.title.x = element_text(size=18,face="bold")) +
  theme(axis.text.y = element_text(size=16,colour = "black"),
        axis.title.y = element_text(size=18,face="bold")) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2))


ggplot(sigGxE, aes(x = abs(true_cov), y = true_GxE_emm))+
  #geom_smooth(method = "glm",method.args = list(family = "binomial"),se = T,colour = "red") + 
  geom_point(alpha = 0.2)+
  geom_smooth(method = "glm")+
  
  ylab("Actual GxE")+xlab(expression("True Cov"["GE"]))+
  theme_bw(base_size = 18, base_family = "Times")+
  theme(axis.text.x = element_text(size=16,colour = "black"),
        axis.title.x = element_text(size=18,face="bold")) +
  theme(axis.text.y = element_text(size=16,colour = "black"),
        axis.title.y = element_text(size=18,face="bold")) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2))

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

######################################
##          Sanity Checks           ##
######################################

# Do we have good parameter coverage? Hex plot

(hexy = ggplot(dat_csv, aes(x = covariance, y = GxE_emm)) + 
  geom_hex()+ theme_classic() + xlab("Covariance Estimate") + ylab("GxE Estimate") +
  ggtitle("HexPlot") + facet_grid(sample_size~n_pop))

(hexy2 = ggplot(dat_dub, aes(x = covariance, y = GxE_emm)) + 
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
    geom_errorbar(aes(x = reorder(row,true_cov), ymin = covariance_lwrCI, ymax = covariance_uprCI),alpha = 0.5)+ 
    geom_point(aes(x = reorder(row,true_cov), y = covariance), shape = 24,colour = "black", fill = "#7570B3",size = 1.2)+
    geom_point(aes(x = reorder(row,true_cov), y = true_cov), shape = 21, colour = "black", fill = "#E6AB02",size = 1.2)+
    theme_classic() + geom_hline(aes(yintercept = 0))+ ggtitle("Covariance - Permutation Check")+
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
    geom_errorbar(aes(x = reorder(row,true_cov), ymin = covariance_lwrCI, ymax = covariance_uprCI),alpha = 0.5)+ 
    geom_point(aes(x = reorder(row,true_cov), y = covariance), shape = 24,colour = "black", fill = "#7570B3",size = 1.2)+
    geom_point(aes(x = reorder(row,true_cov), y = true_cov), shape = 21, colour = "black", fill = "#E6AB02",size = 1.2)+ ggtitle("Covariance - Bootstrap Check")+
    theme_classic() + geom_hline(aes(yintercept = 0))+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())+ 
      facet_wrap(~Covconfintboot,ncol = 2))
 
# Bin Covariance Bootstrap and See whats driving false/true pos's and neg's
dat_csv1$binCov = "NA"#abs(round(dat_csv1$covariance,1))
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
    geom_errorbar(aes(x = reorder(row,GxE_emm), ymin = GxE_emm_lwrCI, ymax = GxE_emm_uprCI),alpha = 0.5)+ 
    geom_point(aes(x = reorder(row,GxE_emm), y = GxE_emm), shape = 24,colour = "black", fill = "#7570B3",size = 1.2)+
    geom_point(aes(x = reorder(row,GxE_emm), y = true_GxE_emm), shape = 21, colour = "black", fill = "#E6AB02",size = 1.2)+
    ggtitle("GxE - Permutation Check")+
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
    geom_errorbar(aes(x = reorder(row,GxE_emm), ymin = GxE_emm_lwrCI, ymax = GxE_emm_uprCI),alpha = 0.5)+ 
    geom_point(aes(x = reorder(row,GxE_emm), y = GxE_emm), shape = 24,colour = "black", fill = "#7570B3",size = 1.2)+
    geom_point(aes(x = reorder(row,GxE_emm), y = true_GxE_emm), shape = 21, colour = "black", fill = "#E6AB02",size = 1.2)+ggtitle("GxE - Bootstrap Check")+
    theme_classic() + geom_hline(aes(yintercept = 0))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+ 
    facet_wrap(~GxEconfintboot,ncol = 2))

# Cov Perm check MEANS
dat_csv1$meansCovconfintperm = rep("NA",nrow(dat_csv1))

for(i in 1:nrow(dat_csv1)){
  if(dat_csv1$true_cov_means[i] != 0 && dat_csv1$covariance_pvalue[i] <= 0.025){dat_csv1$meansCovconfintperm[i] = "True Positive"
  }else if(dat_csv1$true_cov_means[i] == 0 & dat_csv1$covariance_pvalue[i] <= 0.025){dat_csv1$meansCovconfintperm[i] = "False Positive"
  }else if(dat_csv1$true_cov_means[i]!= 0 & dat_csv1$covariance_pvalue[i] > 0.025){dat_csv1$meansCovconfintperm[i] = "False Negative"
  }else if(dat_csv1$true_cov_means[i] == 0 & dat_csv1$covariance_pvalue[i] > 0.025){dat_csv1$meansCovconfintperm[i] = "True Negative"
  }else{dat_csv1$meansCovconfintperm[i] = "None"}
}

(means_cov_perm = ggplot(transform(dat_csv1, meansCovconfintperm = factor(meansCovconfintperm, levels = c("True Positive", "True Negative", "False Positive", "False Negative")))) +
    geom_point(aes(x = reorder(row,true_cov_means), y = cov_means), colour = "firebrick",alpha = 1)+
    geom_point(aes(x = reorder(row,true_cov_means), y = true_cov_means), colour = "springgreen4",alpha = 1)+
    geom_errorbar(aes(x = reorder(row,true_cov_means), ymin = cov_means_lwrCI, ymax = cov_means_uprCI))+ ggtitle("Means: Covariance - Permutation Check")+
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
    geom_errorbar(aes(x = reorder(row,true_cov_means), ymin = cov_means_lwrCI, ymax = cov_means_uprCI),alpha = 0.5)+ 
    geom_point(aes(x = reorder(row,true_cov_means), y = cov_means), shape = 24,colour = "black", fill = "#7570B3",size = 1.2)+
    geom_point(aes(x = reorder(row,true_cov_means), y = true_cov_means), shape = 21, colour = "black", fill = "#E6AB02",size = 1.2)+
   ggtitle("Means: Covariance - Bootstrap Check")+
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
    geom_errorbar(aes(x = reorder(row,true_GxE_means), ymin = GxE_means_lwrCI, ymax = GxE_means_uprCI),alpha = 0.5)+ 
    geom_point(aes(x = reorder(row,true_GxE_means), y = GxE_means), shape = 24,colour = "black", fill = "#7570B3",size = 1.2)+
    geom_point(aes(x = reorder(row,true_GxE_means), y = true_GxE_means), shape = 21, colour = "black", fill = "#E6AB02",size = 1.2)+
    ggtitle("Means: GxE - Bootstrap Check")+
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
    geom_errorbar(aes(x = reorder(row,true_GxE_means), ymin = GxE_means_lwrCI, ymax = GxE_means_uprCI),alpha = 0.5)+ 
    geom_point(aes(x = reorder(row,true_GxE_means), y = GxE_means), shape = 24,colour = "black", fill = "#7570B3",size = 1.2)+
    geom_point(aes(x = reorder(row,true_GxE_means), y = true_GxE_means), shape = 21, colour = "black", fill = "#E6AB02",size = 1.2)+
    ggtitle("Means: GxE - Permutation Check")+
    theme_classic() + geom_hline(aes(yintercept = 0))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+ 
    facet_wrap(~meansGxEconfintperm,ncol = 2))              

# ANOVA vs. GxE_EMM comparison (to show that prop variance explained doesn't equal effect size of GxE)
ggplot(dat_csv1, aes(x = true_GxE_emm, y = GxE_omega, group = factor(std_dev))) + 
    geom_point(aes(shape = factor(n_pop), fill = factor(std_dev)))+
    
    geom_abline(slope = 1, intercept = 0,colour = "black")+theme_classic()+
    
    scale_fill_manual(values = c("0.5" = "#1B9E77", "1" = "#7570B3"),
                      breaks = c("0.5", "1"))+
  scale_shape_manual(values = c("2" = 21, "4" = 22,"8" = 23))+
  geom_smooth(method = "lm", formula = y ~ splines::bs(x, 3),se = F,colour = "black") + 
    xlab("Actual GxE - Estimated Marginal Means")+ylab("GxE - Omega Squared")+
    labs(fill = "Standard Deviation",shape = "Number of Populations")


# Do confidence intervals for each overlap with their true values?  
dat_csv1$probgxe = NULL
dat_csv1$probcov = NULL

for(i in 1:nrow(dat_csv1)){
  if((dat_csv1$true_GxE_emm[i] < dat_csv1$GxE_emm_lwrCI[i]) | (dat_csv1$true_GxE_emm[i] >dat_csv1$GxE_emm_uprCI[i])){dat_csv1$probgxe[i] = TRUE
  }else{dat_csv1$probgxe[i] = FALSE}
  if((dat_csv1$true_cov[i] < dat_csv1$covariance_lwrCI[i]) | (dat_csv1$true_cov[i] >dat_csv1$covariance_uprCI[i])){dat_csv1$probcov[i] = TRUE
  }else{dat_csv1$probcov[i] = FALSE}
}

GxEanoms = dat_csv1 %>%
  filter(probgxe == TRUE)
dim(GxEanoms)
length(which(GxEanoms$GxE_emm_pvalue<=0.05)) 

cov_anom = dat_csv1 %>%
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


# Do means estimates match raw estimates? If yes, should fall along 1:1 line

dat_csv1$meancoverror = abs(dat_csv1$cov_means_uprCI - dat_csv1$cov_means_lwrCI)
dat_csv1$coverror = abs(dat_csv1$covariance_uprCI - dat_csv1$covariance_lwrCI)
dat_csv1$meangxeerror = dat_csv1$GxE_means_uprCI - dat_csv1$GxE_means_lwrCI
dat_csv1$gxeerror = dat_csv1$GxE_emm_uprCI - dat_csv1$GxE_emm_lwrCI

# Covariance
(covmeancheck = ggplot(dat_csv1,aes(x = covariance, y = cov_means))+
    geom_point()+theme_classic()+ylab("Covariance from Means")+xlab("Covariance from Raw")+
    geom_abline(slope = 1, intercept = 0,colour = "red"))

(covcheck = ggplot(dat_csv1,aes(x = true_cov, y = covariance, colour = Covconfintboot))+ # Cov = 1 happens when delta_env = delta_gen, they are both positive, and there is no interaction
    geom_point()+theme_classic()+
    ylab("Covariance")+xlab("True Covariance")+scale_color_viridis(discrete = TRUE)+ # In the no error scenarios, this makes perfectly parallel lines that always have true_cov = 1
    labs(colour = "Confusion Matrix based on CIs")+geom_abline(slope = 1, intercept = 0,colour = "red"))

(coverrorcheck = ggplot(dat_csv1, aes(x = coverror,y = meancoverror)) + 
    geom_point(alpha = 0.5) + theme_classic() + ylab("Length of 95% CI for CovGE means") + xlab("Length of 95% CI for CovGE raw")+
    geom_abline(slope = 1, intercept = 0,colour = "red")+scale_colour_identity())

suspect = dat_csv %>% filter(true_cov == 1) %>% filter(covariance_pvalue <0.05) # Check on those weirdos

# GxE
(gxemeancheck = ggplot(dat_csv1,aes(x = true_GxE_emm, y = true_GxE_means))+
    geom_point()+theme_classic()+ylab("GxE from Means")+xlab("GxE from Raw")+
    geom_abline(slope = 1, intercept = 0,colour = "red"))

(gxecheck = ggplot(dat_csv1,aes(x = true_GxE_emm, y = GxE_emm, colour = factor(GxEconfintperm)))+
    geom_point()+#(aes(colour = col))+#geom_smooth(aes(group = n_pop,colour = "black",linetype = factor(n_pop)),se = FALSE,method = "glm")+
    theme_classic()+ylab("GxE with error")+xlab("True GxE") +scale_color_viridis(discrete = TRUE)+
    labs(colour = "Confusion Matrix based on p-value")+
    geom_abline(slope = 1, intercept = 0,colour = "red"))#+facet_wrap(~n_pop)

(gxeerrorcheck = ggplot(dat_csv1, aes(x = gxeerror,y = meangxeerror)) + 
    geom_point(alpha = 0.3) + theme_classic() + ylab("Length of 95% CI for GxE means") + xlab("Length of 95% CI for GxE raw")+
    geom_abline(slope = 1, intercept = 0,colour = "red")+scale_colour_identity())


# Do GxE pvalues from permutation match pvalues from anova? 
dat_csv1$sanityColor = NULL
for(i in 1:nrow(dat_csv1)){
  if(dat_csv1$true_GxE_emm[i] == 0){dat_csv1$sanityColor[i] = "#FFBB55"}else{dat_csv1$sanityColor[i] ="#660099"}
}


(ggplot(dat_csv1,aes(x = GxE_Anova, y = GxE_emm_pvalue,shape = factor(n_pop), colour = sanityColor))+ #raw data 
    geom_jitter()+theme_classic()+ylab("GxE EMM Pvalue")+xlab("GxE Anova Pvalue")+
    scale_colour_identity()+labs(shape = "Number of Populations")+
    geom_vline(xintercept = 0.05,colour = "red")+
    geom_hline(yintercept = 0.05,colour = "red")+labs(shape = "N_pop"))

(ggplot(dat_csv1,aes(x = GxE_Anova, y = GxE_means_pvalue,colour = GxE_means,shape = factor(n_pop)))+ # means data
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

# Make Tables of Confusion Matrix
(cov_confint = dat_csv1 %>%
  group_by(Covconfintperm) %>%
  summarize(frequency = n()))

(cov_boot = dat_csv1 %>%
  group_by(Covconfintboot) %>%
  summarize(frequency = n()))

(gxe_confint = dat_csv1 %>%
  group_by(GxEconfintperm) %>%
  summarize(frequency = n()))

(GxE_boot = dat_csv1 %>%
  group_by(GxEconfintboot) %>%
  summarize(frequency = n()))

#### For Heuristic Plots in paper

shape1 = c("G_1" = 24, "G_2" = 21)
ColorFill = c("G_1" = "#002266", "G_2" = "#CC6600")

#model_df$phen_corrected = rep(c(-1,1,-.95,1.05),each = 5) For Vg, Ve, Vgxe
# Otherwise use simulated data from Cov_GxE_clusterFun.R
ggplot(model_df, aes(x = exp_env_factor, y = phen_corrected, group = gen_factor, shape = gen_factor, fill = gen_factor, colour = gen_factor)) + 
  geom_point(size = 5) + geom_line(size=1,lty = "dashed") + 
  scale_shape_manual(values = shape1,
                     labels = c("G_1" = "Genotype 1", "G_2" = "Genotype 2"))+
  scale_fill_manual(values = ColorFill,
                    labels = c("G_1" = "Genotype 1", "G_2" = "Genotype 2"))+
  scale_colour_manual(values = ColorFill,
                      labels = c("G_1" = "Genotype 1", "G_2" = "Genotype 2"))+
  ylab("Phenotype") + xlab(" ")+
  theme_classic(base_size = 24, base_family = "Times")+ 
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
