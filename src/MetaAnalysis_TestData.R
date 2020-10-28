
# Load functions
source("~/Documents/GitHub/CnGV/src/Cov_GxE_functions.R")

# Molly's age at metamorphosis data
mm = read.csv("~/Desktop/Work/DataSets/Tadpole Plasticity_2017/mortality2017.csv")
mm = mm[-which(is.na(mm$Jul_metamorph)),] # Include only those that metamorphosed
mm$age = mm$Jul_metamorph-mm$Jul_hatch # Calculate age at MM

# Rename variables for analysis
mm1 = mm %>%
  filter(Pop != "BELL") %>% # Exclude bellamy bc missing data at 6ppt
  filter(tad %in% c(0,6)) %>% # Only use 0, 6ppt
  droplevels()

mm1$gen_factor = paste0("G_",as.numeric(mm1$Pop))
mm1$exp_env_factor = paste0("E_",as.numeric(as.factor(mm1$tad))) 
mm1$nat_env_factor = NULL
for(i in 1:nrow(mm1)){
  if(mm1$Pop[i] == "BOD" | mm1$Pop[i] == "CSI" |mm1$Pop[i] == "LH" |mm1$Pop[i] == "DQ"){ mm1$nat_env_factor[i] = "E_2"
  }else{ mm1$nat_env_factor[i] = "E_1"}
}

ma = data.frame("data_type" = rep("raw", nrow(mm1)),
                "gen_factor" = mm1$gen_factor,
                "exp_env_factor" = factor(mm1$exp_env_factor), 
                "nat_env_factor" = factor(mm1$nat_env_factor), # E_2 = coastal; E_1 = inland
                "phen_data"= mm1$age)

# Make into means data for testing: 
ma_means = ma %>%
  group_by(gen_factor,exp_env_factor) %>%
  summarize("phen_data_mean" = mean(phen_data),
            "phen_sd" = sd(phen_data),
            "sample_size" = n())
ma_means$phen_se = ma_means$phen_sd/(sqrt(sample_size))
ma_means$phen_data = ma_means$phen_data_mean
ma_means$nat_env_factor = mm1$nat_env_factor[match(ma_means$gen_factor,factor(mm1$gen_factor))]


## Geoff's Shell mass growth (mg) data (Fig. 3)
gt = read.csv("~/Documents/GitHub/CnGV/data/GeoffMeansData.csv")
gt$data_type = rep("means",nrow(gt))

## Berven 1979 -- Table sumthin.
berven = data.frame("Temp" = c(18,23,28,18,23,28),
                    "Population" = c("Montane","Montane","Montane","Lowland","Lowland","Lowland"),
                    "SizeatMM" = c(3.61,1.98,3.27,4.63,2.33,2.22),
                    "SizeatMM_se" = c(.14,0.05,0.13,0.1,0.11,0.1),
                    "gen_factor" = c("G_1","G_1","G_1","G_2","G_2","G_2"),
                    "nat_env_factor" = c("E_1","E_1","E_1","E_2","E_2","E_2"),
                    "exp_env_factor" = c("E_1","E_3","E_2","E_1","E_3","E_2"),
                    "TimetoMM" = c(347.5,106.1,140.4,419.5,136.6,99.6),
                    "TimetoMM_se" = c(14.9,4.4,15.7,20.2,10.1,5.2),
                    "data_type" = rep("means",6))
berven1 = berven %>%
  filter(exp_env_factor != "E_3") %>%
  droplevels()

size_df = data.frame("phen_data" = berven1$SizeatMM,"phen_se" = berven1$SizeatMM_se, berven1[,c(5:7,10)])
days_df = data.frame("phen_data" = berven1$TimetoMM,"phen_se" = berven1$TimetoMM_se, berven1[,c(5:7,10)])


# Molly and Geoff Datasets  
MollyTest = amarillo_armadillo(ma, 999, ma$data_type) # Molly's raw data
GeoffTest2002 = amarillo_armadillo(gt[gt$Paper == 2002,], 999, data_type = gt$data_type)
GeoffTest2000 = amarillo_armadillo(gt[gt$Paper == 2000,], 999, data_type = gt$data_type)


Berven_size = amarillo_armadillo(size_df,1000)
Berven_time = amarillo_armadillo(days_df,1000)


# Plots for paper
colors = c("E_1" = "#440154FF", "E_2" = "#29AF7FFF")
molly_labels = c("E_1" = "Inland populations", "E_2" = "Coastal populations")
geoff_labels1 = c("E_1" = "Sheltered populations", "E_2" = "Wave-exposed populations")
geoff_labels2 = c("E_1" = "Southern Population", "E_2" = "Northern Population")

# Standardize data for plots 
ma$phen_corrected =  (ma$phen_data - mean(ma$phen_data))/sd(ma$phen_data)
gt1 = gt[gt$Paper == 2002,]
gt1$phen_corrected =  (gt1$phen_data - mean(gt1$phen_data))/sd(gt1$phen_data)

label1 = MollyTest$Covariance.Estimate
label2 = paste0(":  C.I. = ", MollyTest$Covariance.Lower.CI," - ", MollyTest$Covariance.Upper.CI)
label3 = round(MollyTest$GxE.Estimate,3)
label4 = paste0(":  P-value = ",round(MollyTest$GxE.p.value,3))

(mollyplot = ggplot(ma, aes(x = exp_env_factor, y = phen_corrected, group = gen_factor, colour = nat_env_factor)) + 
  geom_point(position = position_dodge(width = 0.1))+geom_smooth(aes(fill=nat_env_factor),method = "glm")+
  xlab("")+ylab("Normalized Phenotype \n(age at metamorphosis)")+
  scale_x_discrete(breaks=c("E_1","E_2"),labels=c("Freshwater", "Saltwater"))+
  scale_colour_manual(values = colors,labels = molly_labels)+
  scale_fill_manual(values = colors,labels = molly_labels)+
  labs(col=" ",fill=" ")+
    theme_classic(base_size = 26, base_family = "Times")+
    theme(axis.text.x=element_text(colour="black"))+
    theme(axis.text.y=element_text(colour="black"))+
    theme(legend.position="none")+
  annotate("text", x = 1, y = 3.5, label = deparse(bquote("Cov"["GE"] ==~.(label1))),size=6, color = "black",hjust = 0,parse = T)+
  annotate("text", x = 1.35, y = 3.5, label = label2 ,size=6, color = "black",hjust = 0,parse = F)+
  annotate("text", x = 1, y = 3, label = deparse(bquote(bar(Delta)*""["GxE"] ==~.(label3))),size=6 , color = "black", hjust = 0,parse = T)+
  annotate("text", x = 1.3, y = 3, label = label4 ,size=6, color = "black",hjust = 0,parse = F))
  #Save as PDF A4

label5 = GeoffTest2002$Covariance.Estimate
label6 = paste0(":  C.I. = ", GeoffTest2002$Covariance.Lower.CI," - ", GeoffTest2002$Covariance.Upper.CI)
label7 = round(GeoffTest2002$GxE.Estimate,3)
label8 = paste0(":  P-value = ",round(GeoffTest2002$GxE.p.value,3))

(geoffplot2002 = ggplot(gt1, aes(x = exp_env_factor, y = phen_corrected, group = gen_factor, colour = nat_env_factor)) + 
  geom_point(size =5)+geom_line(size = 2)+
  geom_errorbar(aes(ymin = (phen_corrected-phen_se),ymax = (phen_corrected+phen_se)),size = 2, width = 0.1)+
  xlab("Environment")+ylab("Normalized Phenotype \n(Shell mass growth)")+
  scale_x_discrete(breaks=c("E_1","E_2"),labels=c("Low Flow", "High Flow"))+
  scale_colour_manual(values = colors,labels = geoff_labels1)+
  scale_fill_manual(values = colors,labels = geoff_labels1)+
  labs(col=" ")+
  theme_classic(base_size = 26, base_family = "Times")+
    theme(axis.text.x=element_text(colour="black"))+
    theme(axis.text.y=element_text(colour="black"))+
    theme(legend.position="none")+
  annotate("text", x = 1, y = 2.5, label = deparse(bquote("Cov"["GE"] ==~.(label5))),size=6, color = "black",hjust = 0,parse = T)+
  annotate("text", x = 1.35, y = 2.5, label = label6 ,size=6, color = "black",hjust = 0,parse = F)+
  annotate("text", x = 1, y = 2, label = deparse(bquote(bar(Delta)*""["GxE"] ==~.(label7))),size=6 , color = "black", hjust = 0,parse = T)+
  annotate("text", x = 1.3, y = 2, label = label8 ,size=6, color = "black",hjust = 0,parse = F))


geoffplot2000 = ggplot(gt[gt$Paper == 2000,], aes(x = exp_env_factor, y = phen_data, group = gen_factor, colour = nat_env_factor)) + 
  geom_point()+geom_line()+
  geom_errorbar(aes(ymin = (phen_data-phen_se),ymax = (phen_data+phen_se)),width = 0.1)+
  xlab("Environment")+ylab("Phenotype (Shell thickness growth (mm))")+
  scale_x_discrete(breaks=c("E_1","E_2"),labels=c("South", "North"))+
  scale_colour_manual(values = colors,labels = geoff_labels2)+
  scale_fill_manual(values = colors,labels = geoff_labels2)+
  labs(col=" ")+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(axis.text.x = element_text(size=14,colour = "black"),
        axis.title.x = element_text(size=16,face="bold")) +
  theme(axis.text.y = element_text(size=14,colour = "black"),
        axis.title.y = element_text(size=16,face="bold")) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2)) +
  annotate("text", x = "E_1", y = 0.75, 
           label = paste0("Cov = ", GeoffTest2000$Covariance.Estimate,", p = ", round(GeoffTest2000$Covariance.p.value,3)), size = 6, hjust = 0)+
  annotate("text", x = "E_1", y = 0.7, 
           label = paste0("GxE = ",round(GeoffTest2000$GxE.Estimate,2),", p = ",round(GeoffTest2000$GxE.p.value,3)), size = 6, hjust = 0)
geoffplot2000      

# Test data on oysters
setwd("~/Desktop/oyster_data/")

oystermass = read.csv("mass.csv")
mass.test = amarillo_armadillo(oystermass, 99, "raw") # Sig CnGV, non-sig GxE

oysterlength = read.csv("length.csv")
length.test = amarillo_armadillo(oysterlength, 99, "raw") # No sig CovGE, no sig GxE

oysterwidth = read.csv("width.csv") 
width.test = amarillo_armadillo(oysterwidth, 99, "raw") # No sig CovGE, no sig GxE

oystershellarea = read.csv("tissueshellarea.csv") 
shellarea.test = amarillo_armadillo(oystershellarea, 99, "raw")

oystershelllength = read.csv("tissueshelllength.csv")
shelllength.test = amarillo_armadillo(oystershelllength, 99, "raw")







