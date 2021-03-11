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
data.df1 <- read.csv("Meta_data_checked.csv")
output <- data.frame()

# Data wrangle
data.df1$Data_file_name <- as.character(data.df1$Study_ID_phenotype)
info.df$Data_file_name <- as.character(info.df$Study_ID_phenotype)
data.df <- data.df1 %>% 
  filter(is.na(Data_file_name) == FALSE) %>%
  filter(Data_file_name != "")


# Test data Singles
temp = data.df = filter(data.df, Data_file_name == "326_flowering_startdate_2017")
#temp = filter(data.df, data.df$Data_file_name == "45_growth_rate")
info.df = filter(info.df,  Data_file_name == "326_flowering_startdate_2017")

# Test data Group
#testies = c("45_ETS_activity","59_male pupal mass","126_fat_content","126_boldness_distance","221_WCHI_green_inner_cell_width")
#data.df = filter(data.df, data.df$Data_file_name %in% testies)

#Sanity plots
(untrans = ggplot(temp, aes(x = exp_env_factor, y = phen_data, group = paste(gen_factor,exp_env_factor), colour = gen_factor))+ 
  geom_boxplot() + geom_smooth()) #+ geom_errorbar(aes(ymin = phen_data-(phen_SD/sqrt(phen_n)),ymax= phen_data+(phen_SD/sqrt(phen_n)))))
#(trans = ggplot(temp, aes(x = exp_env_factor, y = avg_phen_corrected, group = gen_factor, colour = gen_factor))+ 
#    geom_point() + geom_line() + geom_errorbar(aes(ymin = avg_phen_corrected-error,ymax= avg_phen_corrected+error)))
#grid.arrange(untrans,trans)

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
  if(unique(is.na(temp$phen_SD)) == TRUE & unique(is.na(temp$phen_mean_SE)) == TRUE){data_type = "raw"}else{data_type = "means"}
  
  # Standardize, calculate total sample size, convert SD to SE for means data
  temp$error = NULL
  temp$tss = NULL
  if(data_type == "raw"){ 
    temp$group = paste(temp$gen_factor,temp$exp_env_factor,sep = "-")
    temp$phen_corrected = (temp$phen_data - mean(temp$phen_data, na.rm = TRUE))/sd(tapply(temp$phen_data, temp$group, mean, na.rm = TRUE))
    temp$tss = nrow(temp)
  }else{
    temp$avg_phen_corrected = (temp$phen_data - mean(temp$phen_data))/sd(temp$phen_data)
    temp$tss = sum(temp$phen_n)
    if(unique(is.na(temp$phen_SD))==TRUE){
      temp$error = temp$phen_mean_SE
    }else{
      temp$error =  temp$phen_SD/(sqrt(temp$phen_n))
      }
  }
  
  tryCatch(
    {temp2 = amarillo_armadillo(temp, n_boot, data_type)
  
    # Compile Results 
    temp_out <- data.frame("Index" = x, 
                           "Data_file_name" = unique(temp$Data_file_name),
                           "First.Author" =  unique(temp$First_Author),
                           "Data.type" = data_type,
                           "Phylum" = ref$phylum_division,
                           "Genus" = ref$genus,
                           "Species" = ref$species,
                           "Experimental.Design"= unique(temp$Design),
                           "Total_Sample_Size" = unique(temp$tss),
                           "Phenotype" = ref$phenotype,
                           "Phenotype.unit" = ref$phenotype_unit,
                           "Trait.class" = unique(temp$trait.class),
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
    write.csv(exit_df, paste0("/scratch/albecker/Power_analysis/META/fails/fail_",ref$Data_file_name,".csv")) 
    next }else{
    output = rbind(output,temp_out)
    }
}

write.csv(output, "Meta_analysis_results_Feb20.csv")

##################################  
##    Result Visualization      ##
##################################  

library(ggplot2)
library(viridis)

# Meta-analysis results
res <- read.csv("~/Desktop/Meta_analysis_results_Feb20.csv")

res$studyID = str_split_fixed(res$Data_file_name, "_", 2)[,1]
res$phenotype = str_split_fixed(res$Data_file_name, "_", 2)[,2]

# Summary Data
(N_authors = length(unique(res$First.Author)))
(N_phenotypes = length(unique(res$Data_file_name)))
(N_studies = length(unique(res$studyID)))
(range_SampleSize = range(res$Total_Sample_Size))
(Means_studies = res %>% 
  filter_at(vars(Data.type), all_vars(.=="means")) %>% 
  select(studyID) %>% 
  distinct)
(N_means = length(Means_studies[,1])/length(unique(res$studyID)))

### Plotting Specs ###

# Shape according to Phylum
phyShape <- c("Chordata"=15,"Cnidaria"=16,"Arthropoda"=17,"Dinoflagellata"=18,"Bryophyta"=0,"Tracheophyta"=1,"Mollusca"=2, "Coniferophyta" =4,"Annelida"=5)
cols = c("4" = "#1f77b4","2"="#ff7f0e", "3" = "#2ca02c", "1" = "#d62728")

# Sort and Label according to significance 
res$is.GxE.Sig = NULL
res$is.Cov.SigCI = NULL
res$is.Cov.SigPval = NULL

res = res %>%
  mutate(is.GxE.Sig = ifelse(GxE_Pvalue <= 0.05,"Yes", "No")) %>%
  mutate(is.Cov.SigCI = ifelse(Covariance_Pvalue <= 0.025, "Yes",
                             ifelse(Covariance_LCI < 0 & Covariance_UCI < 0, "Yes",
                                    ifelse(Covariance_LCI > 0 & Covariance_UCI > 0, "Yes", "No"))))

res = res %>%
  mutate(is.Cov.SigPval = ifelse(Covariance_Pvalue <= 0.025, "Yes", "No"))

res$colorCI = NULL
res = res %>%
  mutate(colorCI = ifelse(is.GxE.Sig == "Yes" & is.Cov.SigCI == "Yes", "3",
                            ifelse(is.GxE.Sig == "No" & is.Cov.SigCI == "Yes", "1",
                                   ifelse(is.GxE.Sig == "Yes" & is.Cov.SigCI == "No", "2", "4")))) 

res$colorPval = NULL
res = res %>%
  mutate(colorPval = ifelse(is.GxE.Sig == "Yes" & is.Cov.SigPval == "Yes", "3",
                          ifelse(is.GxE.Sig == "No" & is.Cov.SigPval == "Yes", "1",
                                 ifelse(is.GxE.Sig == "Yes" & is.Cov.SigPval == "No", "2", "4")))) 

## Spread of Significant CovGE 
ggplot(filter(res,colorPval == 1), aes(x = Covariance))+geom_histogram()+theme_classic()
ggplot(filter(res,colorCI == 1), aes(x = Covariance))+geom_histogram()+theme_classic()


## Study effects (most studies contained multiple phenotypes)
(cov_study = ggplot(res, aes(x = studyID, y = Covariance))+geom_point()+theme_linedraw())
(gxe_study = ggplot(res, aes(x = studyID, y = GxE_Estimate))+geom_point()+theme_linedraw())

## Covariance Pvalues vs. Total Sample Size
(cov_ss = ggplot(res, aes(x = Covariance, y = Total_Sample_Size,colour = Covariance_Pvalue))+
    geom_point(alpha = 0.75) + geom_hline(yintercept = 256,linetype = "dashed")+geom_vline(xintercept = 0)+
    labs(colour = "P-value")+
    scale_colour_viridis(option = "plasma")+
    theme_classic(base_size = 20, base_family = "Arial")+
    xlab(expression("Cov"["GE"]))+ ylab("Total Sample Size")+
    theme(axis.text = element_text(colour = "black"))+
    facet_wrap(~Phylum,scales = "free", ncol = 3))

## GxE Pvalues vs. Total Sample Size
(gxe_ss = ggplot(res, aes(x = GxE_Estimate, y = Total_Sample_Size, colour = GxE_Pvalue))+
    geom_point(alpha = 1) + geom_hline(yintercept = 128, linetype = "dashed")+ 
    labs(colour = "P-value")+
    scale_colour_viridis(option = "plasma")+
    theme_classic(base_size = 20, base_family = "Arial")+
    xlab(expression(""*bar(Delta)*""["GxE"]))+ ylab("Total Sample Size")+
    theme(axis.text = element_text(colour = "black"))+
    facet_grid(rows = vars(Phylum), cols = vars(Data.type),scales = "free"))

# Overall Patterns in Covariance/GxE
(CI_overall = ggplot(res, aes(x = Covariance, y = GxE_Estimate, shape = colorCI, fill = colorCI))+
    labs(fill = "Significance", shape = "Significance")+
    geom_vline(aes(xintercept = 0))+
    geom_point(size = 4, alpha = 0.75)+
    scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24),
                       labels = c("1" = expression("Cov"["GE"]*""), "2" = expression(""*bar(Delta)*""["GxE"]) , "3" = "Both ","4"= "Neither "))+
    scale_fill_manual(values = cols, labels = c("1" = expression("Cov"["GE"]*""), "2" = expression(""*bar(Delta)*""["GxE"]) , "3" = "Both ","4"= "Neither "))+
    xlab(expression("Cov"["GE"]))+
    ylab(expression(bar(Delta)*""["GxE"]))+
    theme_classic(base_family = "Times", base_size = 20) +
    ggtitle("B   Confidence Intervals")+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(legend.text.align = 0)+
    theme(legend.position = "none")+
    #guides(fill=guide_legend(override.aes=list(shape=21)))+
    theme(axis.text = element_text(colour = "black")))

(Pval_overall = ggplot(res, aes(x = Covariance, y = GxE_Estimate, shape = colorPval, fill = colorPval))+
    labs(fill = "Significance", shape = "Significance")+
    geom_vline(aes(xintercept = 0))+
    geom_point(size = 4, alpha = 0.75)+
    scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24),
                       labels = c("1" = expression("Cov"["GE"]*""), "2" = expression(""*bar(Delta)*""["GxE"]) , "3" = "Both ","4"= "Neither "))+
    scale_fill_manual(values = cols, labels = c("1" = expression("Cov"["GE"]*""), "2" = expression(""*bar(Delta)*""["GxE"]) , "3" = "Both ","4"= "Neither "))+
    xlab(expression("Cov"["GE"]))+
    ylab(expression(bar(Delta)*""["GxE"]))+
    theme_classic(base_family = "Times", base_size = 20) +
    ggtitle("A   P-values")+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(legend.text.align = 0)+
    theme(legend.position = "none")+
    theme(axis.text = element_text(colour = "black")))

grid.arrange(Pval_overall, CI_overall,ncol= 2)

# According to various classifications
ggplot(res, aes(x = Covariance, y = GxE_Estimate,shape = colorCI, fill = colorCI))+
  labs(fill = "Significance", shape = "Significance")+
  geom_vline(aes(xintercept = 0))+
  geom_point(size = 4, alpha = 0.75)+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24),
                     labels = c("1" = expression("Cov"["GE"]*""), "2" = expression(""*bar(Delta)*""["GxE"]) , "3" = "Both ","4"= "Neither "))+
  scale_fill_manual(values = cols, labels = c("1" = expression("Cov"["GE"]*""), "2" = expression(""*bar(Delta)*""["GxE"]) , "3" = "Both ","4"= "Neither "))+
  xlab(expression("Cov"["GE"]))+
  ylab(expression(bar(Delta)*""["GxE"]))+
  theme_classic(base_family = "Times", base_size = 20) +
  #ggtitle("Significance according to Confidence Interval")+
  theme(plot.title = element_text(size = 20, face = "bold"))+
  theme(legend.text.align = 0)+
  theme(legend.position = "none")+
  #guides(fill=guide_legend(override.aes=list(shape=21)))+
  theme(axis.text = element_text(colour = "black"))+
  facet_wrap(~factor(Trait.class),ncol=3) # Phylum, Trait.class, Environment.type, Experimental.Design


# According to Trait classification
ggplot(filter(res,is.na(colorPval)!=TRUE), aes(x = Covariance, y = GxE_Estimate, fill = colorPval))+
  labs(fill = "Significance")+
  geom_vline(aes(xintercept = 0))+
  geom_point(shape = 21, size = 4) + 
  scale_fill_viridis(option = "plasma", discrete = TRUE,#values = c("1" = "#20A486FF","2" = "#481F70FF", "3" = "#FDE725FF", "4" = "white"),
                     labels = c("1" = expression("Cov"["GE"]*" Significant"), "2" = expression(""*bar(Delta)*""["GxE"]*" Significant") , "3" = "Both Significant","4"="None Significant"))+
  xlab(expression("Cov"["GE"]))+
  ylab(expression(bar(Delta)*""["GxE"]))+
  theme(axis.text = element_text(colour = "black"))+
  theme_classic(base_family = "Times", base_size = 20)+ facet_wrap(~factor(Trait.class),ncol=3)

## Endo vs. EctoTherms
chords = filter(res, Phylum == "Chordata")
ggplot(filter(chords,is.na(colorPval)!=TRUE), aes(x = Covariance, y = GxE_Estimate, group = Data.type, fill = colorPval))+
  labs(fill = "Significance")+
  geom_vline(aes(xintercept = 0))+
  geom_point(shape = 21, size = 4) + 
  scale_fill_viridis(option = "plasma", discrete = TRUE,#values = c("1" = "#20A486FF","2" = "#481F70FF", "3" = "#FDE725FF", "4" = "white"),
                     labels = c("1" = expression("Cov"["GE"]*" Significant"), "2" = expression(""*bar(Delta)*""["GxE"]*" Significant") , "3" = "Both Significant","4"="None Significant"))+
  xlab(expression("Cov"["GE"]))+
  #geom_label_repel(aes(label = Total_Sample_Size), nudge_x = 0, na.rm = TRUE)+
  ylab(expression(bar(Delta)*""["GxE"]))+
  theme(legend.position="bottom")  +
  theme_linedraw(base_family = "Times", base_size = 20)+ 
  theme(axis.text=element_text(colour = "black"))+
  facet_wrap(~Phylum_Divided)

## Genus
ggplot(res, aes(x = Genus, y = Covariance, fill = colorPval)) + 
  geom_point(shape = 21, size = 2) +
  labs(fill = "Significance")+
  scale_fill_viridis(option = "plasma", discrete = TRUE,#values = c("1" = "#20A486FF","2" = "#481F70FF", "3" = "#FDE725FF", "4" = "white"),
                     labels = c("1" = expression("Cov"["GE"]*" Significant"), "2" = expression(""*bar(Delta)*""["GxE"]*" Significant") , "3" = "Both Significant","4"="None Significant"))+
  ylab(expression("Cov"["GE"]))+
  theme_linedraw(base_family = "Times", base_size = 20)+ 
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1))+
  theme(axis.text=element_text(colour = "black"))+
  facet_wrap(~factor(Phylum),scales = "free", ncol=3)

ggplot(res, aes(x = Genus, y = GxE_Estimate, fill = colorPval)) + 
  geom_point(shape = 21, size = 2) +
  labs(fill = "Significance")+
  scale_fill_viridis(option = "plasma", discrete = TRUE,#values = c("1" = "#20A486FF","2" = "#481F70FF", "3" = "#FDE725FF", "4" = "white"),
                     labels = c("1" = expression("Cov"["GE"]*" Significant"), "2" = expression(""*bar(Delta)*""["GxE"]*" Significant") , "3" = "Both Significant","4"="None Significant"))+
  ylab(expression("Cov"["GE"]))+
  theme_linedraw(base_family = "Times", base_size = 20)+ 
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1))+
  theme(axis.text=element_text(colour = "black"))+
  facet_wrap(~factor(Phylum),scales = "free", ncol=3)
  


