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


## Geoff's Shell mass growth (mg) data (Fig. 3)
gt = read.csv("~/Documents/GitHub/CnGV/data/GeoffMeansData.csv")
gt$data_type = rep("means",nrow(gt))

## Berven 1979
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






