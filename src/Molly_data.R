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


