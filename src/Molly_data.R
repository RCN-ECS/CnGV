# Molly's age at metamorphosis data
mm = read.csv("~/Desktop/Work/DataSets/Tadpole Plasticity_2017/mortality2017.csv")

# Massage data to just what I need (Genotype/Environment/PHenotype)
mm = mm[-which(is.na(mm$Jul_metamorph)),] # Include only those that metamorphosed
mm$age = mm$Jul_metamorph-mm$Jul_hatch # Calculate age at MM
mm = mm %>%
  filter(tad %in% c(0,6)) %>% # Exclude 8ppt tadpole treatment bc only 2 from coastal mm'd
  
# Rename variables for analysis
mm$gen_factor = paste0("G_",as.numeric(mm$loc))
mm$exp_env_factor = paste0("E_",as.numeric(as.factor(mm$tad))) # G_1 native to E_2; G_2 native to E_1
mm$nat_env_factor_corrected = NULL
for(i in 1:nrow(mm)){
  if(mm$gen_factor[i] == "G_1"){ mm$nat_env_factor_corrected[i] = "E_2"
    }else{ mm$nat_env_factor_corrected[i] = "E_1"}
}

input_df = data.frame("index" = rep(1, nrow(mm)),
                      "gen_factor" = mm$gen_factor,
                      "exp_env_factor" = mm$exp_env_factor,
                      "nat_env_factor" = mm$nat_env_factor_corrected,
                      "phen_data"= mm$age)
