sim_means_se <- function(data_generation_output){
  
  # Output
  bigsumdat = data.frame()
  
  # Multiple indexes
  for(z in 1:length(unique(data_generation_output$index))){
    
    subdat = filter(data_generation_output, index == (unique(data_generation_output$index)[z]))
    bigsumdat. = data.frame()
    
    for(j in 1:length(unique(subdat$gen))){
      for (k in 1:length(unique(subdat$exp_env_factor))){
        dat. = filter(subdat, gen_factor == unique(subdat$gen_factor)[j])
        dat = filter(dat., exp_env_factor == unique(subdat$exp_env_factor)[k])
        average = mean(dat$phen_data)
        stddev = sd(dat$phen_data)
        bigsumdat.. = data.frame("index" = unique(dat$index),
                                 "phen_data" = average,
                                 "std_dev" = stddev,
                                 "gen_factor" = unique(subdat$gen_factor)[j],
                                 "exp_env_factor" = unique(subdat$exp_env_factor)[k],
                                 "phen_n" = length(dat$phen_data))
        bigsumdat. = rbind(bigsumdat., bigsumdat..)
      }
    }
    
    bigsumdat = rbind(bigsumdat, bigsumdat.)
  }
    bigsumdat$phen_mean_SE = bigsumdat$std_dev/sqrt(bigsumdat$phen_n)
    #bigsumdat.$nat_env_factor = as.factor(c("E1","E1","E2","E2"))
    #bigsumdat.$exp_env_factor = as.factor(c("E1","E2","E1","E2"))
    return(bigsumdat)
  }



