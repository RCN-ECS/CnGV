sim_means_se <- function(data_generation_output){
  
  # Categorical means and SE
  if(data_generation_output$data_type == "categorical"){
    
    bigsumdat = data.frame()
    for(j in 1:length(unique(data_generation_output$gen))){
      for (k in 1:length(unique(data_generation_output$env))){
        gentemp = unique(data_generation_output$gen)[j]
        envtemp = unique(data_generation_output$env)[k]
        dat. = filter(data_generation_output, gen == gentemp)
        dat = filter(dat., env == envtemp)
        average = mean(dat$phen)
        stddev = sd(dat$phen)
        bigsumdat. = data.frame("phen_data" = average,
                                "std_dev" = stddev,
                                "gen_factor" = gentemp,
                                "exp_env_cat" = envtemp,
                                "phen_n" = length(dat$phen))
        bigsumdat = rbind(bigsumdat, bigsumdat.)
      }
    }
    bigsumdat$phen_mean_SE = bigsumdat$std_dev/sqrt(bigsumdat$phen_n)
    bigsumdat$nat_env_factor = as.factor(c("E1","E1","E2","E2"))
    bigsumdat$exp_env_factor = as.factor(c("E1","E2","E1","E2"))
    
  }else{
    
    ## Continuous means and SE
    bigsumdat = data.frame()
    for(j in 1:length(unique(data_generation_output$gen))){
      for (k in 1:length(unique(data_generation_output$env))){
        gentemp = unique(data_generation_output$gen)[j]
        envtemp = unique(data_generation_output$env)[k]
        dat. = filter(data_generation_output, gen==gentemp)
        dat = filter(dat., env == envtemp)
        average = mean(dat$phen)
        stddev = sd(dat$phen)
        bigsumdat. = data.frame("phen_data" = average,
                                "std_dev" = stddev,
                                "gen_factor" = gentemp,
                                "exp_env_cont" = envtemp,
                                "phen_n" = length(dat$phen))
        bigsumdat = rbind(bigsumdat, bigsumdat.)
      }
    }
    bigsumdat$phen_mean_SE = bigsumdat$std_dev/sqrt(bigsumdat$phen_n)
    bigsumdat$nat_env_cont = c(-2,-2,2,2)
  }
  return(bigsumdat)
}
