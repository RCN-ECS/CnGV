
data_generation <- function(input_df){
  
  dat_temp = data.frame()
  row=0 #for indexing
  
  # Starting Parameters
  type = input_df$type
  intercept_G1 = input_df$intercept_G1
  slope_G1 = input_df$slope_G1
  intercept_G2 = input_df$intercept_G2
  slope_G2 = input_df$slope_G2
  sd = input_df$sd
  sample_size = input_df$sample_size
  G1_env = input_df$G1_env
  G2_env = input_df$G2_env
  env_num = input_df$env_num
  env = input_df$env
  true_covGE = input_df$true_covGE
  is.GxE = input_df$is.GxE
  slope_diff = input_df$slope_diff
  
  
  for(a in 1:length(unique(env_num))){ 
    for(b in 1:length(unique(type))){
      
      env_num_temp = unique(env_num)[a]
      type_temp = type[b]
      
      # Set environmental parameters (x-axis)
      if(env_num_temp == 2){env = c(-2,2)
      }else if (env_num_temp == 5){env = c(-2,-1,0,1,2)
      }else{env = seq(from = -2.5, to = 2, by = 0.5)}
      
      # Match genotypes to phenotypes
      if(type_temp == "cngv"){((G1_env= 2) & (G2_env = -2)) 
      }else{ ((G1_env= -2) & (G2_env = 2))}
      
      for(d in 1:length(slope_G1)){
        for(e in 1:length(intercept_G2)){
          for(f in 1:length(slope_G2)){
            for(g in 1:length(sd)){
              for(h in 1:length(sample_size)){
                
                row=row+1 #index 
                
                # GxE logical
                if(slope_G1[d] == slope_G2[f]){is.GxE = "No_GxE"} else {is.GxE = "Yes_GxE"}
                
                # Difference in slope and intercept
                slope_diff = abs(slope_G2[f] - slope_G1[d])
                intercept_diff = (intercept_G2[e]-intercept_G1)
                
                # Generate data 
                phen1 = c(replicate(sample_size[h],(intercept_G1 + slope_G1[d] * env + (rnorm(env_num_temp, 0, sd[g])))))
                temp1 = data.frame("index" = row,
                                   "type" = rep(type_temp, length(env_num_temp*sample_size[h])),
                                   "is.GxE" = rep(is.GxE, length(env_num_temp*sample_size[h])),
                                   "env_num" = rep(env_num_temp, length(env_num_temp*sample_size[h])),
                                   "slope" = rep(slope_G1[d], length(env_num_temp*sample_size[h])),
                                   "slope_diff" = rep(slope_diff, length(env_num_temp*sample_size[h])),
                                   "intercept" = rep(intercept_G1, length(env_num_temp*sample_size[h])),
                                   "intercept_diff" = rep(intercept_diff, length(env_num_temp*sample_size[h])),
                                   "sample_size" = rep(sample_size[h], length(env_num_temp*sample_size[h])),
                                   "stdev" = rep(sd[g], length(env_num_temp*sample_size[h])),
                                   "native_env" = rep(G1_env,length(env_num_temp*sample_size[h])),
                                   "gen" = rep("G1", length(env_num_temp*sample_size[h])),
                                   "env" = rep(env, length(env_num_temp*sample_size[h])),
                                   "phen" = phen1)
                
                phen2 = c(replicate(sample_size[h], (intercept_G2[e] + slope_G2[f] * env + (rnorm(env_num_temp, 0, sd[g])))))
                temp2 = data.frame("index" = row,
                                   "type" = rep(type_temp, length(env_num_temp*sample_size[h])),
                                   "is.GxE" = rep(is.GxE, length(env_num_temp*sample_size[h])),
                                   "env_num" = rep(env_num_temp, length(env_num_temp*sample_size[h])),
                                   "slope" = rep(slope_G2[f], length(env_num_temp*sample_size[h])),
                                   "slope_diff" = rep(slope_diff, length(env_num_temp*sample_size[h])),
                                   "intercept" = rep(intercept_G2[e], length(env_num_temp*sample_size[h])),
                                   "intercept_diff" = rep(intercept_diff, length(env_num_temp*sample_size[h])),
                                   "sample_size" = rep(sample_size[h], length(env_num_temp*sample_size[h])),
                                   "stdev" = rep(sd[g], length(env_num_temp*sample_size[h])),
                                   "native_env" = rep(G2_env,length(env_num_temp*sample_size[h])),
                                   "gen" = rep("G2", length(env_num_temp*sample_size[h])),
                                   "env" = rep(env, length(env_num_temp*sample_size[h])),
                                   "phen" = phen2)
                dat_temp. = rbind(temp1, temp2)
                dat_temp = rbind(dat_temp, dat_temp.)
              }
            }
          }
        }
      }
    }
  }
  return(dat_temp)
}