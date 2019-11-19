a=b=c=d=e=f=g=h=i=j=k=l=m=1


data_generation <- function(input_df){
  
  dat_temp = data.frame()
  row = 0 #for indexing
  
  # Starting Parameters
  data_type = input_df$data_type
  intercept_G1 = input_df$intercept_G1
  slope_G1 = input_df$slope_G1
  intercept_G2 = input_df$intercept_G2
  slope_G2 = input_df$slope_G2
  sd = input_df$sd
  sample_size = input_df$sample_size
  type = NULL
  G1_env = NULL
  G2_env = NULL
  env_num = NULL
  env = NULL
  true_covGE = NULL
  is.GxE = NULL
  slope_diff = NULL
  env_cont = NULL
  
  for(a in 1:length(unique(data_type))){ 
    for(b in 1:length(intercept_G1)){
      for(d in 1:length(slope_G1)){
        for(e in 1:length(intercept_G2)){
          for(f in 1:length(slope_G2)){
            for(g in 1:length(sd)){
              for(h in 1:length(sample_size)){
                
                row=row+1 #index 
                
                env_num_temp = 2
                data_type_temp = unique(data_type)[a]
                #type_temp = unique(type)[b]
                
                # Assign genotypes to native environments
                if(data_type_temp=="continuous"){
                  G1_env = -2
                  G2_env = 2
                  env = c(-2,2)
                  env_cont = c(-2,2)
                }else{
                  G1_env = "E1"
                  G2_env = "E2"
                  env = c("E1","E2")
                  env_cont = c(-2,2)
                  }
              
                # Assign covariance type
                if((intercept_G2[e] == -intercept_G1[b]) & (slope_G1[d] == -slope_G2[f])){
                  type = "pure_GxE"
                }else if((slope_G1[d]<0) & (intercept_G1[b]<intercept_G2[e])){
                  type = "cngv"
                }else if((slope_G1[d]>0) & (intercept_G1[b]<intercept_G2[e])){
                  type = "cogv"
                }else if((slope_G1[d]>0) & (intercept_G1[b]>intercept_G2[e])){
                  type = "cngv"
                }else if((slope_G1[d]<0) & (intercept_G1[b]>intercept_G2[e])){
                  type = "cogv"
                }else{
                  type = "check"
                }
                
                # GxE logical
                if(slope_G1[d] == slope_G2[f]){is.GxE = "No_GxE"} else {is.GxE = "Yes_GxE"}
                
                # Difference in slope and intercept
                slope_diff = abs(slope_G2[f] - slope_G1[d])
                intercept_diff = (intercept_G2[e]-intercept_G1[b])
                
                # Generate data 
                phen1 = c(replicate(sample_size[h],(intercept_G1[b] + slope_G1[d] * env_cont + (rnorm(env_num_temp, 0, sd[g])))))
                temp1 = data.frame("index" = row,
                                   "data_type" = rep(data_type, length(env_num_temp*sample_size[h])),
                                   "type" = rep(type, length(env_num_temp*sample_size[h])),
                                   "is.GxE" = rep(is.GxE, length(env_num_temp*sample_size[h])),
                                   "env_num" = rep(env_num_temp, length(env_num_temp*sample_size[h])),
                                   "slope" = rep(slope_G1[d], length(env_num_temp*sample_size[h])),
                                   "slope_diff" = rep(slope_diff, length(env_num_temp*sample_size[h])),
                                   "intercept" = rep(intercept_G1[b], length(env_num_temp*sample_size[h])),
                                   "intercept_diff" = rep(intercept_diff, length(env_num_temp*sample_size[h])),
                                   "sample_size" = rep(sample_size[h], length(env_num_temp*sample_size[h])),
                                   "stdev" = rep(sd[g], length(env_num_temp*sample_size[h])),
                                   "native_env" = rep(G1_env,length(env_num_temp*sample_size[h])),
                                   "gen" = rep("G1", length(env_num_temp*sample_size[h])),
                                   "env" = rep(env, length(env_num_temp*sample_size[h])),
                                   "phen" = phen1)
                
                phen2 = c(replicate(sample_size[h], (intercept_G2[e] + slope_G2[f] * env_cont + (rnorm(env_num_temp, 0, sd[g])))))
                temp2 = data.frame("index" = row,
                                   "data_type" = rep(data_type, length(env_num_temp*sample_size[h])),
                                   "type" = rep(type, length(env_num_temp*sample_size[h])),
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
