
data_generation <- function(input_df){
  
  dat_temp = data.frame()
  row = 0 #for indexing
  
  # Starting Parameters
  cat_cont = input_df$cat_cont
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
  
  for(a in 1:length(unique(cat_cont))){ 
    for(b in 1:length(intercept_G1)){
      for(d in 1:length(slope_G1)){
        for(e in 1:length(intercept_G2)){
          for(f in 1:length(slope_G2)){
            for(g in 1:length(sd)){
              for(h in 1:length(sample_size)){
                
                row=row+1 #index 
                
                env_num_temp = 2
                data_type_temp = unique(cat_cont)[a]

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
                
                # GxE logical
                if(slope_G1[d] == slope_G2[f]){is.GxE = "No_GxE"} else {is.GxE = "Yes_GxE"}
                
                # Difference in slope and intercept
                slope_diff = abs(slope_G2[f] - slope_G1[d])
                intercept_diff = abs(intercept_G2[e]-intercept_G1[b])
                
                # Generate data 
                phen1 = c(replicate(sample_size[h],(intercept_G1[b] + slope_G1[d] * env_cont + (rnorm(env_num_temp, 0, sd[g])))))
                temp1 = data.frame("source" = rep("sim", (env_num_temp*sample_size[h])),
                                   "index" = rep(row, (env_num_temp*sample_size[h])),
                                   "cat_cont" = rep(cat_cont, (env_num_temp*sample_size[h])),
                                   "is.GxE" = rep(is.GxE, (env_num_temp*sample_size[h])),
                                   "env_num" = rep(env_num_temp, (env_num_temp*sample_size[h])),
                                   "slope" = rep(slope_G1[d], (env_num_temp*sample_size[h])),
                                   "slope_diff" = rep(slope_diff, (env_num_temp*sample_size[h])),
                                   "intercept" = rep(intercept_G1[b], (env_num_temp*sample_size[h])),
                                   "intercept_diff" = rep(intercept_diff, (env_num_temp*sample_size[h])),
                                   "phen_n" = rep(sample_size[h], (env_num_temp*sample_size[h])),
                                   "stdev" = rep(sd[g], (env_num_temp*sample_size[h])),
                                   "nat_env_factor" = rep(G1_env,(env_num_temp*sample_size[h])),
                                   "gen_factor" = rep("G1", (env_num_temp*sample_size[h])),
                                   "exp_env_factor" = rep(env, (env_num_temp*sample_size[h])),
                                   "phen_data" = phen1)
                
                phen2 = c(replicate(sample_size[h], (intercept_G2[e] + slope_G2[f] * env_cont + (rnorm(env_num_temp, 0, sd[g])))))
                temp2 = data.frame("source" = rep("sim", (env_num_temp*sample_size[h])),
                                   "index" = rep(row, (env_num_temp*sample_size[h])),
                                   "cat_cont" = rep(cat_cont, (env_num_temp*sample_size[h])),
                                   "is.GxE" = rep(is.GxE, (env_num_temp*sample_size[h])),
                                   "env_num" = rep(env_num_temp, (env_num_temp*sample_size[h])),
                                   "slope" = rep(slope_G2[f], (env_num_temp*sample_size[h])),
                                   "slope_diff" = rep(slope_diff, (env_num_temp*sample_size[h])),
                                   "intercept" = rep(intercept_G2[e], (env_num_temp*sample_size[h])),
                                   "intercept_diff" = rep(intercept_diff, (env_num_temp*sample_size[h])),
                                   "phen_n" = rep(sample_size[h], (env_num_temp*sample_size[h])),
                                   "stdev" = rep(sd[g], (env_num_temp*sample_size[h])),
                                   "nat_env_factor" = rep(G2_env,(env_num_temp*sample_size[h])),
                                   "gen_factor" = rep("G2", (env_num_temp*sample_size[h])),
                                   "exp_env_factor" = rep(env, (env_num_temp*sample_size[h])),
                                   "phen_data" = phen2)
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
