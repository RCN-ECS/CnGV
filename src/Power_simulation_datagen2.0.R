#######################################
## Parameter Generation -- NO Groups ##
#######################################


# Starting list of parameters
param_list <- list( 
  reps = c(10), 
  delta_env = NULL, 
  delta_gen = c(-1,0.0,1),
  sample_size = c(4,8,16), 
  n_pop = c(2,4,8,16),
  env_scenario = c(1,2),  # 1 = Recip. Transplant ; 2 = Common Garden
  std_dev= c(.5, 1),
  interaction = NULL) 


parameter_generation <- function(param_list){
  
  ############################################################
  ################# Full Reciprocal Transplant ############### 
  ############################################################
  
  params =  expand.grid("n_pop"=param_list$n_pop,
                        "sample_size" = param_list$sample_size,
                        "std_dev" = param_list$std_dev)
  
  # Full Reciprocal Transplant
  params$n_env = params$n_pop
  
  param_table2 = data.frame()
  
  for(z in 1:param_list$reps){
    
    cat(z,"1st","\n")
    param_table1 = data.frame()
    
    for(i in 1:nrow(params)){
     
       # Slopes
      delta_env = runif(n = 50, min = 0, max = 1)
      
      # Intercepts
      delta_gen = runif(n = 50, min = -1, max = 1) 
      
      #if(params$n_pop[i] == 2){
      #  low_int = runif(50, min = 0, max = 2)
      #  interaction = low_int
     # }else{
        low_int = runif(n = 35, min = 0, max = 1.99)
        high_int = runif(n = 15, min = 2, max = params$n_pop[i])
        pot1 = c(low_int, high_int)
        interaction = sample(pot1, size = 50, replace = FALSE)
      #}
      inter_df1 = data.frame("delta_env" = delta_env, "delta_gen" = delta_gen, "interaction" = interaction)
      inter_df = merge(params[i,],inter_df1)
      
      param_table1 = rbind(param_table1, inter_df)
    }
    
    # Add Replicate
    param_table1$replicate = rep(z, nrow(param_table1))
    
    # Add scenario label
    param_table1$env_scenario = rep(1, nrow(param_table1))
    
    # Add rows to ensure false positives and negatives
    fpr1_gxe <- data.frame("n_pop" = sample(rep(c(2,4,8,16),each = 25), 100, replace = FALSE),
                           "sample_size" = sample(rep(c(4,8,16),each = 34), 100, replace = FALSE), # remove the 2, each = 34
                           "std_dev" = rep(c(0.5,1.0), each = 50),
                           "n_env" = NA, 
                           "delta_env" = runif(n = 100, min = 0, max = 1),
                           "delta_gen" = rep(0, 100),
                           "interaction" = rep(0, 100),
                           "replicate" = rep(z, 100),
                           "env_scenario" = rep(1, 100))
    fpr1_gxe$n_env = fpr1_gxe$n_pop
    fpr1_cov <- data.frame("n_pop" = sample(rep(c(2,4,8,16),each = 25), 100, replace = FALSE),
                           "sample_size" = sample(rep(c(4,8,16),each = 34), 100, replace = FALSE), # remove the 2, each = 34
                           "std_dev" = rep(c(0.5,1.0), each = 50),
                           "n_env" = NA, 
                           "delta_env" = rep(0, 100),
                           "delta_gen" = runif(n = 100, min = -1, max = 1),
                           "interaction" = rep(0, 100),
                           "replicate" = rep(z, 100),
                           "env_scenario" = rep(1, 100))
    fpr1_cov$n_env = fpr1_cov$n_pop
    
    # Bind it all up
    param_table2 = rbind(param_table2, param_table1,fpr1_cov,fpr1_gxe)
  }
  
  ############################################################
  ################## Paired Common Garden #################### 
  ############################################################
  
  params2 =  expand.grid("n_pop"=2*param_list$n_pop[-4],
                         "sample_size" = param_list$sample_size,
                         "std_dev" = param_list$std_dev)
  params2$n_env = 2
  
  param_table4 = data.frame()
  
  for(x in 1:param_list$reps){
    
    cat(x,"2nd","\n")
    
    param_table3 = data.frame()
    
    for(i in 1:nrow(params2)){
      
      # Slopes
      delta_env1 = rnorm(2,mean = 0.2, sd = 0.125)
      delta_env2 = rnorm(2,mean = 0.4, sd = 0.125)
      delta_env3 = rnorm(2,mean = 0.6, sd = 0.125)
      delta_env4 = rnorm(2,mean = 0.8, sd = 0.125)
      delta_env5 = c(delta_env1, delta_env2, delta_env3, delta_env4)
      
      # Intercepts
      delta_gen1 = rnorm(2,mean = -0.5, sd = 0.5)
      delta_gen2 = rnorm(2,mean = 0, sd = 0.5)
      delta_gen3 = rnorm(2,mean = 0.5, sd = 0.5)
      delta_gen5 = c(delta_gen1, delta_gen2, delta_gen3)

      # Interaction term
      int = c(0.5, params2$n_pop[i])
      
      dat = expand.grid(delta_env5, delta_gen5,int)
   
      inter_df1 = data.frame("delta_env" = dat[,1], "delta_gen" = dat[,2], "interaction" = dat[,3])
      inter_df2 = merge(params2[i,],inter_df1)
      
      param_table3 = rbind(param_table3, inter_df2)
    }
    
    # Add replicate
    param_table3$replicate = rep(x, nrow(param_table3))
    
    # Label scenario
    param_table3$env_scenario = rep(2, nrow(param_table3))
    
    # Add rows to ensure false positives and negatives
    fpr2_gxe <- data.frame("n_pop" = sample(rep(c(2,4,8,16),each = 25), 100, replace = FALSE),
                           "sample_size" = sample(rep(c(4,8,16),each = 34), 100, replace = FALSE), #remove the 2, each = 34
                           "std_dev" = rep(c(0.5,1.0), each = 50),
                           "n_env" = rep(2, times = 100), 
                           "delta_env" = runif(n = 100, min = 0, max = 1),
                           "delta_gen" = rep(0, 100) , # 100 at ZERO
                           "interaction" = rep(0, 100), # 100 at ZERO
                           "replicate" = rep(x, 100),
                           "env_scenario" = rep(2, 100))
    fpr2_cov <- data.frame("n_pop" = sample(rep(c(2,4,8,16),each = 25), 100, replace = FALSE),
                           "sample_size" = sample(rep(c(4,8,16),each = 34), 100, replace = FALSE), #remove the 2, each = 34
                           "std_dev" = rep(c(0.5,1.0), each = 50),
                           "n_env" = rep(2, times = 100), 
                           "delta_env" = rep(0, 100), # 100 at ZERO
                           "delta_gen" = runif(n = 100, min = -1, max = 1),
                           "interaction" = rep(0, 100),
                           "replicate" = rep(x, 100),
                           "env_scenario" = rep(2, 100))
                       
    # Bind it all up! 
    param_table4 = rbind(param_table4, param_table3,fpr2_gxe,fpr2_cov)
  }
  
  # All Replicates together
  param_table = data.frame()
  param_table = rbind(param_table2,param_table4)
  
  # Set.seed for each sim
  param_table$seed = round(runif(nrow(param_table), min = 1, max = 100000000))
  
  # Assign Rows to dataframe
  row = seq(1:nrow(param_table))
  param_table = data.frame("row" = row, param_table)
  param_table$total_samples <- param_table$n_env * param_table$n_pop * param_table$sample_size
  
  return(param_table)

  }

df = parameter_generation(param_list)
str(df)
df2 = filter(df, total_samples < 500) 
dim(df2)
write.csv(df2, "~/Desktop/df_rerun.csv")

setwd("~/Documents/GitHub/CnGV/CnGV/results/Sim_12.1.20/")
df=read.csv("df.csv")
dim(df)
args = filter(df, row == 25)

# Pull out test scenarios

df1 %>% filter(env_scenario == 2) #%>% filter(n_pop == 8) %>% filter(std_dev ==1) %>% filter(errpop == 0)
args = df1[940,]

df_sim = read.csv("~/Desktop/df.csv")
args = filter(df_sim, row == 16858)
args = args[,-1]

# 128 total samples
ges.448 = filter(df_sim, row == 10088)
args = ges.448[,-1]
filter(dat_csv, row == 10088)

ges.882 = filter(df_sim, row == 146)
args = ges.882[,-1]
filter(dat_csv, row == 146)

# 256 total samples
ges.884 = filter(df_sim, row == 14712)
args = ges.884[,-1]
filter(dat_csv, row == 14712)

ges.4416 = filter(df_sim, row == 7084)
args = ges.4416[,-1]
filter(dat_csv, row == 7084)

