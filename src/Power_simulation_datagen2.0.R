#######################################
## Parameter Generation -- NO Groups ##
#######################################


# Starting list of parameters
param_list <- list( 
  reps = c(10), 
  delta_env = NULL, 
  delta_gen = c(-1,0.0,1),
  sample_size = c(2,4,8,16), 
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
      
      if(params$n_pop[i] == 2){
        low_int = runif(50, min = 0, max = 2)
        interaction = low_int
      }else{
        low_int = runif(n = 35, min = 0, max = 1.99)
        high_int = runif(n = 15, min = 2, max = params$n_pop[i])
        pot1 = c(low_int, high_int)
        interaction = sample(pot1, size = 50, replace = FALSE)
      }
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
                           "sample_size" = sample(rep(c(2, 4,8,16),each = 25), 100, replace = FALSE), # remove the 2, each = 34
                           "std_dev" = rep(c(0.5,1.0), each = 50),
                           "n_env" = NA, 
                           "delta_env" = runif(n = 100, min = 0, max = 1),
                           "delta_gen" = rep(0, 100),
                           "interaction" = rep(0, 100),
                           "replicate" = rep(z, 100),
                           "env_scenario" = rep(1, 100))
    fpr1_gxe$n_env = fpr1_gxe$n_pop
    fpr1_cov <- data.frame("n_pop" = sample(rep(c(2,4,8,16),each = 25), 100, replace = FALSE),
                           "sample_size" = sample(rep(c(2,4,8,16),each = 25), 100, replace = FALSE), # remove the 2, each = 34
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
    fpr2_gxe <- data.frame("n_pop" = sample(rep(c(4,8,16),each = 34), 100, replace = FALSE),
                           "sample_size" = sample(rep(c(2,4,8,16),each = 25), 100, replace = FALSE), #remove the 2, each = 34
                           "std_dev" = rep(c(0.5,1.0), each = 50),
                           "n_env" = rep(2, times = 100), 
                           "delta_env" = runif(n = 100, min = 0, max = 1),
                           "delta_gen" = rep(0, 100) , # 100 at ZERO
                           "interaction" = rep(0, 100), # 100 at ZERO
                           "replicate" = rep(x, 100),
                           "env_scenario" = rep(2, 100))
    fpr2_cov <- data.frame("n_pop" = sample(rep(c(4,8,16),each = 34), 100, replace = FALSE),
                           "sample_size" = sample(rep(c(2,4,8,16),each = 25), 100, replace = FALSE), #remove the 2, each = 34
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
dim(df)
df1 = filter(df, total_samples <500) 
dim(df1)
df=read.csv("~/Desktop/df.csv")
dim(df)

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



##############################################
##                 OLD METHOD               ##
##############################################

parameter_generation <- function(param_list){
  
  # Initial Starting Point  
  params <- expand.grid(#"delta_env" = param_list$delta_env,
    "delta_gen" = param_list$delta_gen,
    "sample_size" = param_list$sample_size,
    "env_scenario" = param_list$env_scenario,
    "std_dev" = param_list$std_dev)
  
  # Establish range of populations and environments 
  param_temp = data.frame()
  for(i in 1:nrow(params)){
    cat(i,"i","\n")
    if(params$env_scenario[i] == 1){ # Then the number of pops = number of environments
      envpop = data.frame("n_env" = param_list$n_pop, "n_pop" = param_list$n_pop)
      interdf = merge(params[i,],envpop)
    } else { # If n_pop is spread across 2 environments 
      envpop = data.frame("n_env" = 2, "n_pop" = 2*param_list$n_pop[-4])
      interdf = merge(params[i,],envpop)
    } 
    param_temp = rbind(param_temp, interdf)
  }
  
  # Create delta_env
  param_temp2 = data.frame()
  for(i in 1:nrow(param_temp)){
    if(param_temp$n_pop[i]==2 | param_temp$n_pop[i]==4){delta_env = runif(param_temp$n_pop[i]+1,min = -1,max =1)
    inter_df2 = merge(param_temp[i,],delta_env)
    } else {delta_env = runif(5,min = -1,max =1)
    inter_df2 = merge(param_temp[i,],delta_env)
    }
    param_temp2 = rbind(inter_df2,param_temp2)
  }
  colnames(param_temp2)[7]<- "delta_env"
  
  # Generate Interaction terms 
  param_temp3 = data.frame()
  for(j in 1:nrow(param_temp2)){
    
    cat(j,"j","\n")
    low_interaction_term <- seq(from = 0, to = 0.99,length.out = 4)
    med_interaction_term <- seq(from = 1, to = 1.99,length.out = 3)
    
    if(param_temp2$n_pop[j] == 2){
      high_interaction_term = 2
    }else{
      high_interaction_term <- sample(seq(from = 2, to = param_temp2$n_pop[j]),size = 3, replace = FALSE) 
    }
    interaction_term <- c(low_interaction_term,med_interaction_term,high_interaction_term)
    interdf3 <- merge(param_temp2[j,],interaction_term)
    param_temp3 <- rbind(param_temp3, interdf3)
  }
  colnames(param_temp3)[8]<- "interaction"
  
  
  # Add replicates
  reps <- rep(c(1:param_list$reps), each = nrow(param_temp3))
  param_temp4 <- data.frame("replicate" = reps, param_temp3)
  
  # Set.seed for each sim
  param_temp4$seed = round(runif(nrow(param_temp4), min = 1, max = 100000000))
  
  # Assign Rows to dataframe
  row = seq(1:nrow(param_temp4))
  param_table = data.frame("row" = row, param_temp4)
  param_table$total_samples <- param_table$n_env * param_table$n_pop * param_table$sample_size
  
  return(param_table)
}


df = parameter_generation(param_list) 
dim(df)

df1 = df %>%
  filter(replicate %in% c(1:10)) %>%
  filter(total_samples < 500) 

df2 = df1[!(df1[,3]==0 & df1[,4]==0),]
dim(df2) #95,760
ggplot(df2, aes(x = delta_gen))+geom_histogram()


#write.csv(df2,"~/Desktop/df.csv",)
x = data.frame("rep" = c(1:1000),"dist" = rchisq(1000,df=0.1))
ggplot(x, aes(x = dist))+geom_histogram(binwidth = 0.1)

# Estimate time on cluster 
nsims = length(which(df2$replicate==1))
pause = 8
chunksize = 1000
(clusterTime = (nsims/chunksize)*8)
(Days = clusterTime/24)




####### Which rows failed? #######
forthemissing = data.frame(anti_join(df_missing, dat_csv1,by = "row")) # dat_csv from "Sim.Plot.Code.R"
forthemissing = as.numeric(forthemissing[[1]])
repeat_df = filter(df2, row %in% forthemissing)
write.csv(repeat_df, "~/Desktop/repeatdf.csv")

###################################################
## Parameter Generation when Divided Into Groups ##
###################################################


# Starting list of parameters
param_list <- list( 
  reps = c(100), 
  delta_env = c(0.0,0.25,0.5,.75,1.0), 
  delta_gen = c(-1,0.0,1),
  sample_size = c(2,4,8,16), 
  n_pop = c(2,4,8,16),
  env_scenario = c(1,2),  # 1 = npop = n_env; 2 = multiple pops across 2 envs
  std_dev= c(.5, 1), # Scale
  interaction = NULL) # Set as function of n_pop

## Table of parameters

parameter_generation <- function(param_list){
  
  params <- expand.grid("delta_env" = param_list$delta_env,
                        "delta_gen" = param_list$delta_gen,
                        "sample_size" = param_list$sample_size,
                        "env_scenario" = param_list$env_scenario,
                        "std_dev" = param_list$std_dev)
  
  # Establish range of populations and environments 
  param_temp = data.frame()
  for(i in 1:nrow(params)){
    cat(i,"i","\n")
    if(params$env_scenario[i] == 1){ # Then the number of pops = number of environments
      envpop = data.frame("n_env" = param_list$n_pop, "n_pop" = param_list$n_pop)
      interdf = merge(params[i,],envpop)
    } else { # If n_pop is spread across 2 environments 
      envpop = data.frame("n_env" = 2, "n_pop" = 2*param_list$n_pop[-4])
      interdf = merge(params[i,],envpop)
    } 
    param_temp = rbind(param_temp, interdf)
  }
  
  # Generate Interaction term 
  param_temp3 = data.frame()
  for(j in 1:nrow(param_temp)){
    cat(j,"j","\n")
    interaction_term <- seq(from = 0, to = param_temp$n_pop[j],length.out = param_temp$n_pop[j]) 
    interdf2 <- merge(param_temp[j,],interaction_term)
    param_temp3 <- rbind(param_temp3, interdf2)
  }
  colnames(param_temp3)[8]<- "interaction"
  
  # Add replicates
  reps <- rep(c(1:param_list$reps), each = nrow(param_temp3))
  param_temp4 <- data.frame("replicate" = reps, param_temp3)
  
  # Set.seed for each sim
  param_temp4$seed = round(runif(nrow(param_temp4), min = 1, max = 100000000))
  
  # Assign Rows to dataframe
  row = seq(1:nrow(param_temp4))
  param_table = data.frame("row" = row, param_temp4)
  param_table$total_samples <- param_table$n_env * param_table$n_pop * param_table$sample_size
  
  # Assign groups for Slurm
  param_table$group = rep(1:(nrow(param_table)/5), each = 5)
  
  return(param_table)
}


df = parameter_generation(param_list) 
dim(df)

df1 = df %>%
  filter(replicate %in% c(1:10)) %>%
  filter(total_samples < 512) 

df2 = df1[!(df1[,3]==0 & df1[,4]==0),]
dim(df2) #38,080

write.csv(df2,"~/Desktop/df.csv",)
cluster.time <- function(n_params, chunk.size, n_rows, avg_time, allowed.jobs, pause.time){
  chunk.submissions <- (n_params/n_rows)/chunk.size
  timepergroup = n_rows*avg_time
  
  total.time.hours<- (total.submissions/allowed.jobs)*avg_time
  #t#otal.time.hours <- total.time.minutes/60
  total.time.days <- total.time.hours/24
  time_project = data.frame("days" = total.time.days,"hours" = total.time.hours, "minutes" = total.time.minutes, "jobs" = n_params)
}
cluster.time(n_params = nrow(df1), chunk.size = 788, avg_time = 20, pause.time = 600, groupsize = 5, jobs_allowed = 72)

cluster.time.nochunks <- function(n_params, avg_time, allowed.jobs){
  total.submissions <- n_params
  total.time.minutes <- (total.submissions/allowed.jobs)*avg_time
  total.time.hours <- total.time.minutes/60
  total.time.days <- total.time.hours/24
  time_project = data.frame("days" = total.time.days,"hours" = total.time.hours, "minutes" = total.time.minutes, "jobs" = n_params)
}
(time.check = cluster.time.nochunks(n_params = nrow(df1), avg_time=21, allowed.jobs = 72)) 

tail(df1)
len = df[df$replicate == 1,]
(time.check = cluster.time(n_params = nrow(df1), chunk.size = 200, pause.time = 600)) # Njobs = 27000

check <- df %>%
  filter(replicate %in% c(1:10)) #%>%
#filter(n_pop == 2)



###############################
# Old Code - Unused as of 7/1/2020
###############################


ring <- function(param_table, n_boot){
  start_time <- Sys.time()
  
  # Load packages
  library("emmeans")
  library("lme4")
  library("tidyverse")
  library("rlist")
  
  # Output dataframe
  output <- data.frame()
  
  for(i in 1:nrow(param_table)){
    
    # Counter
    cat(i, "\n")
    
    # For reproducibility
    # set.seed = 999
    model_df = NULL
    
    # Set Conditional Parameters
    n_environments = param_table$n_pop[i]
    
    # Approximate Cov(G,E)
    cov_GE_approx = param_table$delta_env[i] * param_table$delta_gen[i]
    
    # Dataframe foundations
    gen <- rep(1:param_table$n_pop[i], each = param_table$sample_size[i]*n_environments)
    env <- rep(1:n_environments, each = param_table$sample_size[i],times = n_environments) 
    
    # Random Noise
    noise <- rnorm(param_table$sample_size[i] * param_table$n_pop[i] * n_environments, 0, sd = param_table$std_dev[i]) 
    
    # Interaction Terms
    int <- rep(rnorm(param_table$n_pop[i] * n_environments, 0, sd = param_table$interaction[i]),each = param_table$sample_size[i]) # interaction term - one for each GE level
    
    # Create the model dataframe 
    model_df <- data.frame(gen, env, noise, int)
    model_df$gen_factor = factor(paste("G", model_df$gen, sep = "_"))
    model_df$exp_env_factor = factor(paste("E", model_df$env, sep = "_"))
    
    # Generate phenotype data using regression equation
    phen = param_table$delta_env[i] * model_df$env + param_table$delta_gen[i] * model_df$gen  + model_df$noise + model_df$int 
    model_df$phen = phen
    
    # Standardize data
    model_df$phen_corrected = (phen-mean(phen))/sd(phen)
    
    # Anova
    test_temp <- lm(phen_corrected ~ exp_env_factor * gen_factor, data = model_df)
    
    # Estimated Marginal Means
    emm_options(msg.interaction = FALSE)
    emm_E = emm_G = emm_GxE = NULL
    emm_E = as.data.frame(emmeans(test_temp,"exp_env_factor"))
    emm_G = as.data.frame(emmeans(test_temp, "gen_factor"))
    emm_GxE = as.data.frame(emmeans(test_temp, ~ exp_env_factor*gen_factor))
    
    # Gmeans
    G_matrix = data.frame()
    gtemp = gmean = tempdat = NULL
    for(h in 1:length(unique(emm_GxE$gen_factor))){
      gtemp <- filter(emm_GxE, gen_factor == unique(emm_GxE$gen_factor)[h])
      gmean <- sum(gtemp[,3])/length(unique(emm_GxE$gen_factor))
      tempdat = data.frame("G_means" = gmean,
                           "gen_factor" = unique(emm_GxE$gen_factor)[h])
      G_matrix = rbind(G_matrix,tempdat)
    }
    
    # Emeans
    E_matrix = data.frame()
    etemp = emean = tempdat. = NULL
    for(j in 1:length(unique(emm_GxE$exp_env_factor))){
      etemp <- filter(emm_GxE, exp_env_factor == unique(emm_GxE$exp_env_factor)[j])
      emean <- sum(etemp[,3])/n_environments
      tempdat. = data.frame("E_means" = emean,
                            "exp_env_factor" = unique(emm_GxE$exp_env_factor)[j])
      E_matrix = rbind(E_matrix,tempdat.)
    }
    
    # Covariance
    if(length(G_matrix[,1]) == 2){cov_est = cov(G_matrix$G_means,E_matrix$E_means)
    }else{cov_est = cor(G_matrix$G_means,E_matrix$E_means)}
    
    # Magnitude of GxE using EMMs
    GxE_emm <- abs(mean(model_df$phen_corrected) - # Overall mean
                     (emm_G$emmean[emm_G$gen_factor=="G_1"])- # G
                     (emm_E$emmean[emm_E$exp_env_factor=="E_1"])+ # E
                     (emm_GxE[1,3])) # GxE
    
    #############################
    ## True Covariance and GxE ##
    #############################
    
    # Generate phenotype data with no error using same regression equation
    no_err_phen = param_table$delta_env[i] * model_df$env + param_table$delta_gen[i] * model_df$gen + model_df$int
    
    # Standardize no error data
    model_df$no_err_phen_corrected = (no_err_phen-mean(no_err_phen))/sd(no_err_phen)
    
    # Anova with no error
    test_temp_noerror <- lm(no_err_phen_corrected ~ exp_env_factor * gen_factor, data = model_df)
    
    # Estimated Marginal Means with no error
    emm_options(msg.interaction = FALSE)
    emm_options(quietly = TRUE) # Will get perfect fit warning
    emm_E_ne = emm_G_ne = emm_GxE_ne = NULL
    emm_E_ne = as.data.frame(emmeans(test_temp_noerror,"exp_env_factor"))
    emm_G_ne = as.data.frame(emmeans(test_temp_noerror, "gen_factor"))
    emm_GxE_ne = as.data.frame(emmeans(test_temp_noerror, ~ exp_env_factor*gen_factor))
    
    # Gmeans with no error
    G_matrix_ne = data.frame()
    gtemp_ne = gmean_ne  = tempdat_ne = NULL
    for(y in 1:length(unique(emm_GxE_ne$gen_factor))){
      gtemp_ne <- filter(emm_GxE_ne, gen_factor == unique(emm_GxE_ne$gen_factor)[y])
      gmean_ne <- sum(gtemp_ne[,3])/length(unique(emm_GxE_ne$gen_factor))
      tempdat_ne = data.frame("G_means" = gmean_ne,
                              "gen_factor" = unique(emm_GxE_ne$gen_factor)[y])
      G_matrix_ne = rbind(G_matrix_ne,tempdat_ne)
    }
    
    # Emeans with no error
    E_matrix_ne = data.frame()
    etemp_ne = emean_ne = tempdat._ne = NULL
    for(f in 1:length(unique(emm_GxE_ne$exp_env_factor))){
      etemp_ne <- filter(emm_GxE_ne, exp_env_factor == unique(emm_GxE_ne$exp_env_factor)[f])
      emean_ne <- sum(etemp_ne[,3])/n_environments
      tempdat._ne = data.frame("E_means" = emean_ne,
                               "exp_env_factor" = unique(emm_GxE_ne$exp_env_factor)[f])
      E_matrix_ne = rbind(E_matrix_ne,tempdat._ne)
    }
    
    # Covariance
    if(length(G_matrix_ne[,1]) == 2){true_cov = cov(G_matrix_ne$G_means,E_matrix_ne$E_means)
    }else{true_cov = cor(G_matrix_ne$G_means,E_matrix_ne$E_means)}
    
    # Magnitude of GxE with no error 
    true_GxE <- abs(mean(model_df$no_err_phen_corrected) - # Overall mean
                      (emm_G_ne$emmean[emm_G_ne$gen_factor=="G_1"])- # G
                      (emm_E_ne$emmean[emm_E_ne$exp_env_factor=="E_1"])+ # E
                      (emm_GxE_ne[1,3])) # GxE
    
    ###############
    ## Bootstrap ##
    ###############
    
    # Output Dataframe
    boot_df = data.frame()
    
    # Resampling Loop
    for(a in 1:n_boot){
      new_phen <- NULL
      shuffle_dat <- data.frame()
      
      # Each genotype and environment
      for (l in 1:nlevels(model_df$gen_factor)){
        for (j in 1:nlevels(model_df$exp_env_factor)){
          cond_G <- filter(model_df, gen_factor == unique(model_df$gen_factor)[l])
          cond_E <- filter(cond_G, exp_env_factor == unique(model_df$exp_env_factor)[j])
          
          # Shuffle data 
          new_phen <- sample(cond_E$phen_corrected, size=nrow(cond_E), replace=TRUE)
          
          # Output    
          shuffle_dat_temp <- data.frame("gen_factor" = cond_E$gen_factor,
                                         "exp_env_factor" = cond_E$exp_env_factor,
                                         "phen_corrected" = new_phen)
          shuffle_dat <- rbind(shuffle_dat, shuffle_dat_temp)
        }
      }
      
      # Anova
      test_boot <- lm(phen_corrected ~ exp_env_factor * gen_factor, data = shuffle_dat)
      
      # Estimated Marginal Means
      emm_options(msg.interaction = FALSE)
      emm_E_boot = as.data.frame(emmeans(test_boot,"exp_env_factor"))
      emm_G_boot = as.data.frame(emmeans(test_boot, "gen_factor"))
      emm_GxE_boot = as.data.frame(emmeans(test_boot, ~ exp_env_factor*gen_factor))
      
      # Gmeans - Bootstrap
      G_matrix_boot = data.frame()
      for(p in 1:length(unique(emm_GxE_boot$gen_factor))){
        g_boot <- filter(emm_GxE_boot, gen_factor == unique(emm_GxE_boot$gen_factor)[p])
        g_mean_boot <- sum(g_boot[,3])/length(unique(emm_GxE_boot$gen_factor))
        tempdat = data.frame("G_means" = g_mean_boot,
                             "gen_factor" = unique(g_boot$gen_factor))
        G_matrix_boot = rbind(G_matrix_boot,tempdat)
      }
      
      # Emeans - Bootstrap
      E_matrix_boot = data.frame()
      for(q in 1:length(unique(emm_GxE_boot$exp_env_factor))){
        e_boot <- filter(emm_GxE_boot, exp_env_factor == unique(emm_GxE_boot$exp_env_factor)[q])
        e_mean_boot <- sum(e_boot[,3])/length(unique(emm_GxE_boot$exp_env_factor))
        tempdat. = data.frame("E_means" = e_mean_boot,
                              "exp_env_factor" = unique(e_boot$exp_env_factor))
        E_matrix_boot = rbind(E_matrix_boot,tempdat.)
      }
      
      # Covariance - Bootstrap
      Cov_matrix_boot = data.frame()
      Cov_matrix_boot <- cbind(G_matrix_boot,E_matrix_boot)
      if(param_table$n_pop[i]==2){cov_est_boot = cov(Cov_matrix_boot$G_means,Cov_matrix_boot$E_means)
      }else{cov_est_boot = cor(Cov_matrix_boot$G_means,Cov_matrix_boot$E_means)}
      
      # Magnitude of GxE - Bootstrap
      GxE_emm_boot <- abs(mean(shuffle_dat$phen_corrected) - # Overall mean
                            (emm_G_boot$emmean[emm_G_boot$gen_factor=="G_1"])- # G
                            (emm_E_boot$emmean[emm_E_boot$exp_env_factor=="E_1"])+ # E
                            (emm_GxE_boot[1,3])) # GxE
      
      boot_dat. <- data.frame("covariance" = cov_est_boot,
                              "GxE_mag" = GxE_emm_boot)
      
      boot_df <- rbind(boot_df,boot_dat.)
    }
    
    # Confidence Intervals
    GxE_CI = quantile(boot_df$GxE_mag, probs=c(0.025, 0.975), type=1) 
    GxE_avg = mean(boot_df$GxE_mag) 
    cov_CI = quantile(boot_df$covariance, probs=c(0.025, 0.975), type=1) 
    cov_avg = mean(boot_df$covariance) 
    
    #################
    ## Permutation ##
    #################
    
    # Output dataframe
    perm_df <- data.frame()
    
    # Shuffling loop
    for(b in 1:n_boot){
      
      # Shuffle data
      null_temp <- sample(model_df$phen_corrected, size=nrow(model_df), replace=FALSE)
      
      perm_dat <- data.frame("gen_factor" = model_df$gen_factor,
                             "exp_env_factor" = model_df$exp_env_factor,
                             "phen_corrected" = null_temp)
      
      # Anova
      test_perm <- lm(phen_corrected ~ exp_env_factor * gen_factor, data = perm_dat)
      
      # Estimated Marginal Means - Permutation
      emm_options(msg.interaction = FALSE)
      emm_E_perm = as.data.frame(emmeans(test_perm,"exp_env_factor"))
      emm_G_perm = as.data.frame(emmeans(test_perm, "gen_factor"))
      emm_GxE_perm = as.data.frame(emmeans(test_perm, ~ exp_env_factor*gen_factor))
      
      # Gmeans - Permutation
      G_matrix_perm = data.frame()
      for(r in 1:length(unique(emm_GxE_perm$gen_factor))){
        g_perm <- filter(emm_GxE_perm, gen_factor == unique(emm_GxE_perm$gen_factor)[r])
        g_mean_perm <- sum(g_perm[,3])/length(unique(emm_GxE_perm$gen_factor))
        tempdat = data.frame("G_means" = g_mean_perm,
                             "gen_factor" = unique(g_perm$gen_factor))
        G_matrix_perm = rbind(G_matrix_perm,tempdat)
      }
      
      # Emeans - Permutation
      E_matrix_perm = data.frame()
      for(s in 1:length(unique(emm_GxE_perm$exp_env_factor))){
        e_perm <- filter(emm_GxE_perm, exp_env_factor == unique(emm_GxE_perm$exp_env_factor)[s])
        e_mean_perm <- sum(e_perm[,3])/length(unique(emm_GxE_perm$exp_env_factor))
        tempdat. = data.frame("E_means" = e_mean_perm,
                              "exp_env_factor" = unique(e_perm$exp_env_factor))
        E_matrix_perm = rbind(E_matrix_perm,tempdat.)
      }
      
      # Covariance - Permutation
      Cov_matrix_perm = data.frame()
      Cov_matrix_perm <- cbind(G_matrix_perm,E_matrix_perm)
      if(param_table$n_pop[i]==2){cov_est_perm = cov(Cov_matrix_perm$G_means,Cov_matrix_perm$E_means)
      }else{cov_est_perm = cor(Cov_matrix_perm$G_means,Cov_matrix_perm$E_means)}
      
      # Magnitude of GxE - Permutation
      GxE_emm_perm <- abs(mean(perm_dat$phen_corrected) - # Overall mean
                            (emm_G_perm$emmean[emm_G_perm$gen_factor=="G_1"])- # G
                            (emm_E_perm$emmean[emm_E_perm$exp_env_factor=="E_1"])+ # E
                            (emm_GxE_perm[1,3])) # GxE
      
      perm_dat. <- data.frame("covariance" = cov_est_perm,
                              "GxE_mag" = GxE_emm_perm)
      perm_df <- rbind(perm_df,perm_dat.)
    }
    
    # Covariance P-value
    ptemp = (rank(c(cov_est,perm_df[,1]))[1])/(n_boot+1) 
    cov_pvalue = NULL
    if(ptemp < 0.5){cov_pvalue = ptemp}else{cov_pvalue = (1-ptemp)} # 2-tailed
    
    # GxE P-value
    ptemp1 = (rank(c(GxE_emm,perm_df[,2]))[1])/(n_boot+1) 
    GxE_pvalue = 1-ptemp1 # Right-tailed
    
    # Generate Outputs
    temp_out <- data.frame("row" = param_table$row[i],
                           "delta_env" = param_table$delta_env[i],
                           "delta_gen" = param_table$delta_gen[i],
                           "sample_size" = param_table$sample_size[i],
                           "n_pop" = param_table$n_pop[i],
                           "std_dev" = param_table$std_dev[i],
                           "interaction" = param_table$interaction[i],
                           "true_cov" = true_cov,
                           "cov_estimate" = cov_est,
                           "cov_lwrCI" = cov_CI[[1]],
                           "cov_uprCI" = cov_CI[[2]],
                           "cov_pvalue" = cov_pvalue,
                           "true_GxE" = true_GxE,
                           "GxE_estimate" = GxE_emm,
                           "GxE_lwrCI" = GxE_CI[[1]],
                           "GxE_uprCI" = GxE_CI[[2]],
                           "GxE_pvalue" = GxE_pvalue)
    output = rbind(output, temp_out)
    
    end_time <- Sys.time()
    
    time_difference = end_time - start_time
  }
  return(list(output,time_difference))
}

test2 = ring(df[df$row == 831,],100) # Parameter table, then number of bootstraps/perms  

test2$col = NULL
for(i in 1:nrow(test2)){
  if(test2$cov_pvalue[i] <= 0.05 & test2$GxE_pvalue[i] <= 0.05){test2$col[i] = "red" # Both significant
  }else if(test2$cov_pvalue[i] <= 0.05 & test2$GxE_pvalue[i] > 0.05){test2$col[i] = "darkgreen" # Cov significant
  }else if(test2$cov_pvalue[i] > 0.05 & test2$GxE_pvalue[i] <= 0.05){test2$col[i] = "dodgerblue" # GxE significant
  }else{test2$col[i] = "grey"} # None significant
}



# dat_csv$GxE_estimate = 1.855347e+00

# True GxE by Covariance
require(RColorBrewer)
ggplot(test2, aes(x = true_cov, y = true_GxE, group = factor(n_pop),colour=col)) + 
  geom_point() + theme_classic() + ylim(0,2)+ xlim(-1,1)+
  xlab("True Covariance") + ylab("True GxE") +
  scale_colour_identity()+
  #scale_color_brewer(palette="Set1",name = "Sample Size") +
  facet_grid(~interaction)

p1=ggplot(test2, aes(x = cov_estimate, y = GxE_estimate, group = factor(n_pop),colour = factor(sample_size))) + 
  geom_point() + theme_classic() + 
  xlab("Covariance Estimate") + ylab("GxE Estimate") +
  scale_color_brewer(palette="Set1",name = "Sample Size")+
  ggtitle("Simulation - runif coding")

p2=ggplot(nint, aes(x = cov_estimate, y = GxE_estimate, group = factor(n_pop),colour = factor(sample_size))) + 
  geom_point() + theme_classic() + 
  xlab("Covariance Estimate") + ylab("GxE Estimate") +
  scale_color_brewer(palette="Set1",name = "Sample Size")+
  ggtitle("Full Simulation - original coding")

require(lattice)
require(gridExtra)
# really want to know trade off between n_pop and sample size.
#write.table(test,"Power_data.txt")
#write.csv(test,"Power_data.csv")



##########
## Plot ##
##########

pdata <- read.csv("~/Desktop/power_output/compiled_simdata.csv")


pdata$tickmark <- NULL
for(i in 1:nrow(pdata)){
  if(pdata$cov_pvalue[i] > 0.05){pdata$tickmark[i]=0}else{pdata$tickmark[i]=1} # do 100 replicates for this power
}
# How many times is null false? Check true covariance and true GxE when noise = 0 for null expectation. (no sense in looking at power for noise = 0)

powerdata = data.frame()
for(i in 1:length(unique(pdata$index))){
  tempdata <- filter(pdata, index == unique(pdata$index)[i])
  num = sum(tempdata$tickmark)
  denom = length(tempdata$tickmark)
  power = num/denom
  powerdata. = data.frame("index" = unique(tempdata$index),
                          "power" = power)
  powerdata = rbind(powerdata,powerdata.)
  power_df = inner_join(powerdata,df,by = "index")
}

str(power_df)
ggplot(power_df,aes(x = sample_size, y = std_dev, fill = power)) + geom_tile() + scale_fill_gradient(low="white", high="blue") +
  theme_classic()

ggplot(power_df,aes(x = sample_size, y = n_genotypes, fill = power)) + geom_tile() + scale_fill_gradient(low="white", high="blue") +
  theme_classic()

ggplot(power_df,aes(x = n_genotypes, y = std_dev, fill = power)) + geom_tile() + scale_fill_gradient(low="white", high="blue") +
  theme_classic()
