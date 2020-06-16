
#choose this for a 2x2 case with GxE
n=2
phen <- c(-1,1,1,-1) 

#choose this for a 2x2 case without GxE
n=2
phen <- c(-1,1,-1,1)+rnorm(n*n,0,0.05)

#choose this for an "n" case with GxE
n=3
phen <- rnorm(n*n, 0, 1) # GxE

#choose this for an "n" case with cov(G,E) but no GxE
# this will illustrate how the permutation creates GxE
phen <- c(1:n+rnorm(n*n,0,0.1)) #Cov GE


########## RUN HERE TO CREATE DATA ###############
input_df <- data.frame(avg_phen = phen,
  avg_phen_corrected=((phen-mean(phen))/sd(phen)),#rnorm(n*n),
                       se <- 0.1,
                       exp_env_factor=rep(paste0("E",1:n),n),
                       gen_factor=rep(paste0("G",1:n), each=n))
boxplot(avg_phen_corrected~exp_env_factor*gen_factor, data=input_df)
input_df


######### FUNCTIONS #####################

##########################################
#### Molly's function for mean.GxE ####
##########################################
mean.GxE <- function(input_df){ # input is mean_df
  
  # Clear outputs
  allGEmeans <- c()
  GxE_mean.temp <- c()
  
  # Means of Means
  E_means <- tapply(input_df$avg_phen_corrected, input_df$exp_env_factor, mean)
  G_means <- tapply(input_df$avg_phen_corrected, input_df$gen_factor, mean)
  Gmean_mat <- data.frame("G_means" = G_means, "gen_factor" = unique(input_df$gen_factor))
  Emean_mat <- data.frame("E_means" = E_means, "exp_env_factor" = unique(input_df$exp_env_factor))
  
  # Match means to native
  #Cov_mean_matrix = Gmean_mat
  #Cov_mean_matrix$exp_env_factor <- input_df$nat_env_factor[match(Cov_mean_matrix$gen_factor,input_df$gen_factor)]
  #Cov_mean_matrix$E_means <- Emean_mat$E_means[match(Cov_mean_matrix$exp_env_factor,Emean_mat$exp_env_factor)]
  
  # Magnitude of GxE -- Loop -- Means
  for (i in 1:nlevels(input_df$gen_factor)){
    for (j in 1:nlevels(input_df$exp_env_factor)){
      G_levels <- levels(input_df$gen_factor)
      E_levels <- levels(input_df$exp_env_factor)
      GxE_mean.temp <- abs(input_df$avg_phen_corrected[input_df$gen_factor == G_levels[i] & input_df$exp_env_factor == E_levels[j]] - # GxE (Phenotype of ith genotype in jth environment)
                             mean(input_df$avg_phen_corrected[input_df$gen_factor == G_levels[i]])- # mean phenotype of ith Genotype
                             mean(input_df$avg_phen_corrected[input_df$exp_env_factor == E_levels[j]])+ # mean phenotype of jth Environment
                             mean(input_df$avg_phen_corrected)) # Overall mean
      allGEmeans <- c(allGEmeans, GxE_mean.temp)
    }
  }
  #hist(allGEmeans)
  GxE_means = mean(allGEmeans)
  
  #return(list(Cov_mean_matrix, GxE_means, allGEmeans))
  return(GxE_means) #KEL edit
}
mean.GxE(input_df)

##########################################
############ Our current function for permutation #########
##########################################
# For GxE, this type of permutation CREATES datasets with GxE
# when you think about it, GxE is random deviations from that 
# expected based on the G_means and E_means. When we shuffle
# the identities of G_means and E_means, we create random
# deviations from that expected based on G_means and E_means.
# In a way, this type of permutation tests creates a null distrubtion
# for the hypothesis that there IS GxE in the data
permutation_means <- function(input_df){ # means dataframe (mean_df)
  # Clear outputs
  perm_means <- data.frame()
  null_means. = null_means = NULL
  
  # Shuffle means data
  # KEL EDITS
  null_means. <- rnorm(nrow(input_df), mean = input_df$avg_phen, sd = input_df$se) # create replicate mean
  null_means <- sample(null_means., size=length(null_means.), replace = FALSE) # shuffle means without replacement
  
  # KEL EDITS
  perm_means <- data.frame("gen_factor" = input_df$gen_factor,
                           "exp_env_factor" = input_df$exp_env_factor,
                           #"nat_env_factor" = input_df$nat_env_factor,
                           "avg_phen" = null_means) #KEL EDIT
  # You were already sampling from the 
  
  # Restandardize
  perm_means$avg_phen_corrected = (perm_means$avg_phen - mean(perm_means$avg_phen))/sd(perm_means$avg_phen)
  
  return(perm_means)
}

null <- replicate(1000, mean.GxE(permutation_means(input_df)))
hist(null, breaks=seq(0,1,0.01))
abline(v=mean.GxE(input_df), col="blue")
rank(c(mean.GxE(input_df), null))[1]/length(null+1)
# ASK MOLLY HOW RANKS ARE HANDELED IN P-VALUES


##########################################
############ A solution? #########
##########################################
##################################
##  the null hypothesis is that the GiEj-Gimean-Ejmean=0
# so the question is how big would we expect an interaction to be 
# by random chance given the variation we see around GiEj?
# the answer is rnorm(GiEj, sd=se)-G1mean-E1mean
# where "se" is the SE observed around the mean of the GiEj level,
# and although not shown, the G1means and E1means also take into
# account the variation around their means in the calc
##################################
permutation_means_GxE <- function(input_df){
  mean(input_df$avg_phen_corrected)
  
  E_means <- tapply(input_df$avg_phen_corrected, input_df$exp_env_factor, mean)
  G_means <- tapply(input_df$avg_phen_corrected, input_df$gen_factor, mean)
  
  allGEmeans <- c()
  
  for (i in 1:nlevels(input_df$gen_factor)){
    for (j in 1:nlevels(input_df$exp_env_factor)){
      G_levels <- levels(input_df$gen_factor)
      E_levels <- levels(input_df$exp_env_factor)
    
      GiEj_mean <- input_df$avg_phen_corrected[input_df$gen_factor == G_levels[i] & input_df$exp_env_factor == E_levels[j]]
      Gi_mean <- mean(input_df$avg_phen_corrected[input_df$gen_factor == G_levels[i]])
      Ej_mean <- mean(input_df$avg_phen_corrected[input_df$exp_env_factor == E_levels[j]])
      
      # Create a sample of the null expectation for the GiEj
      GiEj_null_samp <- rnorm(1, mean=(Gi_mean+Ei_mean),
                         sd=mean(input_df$se[input_df$gen_factor == G_levels[i] & input_df$exp_env_factor == E_levels[j]]))
      
      # Create a sample of mean Gi based on the average se across environmments
      # A rough approximation to what we should do, which is sample each G and then calculate the mean
      #Gi_mean_samp <- rnorm(1, mean=Gi_mean,
      #                   sd= mean(input_df$se[input_df$gen_factor == G_levels[i]]))
      
      # Create a sample of mean Ej based on the average se across genotypes
      # A rough approximation to what we should do, which is sample each G and then calculate the mean
      #Ej_mean_samp <- rnorm(1, mean=Ei_mean,
      #                 sd=mean(input_df$se[input_df$exp_env_factor == E_levels[j]]))
      
      GxE_mean.temp <- abs(GiEj_null_samp-
        # GxE (Phenotype of ith genotype in jth environment)
          Gi_mean- # mean phenotype of ith Genotype
          Ej_mean                   + # mean phenotype of jth Environment
                             mean(input_df$avg_phen_corrected)# Overall mean
        ) 
      allGEmeans <- c(allGEmeans, GxE_mean.temp)
    }
  }
  GxE_means=mean(allGEmeans)
  return(GxE_means)
}

null2 <- replicate(1000, permutation_means_GxE(input_df))
hist(null2, breaks=seq(0,1,0.01))
abline(v=mean.GxE(input_df), col="blue")
(P<-1-(rank(c(mean.GxE(input_df), null2))[1]/(length(null2)+1)))

