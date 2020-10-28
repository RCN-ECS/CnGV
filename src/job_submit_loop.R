######################
## Cluster For Loop ##
######################

require(tidyr)
require(ggplot2)
require(dplyr)

# Read in Parameter table
df1 = read.csv("df.csv")

chunk_size <- 1000

df1$chunks <- ggplot2::cut_interval(1:length(unique(df1$row)), length=chunk_size, labels=FALSE)
#chunks1 = data.frame("chunks" = chunks, group = unique(df1$group))

#job_df = full_join(chunks1, df1, by = "group")

for(j in 1:length(unique(df1$chunks))){
    df = df1[df1$chunks == j,]
   
   # Create and submit job for each row
   for(i in 1:nrow(df)){
     id = unique(df$row)[i]
     filename <- id 
     #fileConn<-file(print(paste0(filename,".bash")))
     fileConn<-file(paste0(filename,".bash"))
     
     writeLines(c("#!/bin/bash",
              	##"#SBATCH --reservation=lotterhos",
                "#SBATCH --nodes=1",
                "#SBATCH --tasks-per-node=1",
                paste0("#SBATCH --job-name=",filename,".txt"),
                "#SBATCH --mem=500Mb",
                "#SBATCH --mail-user=m.albecker@northeastern.edu",
                "#SBATCH --mail-type=FAIL",
                ## "#SBATCH --partition=lotterhos",
                "#SBATCH --partition=short",
  	          	"#SBATCH --time=03:00:00",
                ##paste0("#SBATCH --output=",filename,".output"),
                ##paste0("#SBATCH --error=",filename,".error"),
                paste0("Rscript --vanilla Cov_GxE_clusterFun.R ",
                       df$row[i]," ",
                       df$n_pop[i]," ",
                       df$sample_size[i]," ",
                       df$std_dev[i]," ",
                       df$n_env[i]," ",
                       df$delta_env[i]," ",
                       df$delta_gen[i]," ",
                       df$interaction[i]," ",
                       df$errpop[i]," ",
                       df$replicate[i]," ",
                       df$env_scenario[i]," ",
                       df$seed[i])
   ), fileConn)
     system(paste0("sbatch ",filename,".bash"))
     Sys.sleep(5) 
}
Sys.sleep(14400) # 4 hours
}
