## Testing Unbalanced Designs
library(emmeans) 
library(tibble)
library(dplyr)
library(ggplot2)


###### Simulating Unbalanced Designs  #########

emm1 = filter(emms, iteration == 1)
neworder <- c("Pop_balanc","Samp_balanc", "Pop_1Remov" , "Samp_1Remov", "Pop_2Remov" , "Samp_2Remov","Pop_3Remov" , "Samp_3Remov" )
emm2 <- arrange(transform(emm1,
                          categorySplit=factor(categorySplit,levels=neworder)),categorySplit)
emm1$categorySplit = str_split_fixed(emm1$category, "ed", 2)[,1]
emm1$subgroup = str_split_fixed(emm1$category, "ed", 2)[,2]
emm1$subsub = paste(emm1$subgroup, emm1$group, sep = "_")
emm1$subsub_replaced = NULL
emm1 = emm1 %>%
    mutate(subsub_replaced = ifelse(subsub == "_A", "A Group",
                                    ifelse(subsub == "_B", "B Group",
                                           ifelse(subsub == "_A_A", "A Group \n A Removal",
                                                  ifelse(subsub == "_A_B", "A Group \n B Removal",
                                                         ifelse(subsub == "_B_A", "B Group \n A Removal",  "B Group \n B Removal"))))))
emm1$subsub = toupper(emm1$subsub)

emm1$group = NULL
emm1 <- emm1 %>%
    mutate(group = ifelse(gen_factor %in% c("G_1", "G_2", "G_3", "G_4"), "A", "B"))
ggplot(emm2, aes(x = subsub_replaced, y = emmean, group = subsub, fill = categorySplit))+
    xlab("Grouping")+ylab("Estimated Marginal Mean")+
    geom_boxplot()+theme_bw() + facet_wrap(~categorySplit,ncol =2 ) + theme(legend.position = "none")

Delta_gen = seq(from = -1, to = 1, length.out = 50)

results = results. = emm = emm. = data.frame()

source("Cov_GxE_functions.R")

for(i in 1:length(Delta_gen)){
    for(j in 1:50){

#df1 <- df.foundations(1, Delta_gen[i], 8, 2, 1, N_gen[j], 0, 2, 3)
df1 <- df.foundations(1, Delta_gen[i], 8, 2, 1, 8, 0, i, (i+100)) #(delta_env, delta_gen, sample_size, n_env, std_dev, n_pop, interaction, seed1, seed2)

# Create Data
df <- df1[[1]] # Sample 
df.ne <- df1[[3]] # Population

# Balanced Designs
m1 <- mod.GxE(df, is.perm = FALSE) # Sample
cov_matrix <- m1[[1]]
emm_df <- m1[[12]]
emm_df$category <- rep("Samp_balanced", nrow(emm_df))

cov_corrected1 = round(cov.function(cov_matrix, emm_df, is.sample = TRUE),3)

m1.ne <- mod.GxE(df.ne, is.perm = FALSE) # Population
cov_matrix.ne <- m1.ne[[1]]
emm_df.ne <- m1.ne[[12]]
emm_df.ne$category <- rep("Pop_balanced", nrow(emm_df.ne))
cov_corrected1.ne = round(cov.function(cov_matrix.ne, emm_df.ne, is.sample = FALSE),3)

# Unbalanced Designs - Remove 1 Genotype 

G_lev = levels(df$gen_factor)

df2a <- filter(df, gen_factor != c(G_lev[1]))
m2a <- mod.GxE(df2a, is.perm = FALSE) 
cov_matrix2a <- m2a[[1]]
emm_df2a <- m2a[[12]]
emm_df2a$category <- rep("Samp_1Removed_a", nrow(emm_df2a))
cov_corrected2a = round(cov.function(cov_matrix2a, emm_df2a, is.sample = TRUE),3)

df2b <- filter(df, gen_factor != c(G_lev[5]))
m2b <- mod.GxE(df2b, is.perm = FALSE) 
cov_matrix2b <- m2b[[1]]
emm_df2b <- m2b[[12]]
emm_df2b$category <- rep("Samp_1Removed_b", nrow(emm_df2b))
cov_corrected2b = round(cov.function(cov_matrix2b, emm_df2b, is.sample = TRUE),3)

df2a.ne <- filter(df.ne, gen_factor != c(G_lev[1]))
m2a.ne <- mod.GxE(df2a.ne, is.perm = FALSE) 
cov_matrix2a.ne <- m2a.ne[[1]]
emm_df2a.ne <- m2a.ne[[12]]
emm_df2a.ne$category <- rep("Pop_1Removed_a", nrow(emm_df2a.ne))
cov_corrected2a.ne = round(cov.function(cov_matrix2a.ne, emm_df2a.ne, is.sample = FALSE),3)

df2b.ne <- filter(df.ne, gen_factor != c(G_lev[5]))
m2b.ne <- mod.GxE(df2b.ne, is.perm = FALSE) 
cov_matrix2b.ne <- m2b.ne[[1]]
emm_df2b.ne <- m2b.ne[[12]]
emm_df2b.ne$category <- rep("Pop_1Removed_b", nrow(emm_df2b.ne))
cov_corrected2b.ne = round(cov.function(cov_matrix2b.ne, emm_df2b.ne, is.sample = FALSE),3)

# Unbalanced Designs - Remove 2 Genotypes 

'%notin%' <- Negate(`%in%`)
df3a <- filter(df, gen_factor %notin% c(G_lev[1],G_lev[2]))
m3a <- mod.GxE(df3a, is.perm = FALSE) 
cov_matrix3a <- m3a[[1]]
emm_df3a <- m3a[[12]]
emm_df3a$category <- rep("Samp_2Removed_a", nrow(emm_df3a))
cov_corrected3a = round(cov.function(cov_matrix3a, emm_df3a, is.sample = TRUE),3)

df3b <- filter(df, gen_factor %notin% c(G_lev[5],G_lev[6]))
m3b <- mod.GxE(df3b, is.perm = FALSE) 
cov_matrix3b <- m3b[[1]]
emm_df3b <- m3b[[12]]
emm_df3b$category <- rep("Samp_2Removed_b", nrow(emm_df3b))
cov_corrected3b = round(cov.function(cov_matrix3b, emm_df3b, is.sample = TRUE),3)

df3a.ne <- filter(df.ne, gen_factor %notin% c(G_lev[1],G_lev[2]))
m3a.ne <- mod.GxE(df3a.ne, is.perm = FALSE) 
cov_matrix3a.ne <- m3a.ne[[1]]
emm_df3a.ne <- m3a.ne[[12]]
emm_df3a.ne$category <- rep("Pop_2Removed_a", nrow(emm_df3a.ne))
cov_corrected3a.ne = round(cov.function(cov_matrix3a.ne, emm_df3a.ne, is.sample = FALSE),3)

df3b.ne <- filter(df.ne, gen_factor %notin% c(G_lev[5],G_lev[6]))
m3b.ne <- mod.GxE(df3b.ne, is.perm = FALSE) 
cov_matrix3b.ne <- m3b.ne[[1]]
emm_df3b.ne <- m3b.ne[[12]]
emm_df3b.ne$category <- rep("Pop_2Removed_b", nrow(emm_df3b.ne))
cov_corrected3b.ne = round(cov.function(cov_matrix3b.ne, emm_df3b.ne, is.sample = FALSE),3)

# Unbalanced Designs - Remove 3 Genotypes 

df4a <- filter(df, gen_factor %notin% c(G_lev[1],G_lev[2],G_lev[3]))
m4a <- mod.GxE(df4a, is.perm = FALSE) 
cov_matrix4a <- m4a[[1]]
emm_df4a <- m4a[[12]]
emm_df4a$category <- rep("Samp_3Removed_a", nrow(emm_df4a))
cov_corrected4a = round(cov.function(cov_matrix4a, emm_df4a, is.sample = TRUE),3)

df4b <- filter(df, gen_factor %notin% c(G_lev[5],G_lev[6],G_lev[7]))
m4b <- mod.GxE(df4b, is.perm = FALSE) 
cov_matrix4b <- m4b[[1]]
emm_df4b <- m4b[[12]]
emm_df4b$category <- rep("Samp_3Removed_b", nrow(emm_df4b))
cov_corrected4b = round(cov.function(cov_matrix4b, emm_df4b, is.sample = TRUE),3)

df4a.ne <- filter(df.ne, gen_factor %notin% c(G_lev[1],G_lev[2],G_lev[3]))
m4a.ne <- mod.GxE(df4a.ne, is.perm = FALSE) 
cov_matrix4a.ne <- m4a.ne[[1]]
emm_df4a.ne <- m4a.ne[[12]]
emm_df4a.ne$category <- rep("Pop_3Removed_a", nrow(emm_df4a.ne))
cov_corrected4a.ne = round(cov.function(cov_matrix4a.ne, emm_df4a.ne, is.sample = FALSE),3)

df4b.ne <- filter(df.ne, gen_factor %notin% c(G_lev[5],G_lev[6],G_lev[7]))
m4b.ne <- mod.GxE(df4b.ne, is.perm = FALSE) 
cov_matrix4b.ne <- m4b.ne[[1]]
emm_df4b.ne <- m4b.ne[[12]]
emm_df4b.ne$category <- rep("Pop_3Removed_b", nrow(emm_df4b.ne))
cov_corrected4b.ne = round(cov.function(cov_matrix4b.ne, emm_df4b.ne, is.sample = FALSE),3)

emms1 = rbind( emm_df, emm_df.ne,emm_df2a,emm_df2b, emm_df2a.ne,emm_df2b.ne,
                   emm_df3a, emm_df3b, emm_df3a.ne, emm_df3b.ne,
                  emm_df4a, emm_df4b, emm_df4a.ne,emm_df4b.ne)
emms1$iteration = rep(j, nrow(emms1))
emms1$delta_gen = rep(i, nrow(emms1))

emm. = rbind(emms1, emm.)

results1 = data.frame("iteration" =  j,
                      "Delta_gen" = i, 
                    "CovSample_Balanced" = cov_corrected1,
                    "CovPop_Balanced" = cov_corrected1.ne,
                    "Cov_Samp_Unbalanced-1Genotype.a" = cov_corrected2a,
                    "Cov_Samp_Unbalanced-1Genotype.b" = cov_corrected2b,
                    "Cov_Pop_Unbalanced-1Genotype.a" = cov_corrected2a.ne,
                    "Cov_Pop_Unbalanced-1Genotype.b" = cov_corrected2b.ne,
                    "Cov_Samp_Unbalanced-2Genotype.a" = cov_corrected3a,
                    "Cov_Samp_Unbalanced-2Genotype.b" = cov_corrected3b,
                    "Cov_Pop_Unbalanced-2Genotype.a" = cov_corrected3a.ne,
                    "Cov_Pop_Unbalanced-2Genotype.b" = cov_corrected3b.ne,
                    "Cov_Samp_Unbalanced-3Genotype.a" = cov_corrected4a,
                    "Cov_Samp_Unbalanced-3Genotype.b" = cov_corrected4b,
                    "Cov_Pop_Unbalanced-3Genotype.a" = cov_corrected4a.ne,
                    "Cov_Pop_Unbalanced-3Genotype.b" = cov_corrected4b.ne)
results. = rbind(results1, results.)
}
    results = rbind(results., results)
    emm = rbind(emm., emm)
}

write.csv(results, "results_unbalanced.csv")
write.csv(emm, "emms_unbalanced.csv")


results = read.csv("~/Desktop/Unbalanced_results.csv")
emms = read.csv("~/Desktop/emms_unbalanced.csv")

###### Effect of Unbalanced Designs on Estimated Marginal Means #########

emm1 = filter(emms, iteration == 1)
neworder <- c("Pop_balanc","Samp_balanc", "Pop_1Remov" , "Samp_1Remov", "Pop_2Remov" , "Samp_2Remov","Pop_3Remov" , "Samp_3Remov" )
emm2 <- arrange(transform(emm1,
                          categorySplit=factor(categorySplit,levels=neworder)),categorySplit)
emm1$categorySplit = str_split_fixed(emm1$category, "ed", 2)[,1]
emm1$subgroup = str_split_fixed(emm1$category, "ed", 2)[,2]
emm1$subsub = paste(emm1$subgroup, emm1$group, sep = "_")
emm1$subsub_replaced = NULL
emm1 = emm1 %>%
    mutate(subsub_replaced = ifelse(subsub == "_A", "A Group",
                                    ifelse(subsub == "_B", "B Group",
                                           ifelse(subsub == "_A_A", "A Group \n A Removal",
                                                  ifelse(subsub == "_A_B", "A Group \n B Removal",
                                                         ifelse(subsub == "_B_A", "B Group \n A Removal",  "B Group \n B Removal"))))))
emm1$subsub = toupper(emm1$subsub)

emm1$group = NULL
emm1 <- emm1 %>%
    mutate(group = ifelse(gen_factor %in% c("G_1", "G_2", "G_3", "G_4"), "A", "B"))
ggplot(emm2, aes(x = subsub_replaced, y = emmean, group = subsub, fill = categorySplit))+
    xlab("Grouping")+ylab("Estimated Marginal Mean")+
    geom_boxplot()+theme_bw() + facet_wrap(~categorySplit,ncol =2 ) + theme(legend.position = "none")

###### Effect of Unbalanced Designs on CovGE ###########
annotation <- data.frame(
    x = c(0.35, 0.35, 0.35),
    y = c(-.3,-.4,-.5),
    label = c("Balanced Design", "A-side removal","B-side removal")
)
(unbalancedPop1 = ggplot(results, aes(x = CovPop_Balanced)) + 
        geom_point(aes(y = Cov_Pop_Unbalanced.1Genotype.a), colour = "blue")+
       # geom_text(aes(x = 0.2, y= -0.2), colour = "blue", label="A-side removal")+
        geom_point(aes(y = Cov_Pop_Unbalanced.1Genotype.b, group = Delta_gen), colour = "deepskyblue")+
       # geom_text(aes(x = 0.2, y= -0.3), colour = "deepskyblue", label="B-side removal")+
        geom_point(aes(y = CovPop_Balanced),colour = "black") +
        geom_text(data=annotation, aes( x=x, y=y, label=label), 
                  color=c("black", "blue", "deepskyblue"), 
                  size=5)  +    
        ggtitle("Population Effects: 1 Genotype removal")+
        ylab("CovGE") + xlab("CovGE") + theme_bw())

(unbalancedPop2 = ggplot(results, aes(x = CovPop_Balanced)) + 
        geom_point(aes(y = Cov_Pop_Unbalanced.2Genotype.a), colour = "blue")+
        # geom_text(aes(x = 0.2, y= -0.2), colour = "blue", label="A-side removal")+
        geom_point(aes(y = Cov_Pop_Unbalanced.2Genotype.b, group = Delta_gen), colour = "deepskyblue")+
        # geom_text(aes(x = 0.2, y= -0.3), colour = "deepskyblue", label="B-side removal")+
        geom_point(aes(y = CovPop_Balanced),colour = "black") +
        geom_text(data=annotation, aes( x=x, y=y, label=label), 
                  color=c("black", "blue", "deepskyblue"), 
                  size=5)  +    
        ggtitle("Population Effects: 2 Genotype removal")+
        ylab("CovGE") + xlab("CovGE") + theme_bw())

(unbalancedPop3 = ggplot(results, aes(x = CovPop_Balanced)) + 
        geom_point(aes(y = Cov_Pop_Unbalanced.3Genotype.a), colour = "blue")+
        # geom_text(aes(x = 0.2, y= -0.2), colour = "blue", label="A-side removal")+
        geom_point(aes(y = Cov_Pop_Unbalanced.3Genotype.b, group = Delta_gen), colour = "deepskyblue")+
        # geom_text(aes(x = 0.2, y= -0.3), colour = "deepskyblue", label="B-side removal")+
        geom_point(aes(y = CovPop_Balanced),colour = "black") +
        geom_text(data=annotation, aes( x=x, y=y, label=label), 
                  color=c("black", "blue", "deepskyblue"), 
                  size=5)  +    
        ggtitle("Population Effects: 3 Genotype removal")+
        ylab("CovGE") + xlab("CovGE") + theme_bw())

grid.arrange(unbalancedPop1, unbalancedPop2, unbalancedPop3, ncol = 3)

# Pop vs. Sample
(balanced = ggplot(results, aes(x = CovPop_Balanced)) + 
         geom_point(aes(y = CovSample_Balanced, group = Delta_gen), colour = "grey")+
         geom_point(aes(y = CovPop_Balanced),colour = "black") +
         ggtitle("Balanced: Population vs. Sample")+
         ylab("CovGE") + xlab("CovGE") + theme_bw())

(PS_unbalanced1 = ggplot(results, aes(x = CovPop_Balanced)) + 
        geom_point(aes(y = Cov_Samp_Unbalanced.1Genotype.a, group = Delta_gen), alpha = 0.4, colour = "blue")+
        geom_point(aes(y = Cov_Samp_Unbalanced.1Genotype.b, group = Delta_gen), alpha = 0.4, colour ="deepskyblue")+
        geom_point(aes(y = CovPop_Balanced),colour = "black") +
        geom_text(data=annotation, aes( x=x, y=y, label=label), 
                  color=c("black", "blue", "deepskyblue"), 
                  size=5)  +  
        ggtitle("Population vs. Sample: 1 Genotype Removal")+
        ylab("CovGE") + xlab("CovGE") + theme_bw())

(PS_unbalanced2 = ggplot(results, aes(x = CovPop_Balanced)) + 
        geom_point(aes(y = Cov_Samp_Unbalanced.2Genotype.a, group = Delta_gen), alpha = 0.4, colour = "blue")+
        geom_point(aes(y = Cov_Samp_Unbalanced.2Genotype.b, group = Delta_gen), alpha = 0.4, colour ="deepskyblue")+
        geom_point(aes(y = CovPop_Balanced),colour = "black") +
        geom_text(data=annotation, aes( x=x, y=y, label=label), 
                  color=c("black", "blue", "deepskyblue"), 
                  size=5)  +  
        ggtitle("Population vs. Sample: 2 Genotype Removal")+
        ylab("CovGE") + xlab("CovGE") + theme_bw())

(PS_unbalanced3 = ggplot(results, aes(x = CovPop_Balanced)) + 
        geom_point(aes(y = Cov_Samp_Unbalanced.3Genotype.a, group = Delta_gen), alpha = 0.4, colour = "blue")+
        geom_point(aes(y = Cov_Samp_Unbalanced.3Genotype.b, group = Delta_gen), alpha = 0.4, colour ="deepskyblue")+
        geom_point(aes(y = CovPop_Balanced),colour = "black") +
        geom_text(data=annotation, aes( x=x, y=y, label=label), 
                  color=c("black", "blue", "deepskyblue"), 
                  size=5)  +  
        ggtitle("Population vs. Sample: 3 Genotype Removal")+
        ylab("CovGE") + xlab("CovGE") + theme_bw())
 grid.arrange(balanced,PS_unbalanced1,PS_unbalanced2,PS_unbalanced3, ncol = 4)         
