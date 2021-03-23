## Testing Unbalanced Designs
library(emmeans) 
library(tibble)
library(dplyr)
library(ggplot2)

Delta_gen = seq(from = -1, to = 1, length.out = 100)


results = results. = data.frame()
source("~/Documents/GitHub/CnGV/CnGV/src/Cov_GxE_functions.R")

for(j in 1:100){
for(i in 1:length(Delta_gen)){

#df1 <- df.foundations(1, Delta_gen[i], 8, 2, 1, N_gen[j], 0, 2, 3)
df1 <- df.foundations(1, Delta_gen[i], 8, 2, 1, 8, 0, i, (i+100)) #(delta_env, delta_gen, sample_size, n_env, std_dev, n_pop, interaction, seed1, seed2)

# Create Data
df <- df1[[1]] # Sample 
df.ne <- df1[[3]] # Population

# Balanced Designs
m1 <- mod.GxE(df, is.perm = FALSE) # Sample
cov_matrix <- m1[[1]]
emm_df <- m1[[12]]
cov_corrected1 = round(cov.function(cov_matrix, emm_df, is.sample = TRUE),3)

m1.ne <- mod.GxE(df.ne, is.perm = FALSE) # Population
cov_matrix.ne <- m1.ne[[1]]
emm_df.ne <- m1.ne[[12]]
cov_corrected1.ne = round(cov.function(cov_matrix.ne, emm_df.ne, is.sample = FALSE),3)

# Unbalanced Designs - Remove 1 Genotype 

G_lev = levels(df$gen_factor)

df2a <- filter(df, gen_factor != c(G_lev[1]))
m2a <- mod.GxE(df2a, is.perm = FALSE) 
cov_matrix2a <- m2a[[1]]
emm_df2a <- m2a[[12]]
cov_corrected2a = round(cov.function(cov_matrix2a, emm_df2a, is.sample = TRUE),3)

df2b <- filter(df, gen_factor != c(G_lev[5]))
m2b <- mod.GxE(df2b, is.perm = FALSE) 
cov_matrix2b <- m2b[[1]]
emm_df2b <- m2b[[12]]
cov_corrected2b = round(cov.function(cov_matrix2b, emm_df2b, is.sample = TRUE),3)

df2a.ne <- filter(df.ne, gen_factor != c(G_lev[1]))
m2a.ne <- mod.GxE(df2a.ne, is.perm = FALSE) 
cov_matrix2a.ne <- m2a.ne[[1]]
emm_df2a.ne <- m2a.ne[[12]]
cov_corrected2a.ne = round(cov.function(cov_matrix2a.ne, emm_df2a.ne, is.sample = FALSE),3)

df2b.ne <- filter(df.ne, gen_factor != c(G_lev[5]))
m2b.ne <- mod.GxE(df2b.ne, is.perm = FALSE) 
cov_matrix2b.ne <- m2b.ne[[1]]
emm_df2b.ne <- m2b.ne[[12]]
cov_corrected2b.ne = round(cov.function(cov_matrix2b.ne, emm_df2b.ne, is.sample = FALSE),3)

# Unbalanced Designs - Remove 2 Genotypes 

'%notin%' <- Negate(`%in%`)
df3a <- filter(df, gen_factor %notin% c(G_lev[1],G_lev[2]))
m3a <- mod.GxE(df3a, is.perm = FALSE) 
cov_matrix3a <- m3a[[1]]
emm_df3a <- m3a[[12]]
cov_corrected3a = round(cov.function(cov_matrix3a, emm_df3a, is.sample = TRUE),3)

df3b <- filter(df, gen_factor %notin% c(G_lev[5],G_lev[6]))
m3b <- mod.GxE(df3b, is.perm = FALSE) 
cov_matrix3b <- m3b[[1]]
emm_df3b <- m3b[[12]]
cov_corrected3b = round(cov.function(cov_matrix3b, emm_df3b, is.sample = TRUE),3)

df3a.ne <- filter(df.ne, gen_factor %notin% c(G_lev[1],G_lev[2]))
m3a.ne <- mod.GxE(df3a.ne, is.perm = FALSE) 
cov_matrix3a.ne <- m3a.ne[[1]]
emm_df3a.ne <- m3a.ne[[12]]
cov_corrected3a.ne = round(cov.function(cov_matrix3a.ne, emm_df3a.ne, is.sample = FALSE),3)

df3b.ne <- filter(df.ne, gen_factor %notin% c(G_lev[5],G_lev[6]))
m3b.ne <- mod.GxE(df3b.ne, is.perm = FALSE) 
cov_matrix3b.ne <- m3b.ne[[1]]
emm_df3b.ne <- m3b.ne[[12]]
cov_corrected3b.ne = round(cov.function(cov_matrix3b.ne, emm_df3b.ne, is.sample = FALSE),3)

# Unbalanced Designs - Remove 3 Genotypes 

df4a <- filter(df, gen_factor %notin% c(G_lev[1],G_lev[2],G_lev[3]))
m4a <- mod.GxE(df4a, is.perm = FALSE) 
cov_matrix4a <- m4a[[1]]
emm_df4a <- m4a[[12]]
cov_corrected4a = round(cov.function(cov_matrix4a, emm_df4a, is.sample = TRUE),3)

df4b <- filter(df, gen_factor %notin% c(G_lev[5],G_lev[6],G_lev[7]))
m4b <- mod.GxE(df4b, is.perm = FALSE) 
cov_matrix4b <- m4b[[1]]
emm_df4b <- m4b[[12]]
cov_corrected4b = round(cov.function(cov_matrix4b, emm_df4b, is.sample = TRUE),3)

df4a.ne <- filter(df.ne, gen_factor %notin% c(G_lev[1],G_lev[2],G_lev[3]))
m4a.ne <- mod.GxE(df4a.ne, is.perm = FALSE) 
cov_matrix4a.ne <- m4a.ne[[1]]
emm_df4a.ne <- m4a.ne[[12]]
cov_corrected4a.ne = round(cov.function(cov_matrix4a.ne, emm_df4a.ne, is.sample = FALSE),3)

df4b.ne <- filter(df.ne, gen_factor %notin% c(G_lev[5],G_lev[6],G_lev[7]))
m4b.ne <- mod.GxE(df4b.ne, is.perm = FALSE) 
cov_matrix4b.ne <- m4b.ne[[1]]
emm_df4b.ne <- m4b.ne[[12]]
cov_corrected4b.ne = round(cov.function(cov_matrix4b.ne, emm_df4b.ne, is.sample = FALSE),3)

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
results. = rbind(results1, results)
}
    results = rbind(results., results)
}

results = read.csv("~/Desktop/Unbalanced_results.csv")

# Balanced vs. 1 Removal
(minus1 = ggplot(results, aes(x = CovPop_Balanced)) + 
         geom_point(aes(y = CovPop_Balanced),colour = "black")+
         ggtitle("Population: Balanced vs. 1 Removal from each genotype group ")+
         geom_point(aes(y = Cov_Pop_Unbalanced.1Genotype.a),colour = "red")+
         geom_point(aes(y = Cov_Pop_Unbalanced.1Genotype.b),colour = "blue")+ 
         geom_text(aes(x=CovPop_Balanced, y=Cov_Pop_Unbalanced.1Genotype.b+0.001), family = "Times", colour = "blue", label="B-side removal")+
         geom_text(aes(x=CovPop_Balanced, y=CovPop_Balanced+0.001), colour = "black", label="Balanced - No Removal")+
         geom_text(aes(x=CovPop_Balanced, y=Cov_Pop_Unbalanced.1Genotype.a+0.001), colour = "red", label="A-side removal")+
         ylab("CovGE of Population") + xlab("CovGE of Population") + theme_bw())

(minus2 = ggplot(results, aes(x = CovPop_Balanced)) + 
    geom_point(aes(y = CovPop_Balanced),colour = "black")+
    ggtitle("Population: Balanced vs. 2 Removals from each genotype group ")+
    geom_point(aes(y = Cov_Pop_Unbalanced.2Genotype.a),colour = "red")+
    geom_point(aes(y = Cov_Pop_Unbalanced.2Genotype.b),colour = "blue")+ 
    geom_text(aes(x=CovPop_Balanced, y=Cov_Pop_Unbalanced.2Genotype.b+0.002), family = "Times", colour = "blue", label="B-side removal")+
    geom_text(aes(x=CovPop_Balanced, y=CovPop_Balanced+0.002), colour = "black", label="'+' is Balanced - No Removal")+
    geom_text(aes(x=CovPop_Balanced, y=Cov_Pop_Unbalanced.2Genotype.a+0.002), family = "Times", colour = "red", label="A-side removal")+
    ylab("CovGE of Population") + xlab("CovGE of Population") + theme_bw())

(minus3 = ggplot(results, aes(x = CovPop_Balanced)) + 
    geom_point(aes(y = CovPop_Balanced),colour = "black")+
    ggtitle("Population: Balanced vs. 3 Removals from each genotype group ")+
    geom_point(aes(y = Cov_Pop_Unbalanced.3Genotype.a),colour = "red")+
    geom_point(aes(y = Cov_Pop_Unbalanced.3Genotype.b),colour = "blue")+ 
    geom_text(aes(x=CovPop_Balanced, y=Cov_Pop_Unbalanced.3Genotype.b+0.002), family = "Times", colour = "blue", label="B-side removal")+
    geom_text(aes(x=CovPop_Balanced, y=CovPop_Balanced+0.002), colour = "black", label="Balanced - No Removal")+
    geom_text(aes(x=CovPop_Balanced, y=Cov_Pop_Unbalanced.3Genotype.a+0.002), family = "Times", colour = "red", label="A-side removal")+
    ylab("CovGE of Population") + xlab("CovGE of Population") + theme_bw())

grid.arrange(minus1, minus2, minus3,ncol = 3)

(samppopminus1 = ggplot(results) + 
    ggtitle("Sample vs.Population: 1 Genotype Removal")+
    geom_boxplot(aes(x = -1, y = Cov_Samp_Unbalanced.1Genotype.a), fill = "red", alpha = 0.5,  colour = "red")+
    geom_boxplot(aes(x = 1, y = Cov_Samp_Unbalanced.1Genotype.b), fill = "blue", alpha = 0.5, colour = "blue")+ 
    geom_point(aes(x = 0, y = CovPop_Balanced),shape = 3, size = 4, colour = "black")+
    geom_text(aes(x = 1, y= -0.45), family = "Times", colour = "blue", label="B-side removal")+
    geom_text(aes(x = -0.2, y = -.2), colour = "black", label="Balanced - No Removal")+
    geom_text(aes(x = -1, y = -0.55),colour = "red", label="A-side removal")+
    xlab(" ")+ylab("CovGE Estimate") + theme_bw()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()))

(samppopminus2 = ggplot(results) + 
    ggtitle("Sample vs.Population: 2 Genotype Removal")+
    geom_boxplot(aes(x = -1, y = Cov_Samp_Unbalanced.2Genotype.a), fill = "red", alpha = 0.5,  colour = "red")+
    geom_boxplot(aes(x = 1, y = Cov_Samp_Unbalanced.2Genotype.b), fill = "blue", alpha = 0.5, colour = "blue")+ 
    geom_point(aes(x = 0, y = CovPop_Balanced),shape = 3, size = 4, colour = "black")+
    geom_text(aes(x = 1, y= -0.45), family = "Times", colour = "blue", label="B-side removal")+
    geom_text(aes(x = -0.2, y = -.2), colour = "black", label="Balanced - No Removal")+
    geom_text(aes(x = -1, y = -0.55),colour = "red", label="A-side removal")+
    xlab(" ")+ylab("CovGE Estimate") + theme_bw()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()))

(samppopminus3 = ggplot(results) + 
    ggtitle("Sample vs.Population: 3 Genotype Removals")+
    geom_boxplot(aes(x = -1, y = Cov_Samp_Unbalanced.3Genotype.a), fill = "red", alpha = 0.5,  colour = "red")+
    geom_boxplot(aes(x = 1, y = Cov_Samp_Unbalanced.3Genotype.b), fill = "blue", alpha = 0.5, colour = "blue")+ 
    geom_point(aes(x = 0, y = CovPop_Balanced),shape = 3, size = 4, colour = "black")+
    geom_text(aes(x = 1, y= -0.45), family = "Times", colour = "blue", label="B-side removal")+
    geom_text(aes(x = -0.2, y = -.2), colour = "black", label="Balanced - No Removal")+
    geom_text(aes(x = -1, y = -0.55),colour = "red", label="A-side removal")+
    xlab(" ")+ylab("CovGE Estimate") + theme_bw()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()))
library(gridExtra)
grid.arrange(samppopminus1, samppopminus2,samppopminus3, ncol = 3)
