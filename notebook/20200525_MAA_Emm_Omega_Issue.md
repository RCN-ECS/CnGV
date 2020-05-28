#### Mo' Problemz. 

When I plot the magnitude of GxE using emms vs omega^2, here is the result. Purple = Omega^2 is significant, grey = not signficant. 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.25.OmegaVsEmm.png)

If I filter this data to just Omega^2 estimates above 0.9 and EMM estimates below 0.1, I get 26 hits. The parameters that led to that outcome are: 

delta_env |delta_gen| sample_size |n_pop |std_dev |interaction
|---|---|---|---|---|---|
   |    0.01   |   0.01   |       10  |  15   |  1.0  |      7.50
   |    0.01  |    0.01    |      20   | 15  |   0.5  |     11.25
   |    0.01   |   0.01      |     5 |   15  |   0.5    |   11.25
   |    0.50   |   0.01     |     20  |  15   |  0.5   |    15.00
   |    0.01    |  0.01     |      5  |  15   |  1.0   |    15.00
   |    0.01     | 0.01     |      5  |  15   |  0.5   |    11.25
   |    0.50     | 0.01     |     10  |  15   |  1.0   |    15.00
   |    0.01     | 0.01     |     10  |  15   |  0.5   |    11.25
   |    0.01     | 0.01     |     10  |  15   |  0.5   |    11.25
   |  0.01  |    0.01       |   20   | 15   |  0.5     |  11.25
   |  0.01   |   0.01       |   20   | 15   |  0.5     |   7.50
   |  0.01    |  0.01       |   20   | 15   |  0.5     |  15.00
   |  0.01  |    0.01       |    5   | 15   |  0.5     |  15.00
   |  0.01   |   0.01       |    5   | 15   |  1.0     |  15.00
   |  0.01   |   0.01       |   10   | 15   |  0.5     |  15.00
   |  0.50   |   0.01       |   20   | 15   |  0.5     |  15.00
   |   0.01  |    0.01      |     5  |  15  |   1.0    |   15.00
   |  0.01   |   0.01       |   10   | 10   |  0.5     |   7.50
   |   1.00  |    0.01      |     5  |  15  |   1.0    |   15.00
   |  0.01    |  0.01       |   20   | 15   |  1.0     |  11.25
   |   0.50   |   0.01      |    20  |  15  |   1.0    |   11.25
   |   0.01   |   0.01      |    10  |  15  |   0.5    |    7.50
   |   0.01   |   0.01      |     5  |  15  |   0.5    |   11.25
   |   0.01   |   0.01      |     5  |  15  |   1.0    |   11.25
   |  0.01     | 0.01       |    5   | 15   |  1.0     |  15.00
   |  0.01     | 0.01       |   20   | 15   |  1.0     |  15.00

Interesting that they are all the same delta_gen... Unsure how that would affect things though. 

Here is the code used to estimate both - based on model: 

```{m1}
    test_temp <- lm(phen_corrected ~ exp_env_factor * gen_factor, data = model_df)

```

Estimated marginal means: 

Formula: GxE_emm = |(Pij - Gi- Ej+ P)| 
```{emm}
# Magnitude of GxE -- EMMs
    GxE_emm <- abs(mean(model_df$phen_corrected) - # Overall mean
                  (emm_G$emmean[emm_G$gen_factor == "G_1"])- # G
                  (emm_E$emmean[emm_E$exp_env_factor == "E_1"])+ # E
                  (emm_GxE[1,3])) # GxE
```
Omega^2: 

Formula: ω2 = (SSeffect – (dfeffect)(MSerror)) / MSerror + SStotal


```{omega}
# Magnitude of GxE -- Omega^2
    w2_GxE = (summary(aov(test_temp))[[1]][3,2] - #(SS_effect -
             (summary(aov(test_temp))[[1]][3,1]*summary(aov(test_temp))[[1]][4,3])) / #(Df_effect * MS_error))/
             (sum(summary(aov(test_temp))[[1]][,2]) + # (SS_total+
             (summary(aov(test_temp))[[1]][4,3])) # MS_error)
```
Here is an example output from an aov for double checking: 

```{output}
summary(aov(test_temp))
                          Df Sum Sq Mean Sq F value  Pr(>F)    
exp_env_factor             1   5.14   5.143   56.18 2.5e-10 ***
gen_factor                 7  66.75   9.536  104.16 < 2e-16 ***
exp_env_factor:gen_factor  7   1.24   0.178    1.94  0.0776 .  
Residuals                 64   5.86   0.092                    
```
### 5/26
As I state above, I narrowed my sims to just those outputs that gave me an Omega^2 >.90 and GxE_Emm <0.1. This yielded 26 results. I picked one to investigate further: 
Parameter set: 
delta_env| delta_gen|sample_size|n_pop|std_dev|interaction
|---|---|---|---|---|---|
| 0.01  | 0.01 |  10 | 15  | 1 | 7.5|

Which resulted in this pattern: 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.26.GxEdiscrepancy.png)
 
 Visually, it LOOKS like there is a GxE. The omega^2 picks this up, but EMM does not. 
 
 I ran the code for this particular set of parameters, but I cannot find any particular issue/error that would drive this discrepancy. I am really at a loss for why there would be such a low EMM GxE. 
 
 My one comforting thought is that out of ~10k sims, there are only 26 that have such an odd discrepancy. Perhaps just a fluke? 
 
 Really uncertain about what to think about this.

## 5.27

Katie recommended to expand to Omega^2 >0.5 and GxE_emm < 0.1. When I do that, 3376 cases are returned. Katie predicted that most would be 10-15 pops. This is partially true: 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.27.numhistogram.png)

Katie also suggested to replot the above case to boxplots and correctly order the X-axis. This yields: 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.27.BoxplotProbPlot.png)

One with n_pop = 5;
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.27.ProbPlot2.png)

## Code to Reproduce results from params that produced above plot (n_pop = 5)
```{code}
# Load packages
  library("emmeans")
  library("lme4")
  library("rlist")
  library("dplyr")
  
  # .csv sent via Slack
  model_df <- read.csv("model_df.csv")
  result <- read.csv("result.csv") # FOR REFERENCE: original parameter set and results from full code
  
  # Anova
    test_temp <- lm(phen_corrected ~ exp_env_factor * gen_factor, data = model_df)
    
    # Estimated Marginal Means
    emm_options(msg.interaction = FALSE)
    emm_E = as.data.frame(emmeans(test_temp,"exp_env_factor"))
    emm_G = as.data.frame(emmeans(test_temp, "gen_factor"))
    emm_GxE = as.data.frame(emmeans(test_temp, ~ exp_env_factor*gen_factor))
    
    # Magnitude of GxE -- EMMs
    GxE_emm <- abs(emm_GxE$emmean[emm_GxE$gen_factor == "G_1" & emm_GxE$exp_env_factor == "E_1"] - # GxE (Phenotype of ith genotype in jth environment)
                  emm_G$emmean[emm_G$gen_factor == "G_1"] - # phenotype of ith Genotype
                  emm_E$emmean[emm_E$exp_env_factor == "E_1"] + # phenotype of jth Environment
                  mean(model_df$phen_corrected)) # Overall mean phenotype # GxE_emm = 0.05991497 
    
    # Magnitude of GxE -- Omega^2
    w2_GxE = (summary(aov(test_temp))[[1]][3,2] - #(SS_effect -
             (summary(aov(test_temp))[[1]][3,1]*summary(aov(test_temp))[[1]][4,3])) / #(Df_effect * MS_error))/
             (sum(summary(aov(test_temp))[[1]][,2]) + # (SS_total+
             (summary(aov(test_temp))[[1]][4,3])) # MS_error) #  = 0.807962
    
    # Magnitude of GxE -- Eta^2
    eta2_GxE = summary(aov(test_temp))[[1]][3,2]/sum(summary(aov(test_temp))[[1]][,2]) # = 0.8260634
 
```
