####Mo' Problemz. 

When I plot the magnitude of GxE using emms vs omega^2, here is the result. Purple = Omega^2 is significant, grey = not signficant. 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/5.25.OmegaVsEmm.png)

If I filter this data to just Omega^2 estimates above 0.9 and EMM estimates below 0.1, I get 26 hits. The parameters that led to that outcome are: 

delta_env |delta_gen| sample_size |n_pop |std_dev |interaction
|---|---|---|---|---|---|
       0.01   |   0.01   |       10  |  15   |  1.0  |      7.50
       0.01  |    0.01    |      20   | 15  |   0.5  |     11.25
       0.01   |   0.01      |     5 |   15  |   0.5    |   11.25
       0.50   |   0.01   |       20  |  15   |  0.5   |    15.00
       0.01    |  0.01     |      5  |  15   |  1.0   |    15.00
       0.01     | 0.01     |      5  |  15   |  0.5   |    11.25
       0.50     | 0.01     |     10  |  15   |  1.0   |    15.00
       0.01     | 0.01     |     10  |  15   |  0.5   |    11.25
       0.01     | 0.01     |     10  |  15   |  0.5   |    11.25
     0.01  |    0.01       |   20   | 15   |  0.5     |  11.25
     0.01   |   0.01       |   20   | 15   |  0.5     |   7.50
     0.01    |  0.01       |   20   | 15   |  0.5     |  15.00
     0.01  |    0.01       |    5   | 15   |  0.5     |  15.00
     0.01   |   0.01       |    5   | 15   |  1.0     |  15.00
     0.01   |   0.01       |   10   | 15   |  0.5     |  15.00
     0.50   |   0.01       |   20   | 15   |  0.5     |  15.00
      0.01  |    0.01      |     5  |  15  |   1.0    |   15.00
     0.01   |   0.01       |   10   | 10   |  0.5     |   7.50
      1.00  |    0.01      |     5  |  15  |   1.0    |   15.00
     0.01    |  0.01       |   20   | 15   |  1.0     |  11.25
      0.50   |   0.01      |    20  |  15  |   1.0    |   11.25
      0.01   |   0.01      |    10  |  15  |   0.5    |    7.50
      0.01   |   0.01      |     5  |  15  |   0.5    |   11.25
      0.01   |   0.01      |     5  |  15  |   1.0    |   11.25
     0.01     | 0.01       |    5   | 15   |  1.0     |  15.00
     0.01     | 0.01       |   20   | 15   |  1.0     |  15.00


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
In the above, the above code with the output listed above produces w2_GxE = 0.007614049
