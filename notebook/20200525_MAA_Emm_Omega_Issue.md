####Mo' Problemz. 

When I plot the magnitude of GxE using emms vs omega^2, here is the result. Purple = Omega^2 is significant, grey = not signficant. 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/Omega%20Vs.%20Emm.png)

Here is the code used to estimate both - based on model: 

```{m1}
    test_temp <- lm(phen_corrected ~ exp_env_factor * gen_factor, data = model_df)

```

Estimated marginal means: 
```{emm}
# Magnitude of GxE -- EMMs
    GxE_emm <- abs(mean(model_df$phen_corrected) - # Overall mean
                  (emm_G$emmean[emm_G$gen_factor == "G_1"])- # G
                  (emm_E$emmean[emm_E$exp_env_factor == "E_1"])+ # E
                  (emm_GxE[1,3])) # GxE
```

```{omega}
# Magnitude of GxE -- Omega^2
    w2_GxE = (summary(aov(test_temp))[[1]][3,2] - #(SS_effect -
             (summary(aov(test_temp))[[1]][3,1]*summary(aov(test_temp))[[1]][4,3])) / #(Df_effect * MS_error))/
             (sum(summary(aov(test_temp))[[1]][,2]) + # (SS_total+
             (summary(aov(test_temp))[[1]][4,3])) # MS_error)
```
