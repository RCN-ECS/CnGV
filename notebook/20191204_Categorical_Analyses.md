# Categorical Analyses Progress

I have spent a few days updating functions to run on simulated data or data from the meta-analysis. These functions are all stored in the github CnGV src folder under "Categorical_analyses".

## Estimate Covariance and GxE

I have created 2 functions. One to handle raw data, and one to handle means/SE. The main difference is that raw data uses `anova` to get estimated marginal means using function `emmeans` while means use the means provided by the study to estimate covariance. 
To calculate GxE using just means, I use the following code (based on method from book)

```#meansGxE
GxE_emm = abs(overall_mean -
              (mean(input_df$phen_corrected[input_df$gen_factor == "G_1"]))- # G1
              (mean(input_df$phen_corrected[input_df$exp_env_factor=="E_1"]))+  # E1
              (mean(input_df$phen_corrected[input_df$exp_env_factor=="E_1" & input_df$gen_factor == "G_1"]))) #G1E1
```
