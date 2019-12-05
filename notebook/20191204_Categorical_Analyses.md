# Categorical Analyses Progress

I have spent a few days updating functions to run on simulated data or data from the meta-analysis. These functions are all stored in the github CnGV src folder under "Categorical_analyses".

## Estimate Covariance and GxE

I have created 2 functions. One to handle raw data, and one to handle means/SE. The main difference is that raw data uses `anova` to get estimated marginal means using function `emmeans` while means use the means provided by the study to estimate covariance. 
To calculate GxE using means, I use the following code (based on method from book)

**Important: Number of genotypes must match number of experimental environments.**
```#meansGxE
GxE_emm = abs(overall_mean -
             (mean(input_df$phen_corrected[input_df$gen_factor == "G_1"]))- # G1 mean
             (mean(input_df$phen_corrected[input_df$exp_env_factor=="E_1"]))+  # E1 mean
             (mean(input_df$phen_corrected[input_df$exp_env_factor=="E_1" & input_df$gen_factor == "G_1"]))) #G1E1
```
The GxE and true covariance are extracted. 

## Bootstrap

I have 2 bootstrap functions for raw or means data as well. 

The raw data are reshuffled with replacement in standard bootstrapping format. 
The means are resampled means using the function below and assuming Central Limit Theorum.

In the below data, `Cond_E` is a dataframe filtered for each genotype and environment. `Phen_data` are the means, while `phen_mean_SE` is the standard error of the means which we employ as standard deviation. 

```#Bootstrap means
# Shuffle data 
new_phen <- rnorm(nrow(cond_E), mean = cond_E$phen_data, sd =  cond_E$phen_mean_SE)      
```

## Permutation

The permutation codes for raw data are similar to the bootstrapping codes, except that the raw data is reshuffled without replacement.

For the means data, rather than resampling, I shuffle the G's and E's without replacement and use this to form the null distribution. **Need to check to make sure this is appropriate**

## Overall
I have 2 final functions. One for simulated data and one for meta analysis data. I'd have prefered to combine them, but there are enough differences between the two data types that required 2. In particular, I have to pipe in meta-data for the meta-analysis to make sense, and I couldn't come up with a way to easily create or ignore that necessity if the data was simulated. 

These functions give me the Covariance and GxE estimates with confidence intervals and P-values. 

```#output
# Example Output frrom simulated data across 4 different slope combinations: 
Index true_covariance  Cov_lowCI Cov_highCI covariance_p_value GxE_magnitude
1     1      -0.4662494 -0.5717200 -0.3534593         0.01960784    0.01828558
2     2      -0.2601160 -0.3654361 -0.1848236         0.01960784    0.31666121
3     3      -0.2283763 -0.3171805 -0.1468530         0.01960784    0.34509544
4     4      -0.2469642 -0.2979271 -0.1799924         0.01960784    0.01164237

    GxE_lowCI GxE_highCI GxE_pvalue
1 0.001680960 0.07126563 0.07843137
2 0.258022895 0.36594640 0.03921569
3 0.306580030 0.38973143 0.00000000
4 0.001289295 0.05018947 0.07843137
```

## Next Steps: 

1. Automate the meta-analysis function so that multiple studies can be analyzed at once. 
2. May need a better reference than "index" for the simulated data
3. Extract significant GxE and Covariance estimates from both simulations and meta-analysis and plot! 
