# Magnitude of GxE Interactions

Genotype by environment interactions are common in nature and therefore an important component of CoGV and CnGV. 
However, finding an effective means of measuring the magnitude of GxE is not straightforward. 

One method involves the complementary metrics, `omega^2` and `eta^2`

`eta^2` is the most simple calculation: 
`SSeffect / SStotal `

`omega^2` is less subject to bias of sample size than `eta^2` and is calculated: 
`SSeffect - (dfeffect* MSerror)/MSerror + SStotal`

Each of these metrics gives an estimate between zero and one. A value of one indicates the GxE interaction explains 100% the variance in the model.

Another method is estimated marginal means (emm). Emms estimate effect sizes based on models. Emms may be preferable since it should incorporate effect size into the calculation rather than simply proportion variance explained. 

There are two ways to calculate GxE using emms, manually vs. function `emmeans` in R.

To calculate manually using simulated data: 
(this is based on Katie's design found [here](https://github.com/RCN-ECS/CnGV/edit/master/notebook/20191115_KEL_compareOmega2_effectsize_GxE.md): 

Only one calculation is shown below, but all calculations should be equal regardless of comparison (G11, G12, G21, G22). 

```R
overall_mean <- mean(ind_dat$phen_corrected) 

G1M = abs(overall_mean -
         (mean(ind_dat$phen_corrected[ind_dat$gen == "G1"])) - # G1
         (mean(ind_dat$phen_corrected[ind_dat$env == "E1"])) +  # E1
         (mean(ind_dat$phen_corrected[ind_dat$env == "E1" & ind_dat$gen == "G1"]))) 
```

To calculate using the `emmeans` function in R: 

**Important Note on `emmeans` function:** All predictors in model must be factors to be treated properly. If the predictor is a numeric, the function will take the average value. 

Again, only one is shown, but should be same for all comparisons

```
emm_E = emmeans(test_temp,"env") # test_temp is the model
emm_G = emmeans(test_temp, "gen")
emm_GxE = as.data.frame(emmeans(test_temp, ~ env*gen))

G1E1_emm = abs(overall_mean -
                       (emm_G$emmean[emm_G$gen=="G1"])- # G1
                       (emm_E$emmean[emm_E$env=="E1"])+ # E1
                       (emm_GxE[1,3])) # G1E1
```

## Test of Emmeans vs. `eta^2` and `omega^2`

I'll compare the performance of each using simulated data with 4 scenarios:

Intercepts of G1 and G2 = 0
Standard Deviation = 0

|Scenario| slope G1 | slope G2 | environment type | `eta^2` | `omega^2` | Emmeans Manual | Emmeans R |
|----------|----------|----------|----------|------------- | ------------ | ------------- | ------------- |
|1 | 1 | -1 | continuous (-2, 2)|  1 | 1 | 0.9746794 | 0.9746794 |
|2 | 1 | -1 | categorical (E1, E2)| 1 | 1 | 0.9746794 | 0.9746794 |
|3 | 0.25 | -0.25 | continuous (-2, 2)| 1 | 1 | 0.9746794 | 0.9746794 |
|4 | 0.25 | -0.25 | categorical (E1, E2)|1 | 1 | 0.9746794 | 0.9746794 |

The emmeans give slightly lower estimates, but there is no difference between the different slopes. 

I'm going to repeat the same steps with some error (Std. dev = 0.05) and bootstrap these estimates 100x to get means and CIs: 

|Scenario | `eta^2` (CI) | `omega^2` (CI) | Emmeans Manual (CI) | Emmeans R (CI)|
| ------------ | ------------- | ------------ | ------------- | ------------- |
| 1 |  0.99 (0.9990-0.9997) | 0.99 (0.9989-0.9997) | 0.97 (0.9742-0.9745) | 0.97 (0.9742-0.9745)|
| 2 | 0.99 (0.9988-0.9997) | 0.99 (0.9988-0.999) | 0.97 (0.9741-0.9745) | 0.97 (0.9741-0.9745)|
| 3 | 0.99 (0.9845-0.9957) | 0.98 (0.9827-0.9953) | 0.97 (0.9671-0.9726) | 0.97 (0.9671-0.9726)|
| 4 | 0.99 (0.985-0.995) | 0.99 (0.9840-0.9949) | 0.97 (0.9677-0.9724) | 0.97 (0.9677-0.9724)|


#### Conclusions:
The prediction that emms would incorporate effect size more than `omega^2` and `eta^2` is not supported by these initial simulations. All methods seems to perform similarly because the data are all standardized ((data - mean)/std. dev). 

Moving forward we will use `emmeans` function in R to estimate the magnitude of GxE.

## Simulating Emmeans with Covariance 

Cogradient variation (CoGV) gives positive covariance, while countergradient variation (CnGV) gives negative covariance. As the magnitude of GxE goes up, the covariance among genotypes should decrease. 

### Proof of Concept: 
```# Categorical Starting parameters
Diff_means_cat <- list(
  "data_type" = c("categorical"), 
  "intercept_G1" = seq(from = -5, to = 5, by = 2),
  "slope_G1" = seq(from = -1, to = 1, by = 0.5),
  "intercept_G2" = seq(from = -5, to = 5, by = 2),
  "slope_G2" = seq(from = -1, to = 1, by = 0.5), 
  "sd" = 0.5, #seq(from = 0, to = 1, by = 0.5),
  "sample_size" = c(5)) 
```
Data are simulated (`data_generation`) and processed through the appropriate covariance generating function (`Cov_matrix_sim_cat`). 

Plotting the results gives the following:
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/GxE_Emeans_Covariance.png)
