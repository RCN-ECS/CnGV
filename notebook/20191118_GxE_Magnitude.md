Genotype by environment interactions are common in nature and therefore an important component of CoGV and CnGV. 
However, finding an effective means of measuring the magnitude of GxE is not straightforward. 

One method involves the complementary metrics, \\omega^{2} and \\eta^{2}.

`\eta^{2}` is the most simple calcuation: `SS_{effect} / SS_{total} `

`\omega^{2}` is less subject to bias of sample size than \eta^{2} and is calculated: SS_{effect} - (df_{effect}* MS_{error})/MS_{error} + SS_{total}

Each of these metrics gives an estimate between zero and one. A value of one indicates the GxE interaction explains 100% the variance in the model.

Another method is estimated marginal means (emm). Emms estimate effect sizes based on models. Emms may be preferable since it should incorporate effect size into the calculation rather than simply proportion variance explained. 

There are two ways to calculate GxE using emms. Manual calculation vs. function `<emmeans>` in R.

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

Again, only one is shown, but should be same for all comparisons

```
emm_E = emmeans(test_temp,"env") 
emm_G = emmeans(test_temp, "gen")
emm_GxE = as.data.frame(emmeans(test_temp, ~ env*gen))

G1E1_emm = abs(overall_mean -
                       (emm_G$emmean[emm_G$gen=="G1"])- # G1
                       (emm_E$emmean[emm_E$env=="E1"])+ # E1
                       (emm_GxE[1,3])) # G1E1
```

# Test of Emmeans vs. `\eta^{2}` and `\omega^{2}`

I'll compare the performance of each using simulated data with 4 scenarios:

Intercepts of G1 and G2 = 0
Standard Deviation = 0

1: G1 slope = 1; G2 slope = -1; environment = continuous (-2, 2) 
2: G1 slope = 1; G2 slope = -1; environment = categorical (E1, E2)
3: G1 slope = 0.25; G2 slope = -0.25; environment = continuous (-2, 2) 
4: G1 slope = 0.25; G2 slope = -0.25; environment = categorical (E1, E2)

Results: 

| Scenario | `\eta^{2}` | `\omega^{2}` | Emmeans Manual | Emmeans R |
| ------------ | ------------- | ------------ | ------------- |
| 1,-1,continuous | 1 | 1 | 0.9746794 | 0.9746794 |
| 1,-1,categorical | 1 | 1 | 0.9746794 | 0.9746794 |
| 0.25,-0.25,continuous | 1 | 1 | 0.9746794 | 0.9746794 |
| 0.25,-0.25, categorical | 1 | 1 | 0.9746794 | 0.9746794 |

They are all identical. 

I'm going to repeat the same steps with some error (Std. dev = 0.05) and bootstrap these estimates 100x to get CIs: 

|Scenario | \eta^{2} (CI) | \omega^{2} (CI) | Emmeans Manual (CI) | Emmeans R (CI)|
| ------------ | ------------- | ------------ | ------------- |
| 1,-1,continuous |  0.9994531 (0.9990557-0.9997713) | 0.9993922 (0.9989492-0.9997471) | 0.9744129 (0.9742191-0.974568) | 0.9744129 (0.9742191-0.974568)|
|1,-1,categorical | 0.9994139 (0.9988966-0.9997277) | 0.9993491 (0.9988116-0.999695) | 0.9743938 (0.9741416-0.9745467) | 0.9743938 (0.9741416-0.9745467)|
| 0.25,-0.25,continuous | 0.990741 (0.9845162-0.9957977) | 0.9897063 (0.9827389-0.9953629) | 0.9701556 (0.9671041-0.9726293) | 0.9701556 (0.9671041-0.9726293)|
| 0.25,-0.25, categorical | 0.9910729 (0.985733-0.995461) | 0.9900793 (0.9840536-0.9949183) | 0.9703182 (0.9677016-0.9724649) | 0.9703182 (0.9677016-0.9724649)|


# Conclusions: 
The prediction that emmeans would better estimate effect size is not supported by these initial simulations. Each method predicts the magnitude of the interaction relative to the other factors in the model (environment and genotype alone).
