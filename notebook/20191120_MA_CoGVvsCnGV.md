# Predicting CoGV and CnGV based on G1E1 and G2E2

## TL:DR - Just use covariance estimate to define whether a scenario is CoGV or CnGV

Cogradient variation can be defined as positive covariance, while countergradient variation can be defined as negative covariance. Defining which scenario requires knowledge of each genotype's native environment. In our simulations, G1 is always native to E1, while G2 is always native to E2. 

The following two images are examples of *CoGV*. Native environment is denoted by the X. If slope is positive, then the intercept of G2 should be greater than G1. If the slope is negative, the intercept of G1 should be greater than G2. 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/plotA.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/PlotC.png)

Conversely, in these next two examples of *CnGV*, if the slope is positive, then the intercept of G1 should be greater than G2. If the slope is negative, the intercept of G2 should be greater than G1. 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/PlotB.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/plotD.png)

In pure GxE scenarios, we should expect the slopes of each genotype to be exactly opposite but have the same intercept: 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/PlotE.png)

I coded these predictions in the simulations using the following statement in the `data_generation` function: 

```#Code
if((intercept_G2[e] == -intercept_G1[b]) & (slope_G1[d] == -slope_G2[f])){
                  type = "pure_GxE"
                }else if((slope_G1[d]<0) & (intercept_G1[b]<intercept_G2[e])){
                  type = "cngv"
                }else if((slope_G1[d]>0) & (intercept_G1[b]<intercept_G2[e])){
                  type = "cogv"
                }else if((slope_G1[d]>0) & (intercept_G1[b]>intercept_G2[e])){
                  type = "cngv"
                }else if((slope_G1[d]<0) & (intercept_G1[b]>intercept_G2[e])){
                  type = "cogv"
                }else{
                  type = "GxE"
                }
```
The "GxE" label in this initial statement is the catch-all bin for the scenarios in which the intercepts are equal or slopes are both positive and negative. 

This approach does a marginal job at predicting CoGV and CnGV but there are still some positive CnGVs and negative CoGVs and it does not define GxE scenarios well.

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/GxE_better.png)

A more general way to think about this that better includes GxE interactions is by describing the relative location of each genotype's native phenotype on a plot. To allow for GxE, the location of G2E1 can vary and so I leave that square blank.

H = High, L = Low

**CoGradient Variation:**

-- | E1 | E2
---|---|---
G1 | L | H
G2 | - | H

-- | E1 | E2
---|---|---
G1 | H | L
G2 | - | L

**CounterGradient Variation:** 

-- | E1 | E2
---|---|---
G1 | L | H
G2 | - | L

-- | E1 | E2
---|---|---
G1 | H | L
G2 | - | H

These show that you can predict CoGV or CnGV based on the spatial relationship of G1E1 to G2E2. If the phenotypes are similar, it is CnGV. If they are not similar, it is likely CoGV. This is not new, it is just a different way of describing the above plots that Molly typically uses. However this approach can become problematic when phenotypes do not perfectly match. 

This begs the question: 
**When phenotypes from G1E1 and G2E2 do not perfectly match, how different must the native phenotypes (G1E1 and G2E2) be in order to flip from negative to positive covariance?** 

Answering this question may allow us to predict when significant covariance is occurring in nature. Additionally, since the magnitude of GxE is also related to the strength of covariance (see the triangle plot above), at what point does a GxE interaction weaken the ability to detect covariance? 

I used simulated data based on the following parameters: 
```#starters
Diff_means_cat <- list(
  "data_type" = c("categorical"), 
  "intercept_G1" = seq(from = -5, to = 5, by = 2),
  "slope_G1" = seq(from = -1, to = 1, by = 0.5),
  "intercept_G2" = seq(from = -5, to = 5, by = 2),
  "slope_G2" = seq(from = -1, to = 1, by = 0.5), 
  "sd" = 0.05, #seq(from = 0, to = 1, by = 0.5),
  "sample_size" = c(5)) 
```
Here is the resultant plot: 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/cogvcngv_plot.png)

One small tidbit gleaned from this day's worth of wasted time (see TL:DR above) is that scenarios with a perfect GxE are the switchpoints between CoGV and CnGV. But "perfect" GxEs are probably non-existent in nature.
Carry on.
