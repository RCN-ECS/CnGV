# Predicting CoGV and CnGV based on G1E1 and G2E2

## TL:DR- Molly discovers the point of power analyses again 

Cogradient variation can be defined as positive covariance, while countergradient variation can be defined as negative covariance. Defining which scenario requires knowledge of each genotype's native environment. In our simulations, G1 is always native to E1, while G2 is always native to E2. 

The following two images are examples of *CoGV*. Native environment is denoted by the X. If slope is positive, then the intercept of G2 should be greater than G1. If the slope is negative, the intercept of G1 should be greater than G2. 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/plotA.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/PlotC.png)

Conversely, in these next two examples of *CnGV*, if the slope is positive, then the intercept of G1 should be greater than G2. If the slope is negative, the intercept of G2 should be greater than G1. 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/PlotB.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/plotD.png)

In pure GxE scenarios, we should expect the slopes of each genotype to be exactly opposite but have the same intercept: 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/PlotE.png)

I have coded these predictions in the simulations using the following statement in the `data_generation` function: 

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

A more general way to think about this that better includes GxE interactions is by describing the relative location of each genotype's phenotype on a plot. To allow G2's slope to vary, I leave it blank.
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

These show that you can predict CoGV or CnGV based on the spatial relationship of G1E1 to G2E2. This is not new, it is just a different way of describing the above plots that have been Molly's typical approach. Thus, if phenotypes in G1E1 equal phenotypes in G2E2, then we have CnGV. However this can become difficult when phenotypes do not perfectly match. 

This begs the question: 
**When phenotypes from G1E1 and G2E2 do not perfectly match, how different must the phenotypes (G1E1 vs. G2E2) be in order to flip from negative to positive covariance?** 

Answering this question may allow us to predict when significant covariance is occurring in nature. Additionally, since the magnitude of GxE is also related to the strength of covariance (see the triangle plot above), at what point does a GxE interaction weaken the ability to detect covariance? To answer these questions, we can use simulations and power analysis (7 months later she gets it). <pause while I run some analyses>
