# Updated Plot designs

First of all, counts so we can identify good binning windows for False Pos and False Neg Rates. Looks like binning the upper window at 0.75-1.0 gives good sampling counts. 

| --- | Full Reciprocal Transplant: CovGE | Full Reciprocal Transplant: GxE | Common Garden: CovGE | Common Garden: GxE |
|---|---|---|---|---|
| Population (true) Value = 1 exactly | 2 | 0 | 0 | 0 | 
| between 0.9 and 1.0| 697| 0|110|1|
| between 0.8 and 1.0|1521| 357| 477| 195 | 
| between 0.75 and 1.0 |2005 | 674| 787| 571| 
| between 0.7 and 1.0|2491 |1285| 1192| 1413| 
| Population Value = 0 exactly |2058 |2056 | 3777|3610 | 


After discussing various different ways to visualize the data, I have come up with the following plots: 

**Molly: Make sure all axes labeled with SAMPLE or POPULATION**

### Figure 3.
Phenotype plot with just 2x2 and 4x4s. Bigger text. 
![image]()

### Figure 4. 
GxE Panel. Three plots. Show omega^2 (how these metrics differ), FPR for Anova vs. Permutation, and Powerheatmap. 

### Figure 5.
Cov Panel. FalsePositive vs. Power... Pick 5-6 designs. Then use Line plots to show FPRs according to effect size, FPR according to designs, and heatmap to show. 

### Figure 6. 
Geoff and my data. 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/12.21.RealData.panel.png)

### Figure 7. 
Tradeoff. (candidate for supplement if we want to cut to 6 plots?) 
![image]()

## Supplemental Plots

### Supplemental Figure 1. 
hex plot
![image]()

### Supplemental Figure 2. 
Same panel for CovGE but with Means? (or do 1:1 line test on Power to see if they differ) 

### Supplemental Figure 3. 
Variance Partitioning - Katie has these plots

### Supplemental Figure 4. 
Raw vs. Means - update to lwr/upr ci and change axes. 

### Supplemental Figure 5. 
PL statistic Panels - Katie has these plots.


