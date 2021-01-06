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
In the meeting we discussed a few different ways of plotting false positives. If we want to show 1) Here's how the metric works, 2) Here's how accurate it is, nad 3) Here's the power. 

In light of that, I was advised to pull a few scenarios and plot the false positive rates.  I've created a few different ways here. In each plot, panel A, C, and D are the same. A shows the relationship between omega^2 and DeltaGxE colored according to significance (goal 1). C and D show power (1-FNR) for Full Reciprocal Tranplant and Common Garden designs (goal 3). The plot that changes out is B. We need to find the best way to show the accuracy of the methods.

Permutation FPR is 0 in basically every design. The overall FPR for permutation across all sample sizes and # genotypes is 0.01, so its not totally surprising. The only design with FPR of 0.03 is 4 sample size/ 8 genotype for FRT which is only shown in the 2nd & 4th mockup.  

The first mock-up shows false positive rates for just  5 scenarios with 128 total samples (3 for Common Garden, 2 for Recip. Transplant). 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/1.6.GxEPanels1.png)

The second mock-up shows false positive rates for the 9 scenarios with power above 80% ordered by total sample size.
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/1.6.GxEPanels2.png)

The next two show the same data but now I've made the designs into the bars and grouped according to method (ANOVA vs. Permutation)

The third mockup shows the flipped grouping with just 128 total samples (same as the first one above)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/1.6.GxEPanels3.png)

The final mockup shows the flipped grouping with the 9 high powered designs (same as the second one above)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/1.6.GxEPanels4.png)

### Figure 4.
Same as GxE, I've made some different mockups to figure out how to best show the Covariance data. Cov Panel. These only show False Positives and Power. The bottom plots remain the same but the top plots (A) differ. 

The first mock up shows False Positive rates in just 128 total samples
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/1.6.CovariancePanels1.png)

The second one shows False Positive rates in the 4 scenarios with power above 80%. 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/1.6.CovariancePanels2.png)

The third is the same as the first but flipping to show design as the bars and method as the x axis group 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/1.6.CovariancePanels3.png)

The final is showing the false positives for higher power designs but again with the flipping.
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/1.6.CovariancePanels4.png)

### Figure 5.
Phenotype plot with just 2x2 and 4x4s. Bigger text. 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/1.4.21.PhenotypePlotPanel.png)

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


