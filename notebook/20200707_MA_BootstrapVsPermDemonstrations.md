# Bootstrap vs. Permutation Comparisons.

The previous [post](https://github.com/RCN-ECS/CnGV/blob/master/notebook/20200702_MAA_Discouraged.md) helped identify differences in bootstrap vs. permutation approaches for determining significance. 

### Raw Data Confusion Matrices :)

I have updated these plots (removed rounding and reordered panels). *Green = True value, Red = Estimate* 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.7.BootstrapCovariance.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.7.CovariancePermutation.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.7.GxEBootstrap.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.7.PermutationGxE.png)

For covariance, we found that permutation is an exceedingly conservative approach resulting in many false negatives. Therefore bootstrap may be a more reliable method for determining significance (assuming CovGE != 0 and the 95% CIs don't overlap with zero). 

For GxE, bootstrap randomization creates levels of GxE that may not accurately represent significance of the GxE emm. Therefore for GxE, permutation may be the better approach. 

### Where are false positives/negatives for the bootstrap? 
If we use Bootstrap for the covariance, we want to know what levels of covariance are giving us false positives and negatives. Here is a barchart showing the proportion of instances for the absolute value of binned covariance.

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.8.CovarianceBootstrapBarPlot.png)


### Means Data Confusion Matrices :) 

Thankfully, these follow the same patterns as the raw data, give or take a few. 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.7.MeansCovarianceBootrap.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.7.MeansCovariancePermutation.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.7.MeansGxEBootstrap.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.7.MeansGxEPermutation.png)

## Why not just use ANOVA for GxE hypothesis testing?

Some reviewers may question the decision to use permutation rather than taking the p-value from the ANOVA. The reasons for using permutation are 2-fold: 
1. We want to be able to accurately measure covariance with means data, which does not use ANOVA. 
2. ANOVA p-values are based on the proportion of variance explained by fixed effects (in this case, the interaction). However, as seen in the graph below, GxE EMM effect sizes do not directly scale with the proportion of variance explained. ANOVA produces a lower true GxE magnitude than EMM until high GxE values (~0.7). (lines are b-spline curves)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.8.TrueGxEv.Omega.png)

*Below shows the original with GxE Estimate, not the true value. Retained for comparison*
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.7.GxEEmmVsAnova.png)


