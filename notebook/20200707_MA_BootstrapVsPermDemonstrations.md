# Bootstrap vs. Permutation Comparisons.

The previous [post](https://github.com/RCN-ECS/CnGV/blob/master/notebook/20200702_MAA_Discouraged.md) helped identify differences in bootstrap vs. permutation approaches for determining significance. 

I have updated these plots (removed rounding and reordered panels). *Green = True value, Red = Estimate* 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.7.BootstrapCovariance.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.7.CovariancePermutation.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.7.GxEBootstrap.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.7.PermutationGxE.png)

For covariance, we found that permutation is an exceedingly conservative approach resulting in many false negatives. Therefore bootstrap may be a more reliable method for determining significance (assuming CovGE != 0 and the 95% CIs don't overlap with zero). 

For GxE, bootstrap randomization creates levels of GxE that may not accurately represent significance of the GxE emm. Therefore for GxE, permutation may be the better approach. 

Some reviewers may question why we use permutation rather than taking the p-value from the ANOVA. The reasons are 2-fold. 
1. We want to be able to accurately measure covariance with means data, which does not use ANOVA. 
2. ANOVA p-values are based on the proportion of variance explained by fixed effects (in this case, the interaction). Therefore, it is conceivable to have different magnitudes of GxE based on estimated marginal means, but they explain same proportions of variation (same omega^2). 


