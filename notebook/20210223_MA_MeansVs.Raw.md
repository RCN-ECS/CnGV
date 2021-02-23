# Comparison of performance of means vs. raw data 

The following plots show false positive rates for CovGE and GxE for raw data and means data. These should give an idea of whether our approach to estimate significance and confidence intervals for means is doing what we think it should be doing. 
For the following plots, the legend is the same: (Yellow = 1; Purple = 0; Green/blue = 0.5)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/2.23.Heatmap_LEGEND.pdf)

## CovGE False Positive Comparison
False Positive rates according to bootstrap (upper tiles) and permutation (lower tiles) for CovGE for raw (LEFT) and means (RIGHT) data for Full reciprocal transplant (TOP) and paired common garden (BOTTOM) scenarios. 
There is good agreement for full reciprocal transplant designs, but differences in the 4 and 8 genotype plots. I'll explore further, not sure what is causing that. 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/2.23_CovGE_FP_MeansVsRaw.pdf)

## GxE False Positive Comparison
False Positive rates for GxE for raw (LEFT) and means (RIGHT) for full reciprocal transplant (TOP) and paired common garden (BOTTOM). Since we don't run ANOVAs on the means data, I only show Permutation for means data, but show ANOVA (lower tile) and permutation (upper tile)for raw data. Again, good agreement. 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/2.23_GxE_FP_MeansVsRaw.pdf)

## CovGE False Negative Comparison
Same layout as above, just now showing false negative rates. Power would be 1-value shown in the tile. These show good agreement. 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/2.23_CovGE_FN_MeansVsRaw.pdf)

## GxE False Negative Comparison
Again, same layout. Since we do not perform ANOVAs on means data, I only show FNRs from permutation. Power would be 1-FNR. 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/2.23.GxE_FN_MeansVsRaw.pdf)
