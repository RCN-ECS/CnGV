# Sim Results -- Updated 10/30/2020

I ran a full round (10 replicates) for the following starting parameters. [Recall I am using a different way to set up parameters.](https://github.com/RCN-ECS/CnGV/blob/master/notebook/20201006_NewParamGeneration.md)

I have made the following plots, ordered in terms of where they will go in the paper or supplementary material. Figures 1 and 2 are handmade heuristics. They can be found [here.](https://docs.google.com/document/d/1CcoJFTX6I8zPptLAzBvITiTaTkNtdgGxNNdQlrcxlSY/edit#heading=h.68abrs8dfl95)
```{params}
param_list <- list( 
  reps = c(10), 
  delta_env = NULL, 
  delta_gen = NULL,
  sample_size = c(2,4,8,16), 
  n_pop = c(2,4,8,16),
  env_scenario = c(1,2),  # 1 = npop = n_env; 2 = multiple pops across 2 envs
  std_dev= c(.5, 1), # Scale
  interaction = NULL) 
```

Note that some of these images are not great quality because they are saved as PNGs from powerpoint. I have high quality pdfs that will be submitted for review. Figure legends are written in the document linked to above. 

## Figure 3

May need to change if we want to add variance partitioning, I might also make the text bigger.

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_10.30.2020/PhenPanels.png)


## Figure 4
Confusion Panel to show full extent of differences between full reciprocal transplant (RT) and paired common garden (CG)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_10.30.2020/cogvcngv_plot.png)


## Figure 5

Heatmap: Top number = total samples; bottom number = Power
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_10.30.2020/Heatmap_Panel_HighStdDev.png)


## Figure 6

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_10.30.2020/MollyGeoff_DataPanel.png)


## Figure 7

Filtered to remove all 0's, false positives, and keep only data that are significant CovGE or GxE. 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_10.30.2020/10.30.CovGxE_Tradeoff.png)


# Supplementary Materials:

## Confusion Matrix 

### Full reciprocal transplant Design (n = 10000) 

**Raw Data**

| | CovGE Permutation | CovGE Bootstrap | GxE Bootstrap | GxE Permutation |
| ---| ---| ---| ---| ---|
| True Positive | 560 - 5.6% | 5595 - 55.95% |  9496 - 94.96% | 4730 - 47.3%|
| True Negative | 143 - 1.43% | 127 - 1.27%|  16 - 0.16%| 65 - 0.65%|
| False Positive | 0 |  16 - 0.16%| 49 - 0.49%| 0 | 
| False Negative | 9297 - 92.97%| 4262 - 42.62% |439 - 4.39% | 5205 - 52.05% | 
|---|---|---|---|---|
| False Negative Rate | 0.943 | 0.432 |  0.044 |0.524 |
| False Positive Rate | 0| 0.112 | 0.754 | 0 | 

**Means Data**

| | CovGE Permutation | CovGE Bootstrap | GxE Bootstrap | GxE Permutation |
| ---| ---| ---| ---| ---|
| True Positive | 485 - 4.85%| 5463 - 54.63%| 9515 - 95.15% |7264 - 72.64% |
| True Negative |143 - 1.43%| 133 - 1.33%| 19 - 0.19%| 62 - 0.62%|
| False Positive |0| 10 - 0.1%| 46 - 0.46% | 3 - 0.03%|
| False Negative |9372 - 93.72%|4394 - 43.94% | 420 - 4.2% |2671 - 26.71% |
|---|---|---|---|---|
| False Negative Rate |0.95| 0.445|0.042 |0.269 |
| False Positive Rate |0| 0.069|0.708 |0.046 |

### Paired Common Garden Design (n = 11000) 

**Raw Data**

| | CovGE Permutation | CovGE Bootstrap | GxE Permutation | GxE Bootstrap |
| ---| ---| ---| ---| ---|
| True Positive |706 - 6.42% | 4932 - 44.8%| 6583 - 59.8%| 10982 - 99.8%|
| True Negative |779 - 7.08% | 757 - 6.88%| 18 - 0.164%| 0 |
| False Positive |1 - 0.009% |23 - 0.209% | 0|18 - 0.164% |
| False Negative |9514 - 86.5%| 5288 - 48.1%| 4399 - 40.0% |0 |
|---|---|---|---|---|
| False Negative Rate |0.93|0.517|0.40| 0 |
| False Positive Rate |0.001|0.029 | 0| 1.0 |  

**Means Data**

| | CovGE Permutation | CovGE Bootstrap | GxE Permutation | GxE Bootstrap |
| ---| ---| ---| ---| ---|
| True Positive | 637 - 5.79%|4798 - 49.3% |8873 - 80.7% | 10982 - 99.8% |
| True Negative | 778 - 7.07%| 765 - 6.95%| 18 - 0.164%| 0|
| False Positive |1 - 0.009% | 14 - 0.127%| 0 | 18 - 0.164%|
| False Negative |9584 - 87.1% | 5423 - 49.3%| 2109 - 19.2%| 0|
|---|---|---|---|---|
| False Negative Rate |0.938 |0.531 |0.192 | 0|
| False Positive Rate |0.001 | 0.018| 0 | 1.0|

## Confusion Plots for Means Data

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_10.30.2020/ConfusionPanel_means.png)

## Parameter Coverage

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_10.30.2020/10.30.HexPlot_supplement.png)

## CovGE vs. GxE

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_10.30.2020/CovGxE_Pyramid.png)

## Heatmap for Std. Dev = 0.5

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_10.30.2020/Heatmap_Panel_lowStdDev.png)

## Omega^2 vs. GxE with Estimated Marginal Means

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_10.30.2020/10.30.EmmVsAnova.png)

## Confusion Panels - Raw Data

**Full Reciprocal Transplant**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_10.30.2020/Supp_confmat_Raw1.png)

**Paired Common Garden**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_10.30.2020/Supp_confmat_Raw2.png)

## Confusion Panels - Mean Data

**Full Reciprocal Transplant**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_10.30.2020/Supp_confmat_Means1.png)

**Paired Common Garden**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_10.30.2020/Supp_confmat_Means2.png)

## Covariance: Means Vs. Raw Checks - Reciprocal Transplant 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_10.30.2020/11.1.Cov_RawVsMeans.png)

## Covariance: Means Vs. Raw Checks - Paired Common Garden

Yes I see it, No I haven't figured it out. 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_10.30.2020/11.1.Cov_RawVsMeans_Env2.png)

## GxE: Means Vs. Raw Checks - Reciprocal Transplant

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_10.30.2020/11.1.GxE_RawVs.Means.png)

## GxE: Means Vs. Raw Checks - Common Garden

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_10.30.2020/11.1.GxE_MeansVsRaw_Scen2.png)

## GxE Permutation P-values vs. Anova P-values

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_10.30.2020/11.1.GxE_PvalComparison.png)







