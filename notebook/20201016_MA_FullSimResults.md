# Sim Results

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

~~KEL: change to "Non-significant GxE" and "Significant GxE"~~

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Simulation_10.10.2020/PhenotypePanel.png)


## Figure 4
Confusion Panel to show full extent of differences between full reciprocal transplant (RT) and paired common garden (CG)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Simulation_10.10.2020/Confusion2Panel.png)


## Figure 5

~~I can't read the numbers on the dark blue background~~
Heatmap: Top number = total samples; bottom number = Power
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Simulation_10.10.2020/Heatmap_Panel_HighStdDev.png)


## Figure 6

~~Both of these have the same COV-Ge and delta-GxE in the text~~

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Simulation_10.10.2020/Real_DataPanel.png)


## Figure 7
~~What is the filtering on this plot? Did you exclude cases in which one or both equaled zero?~~

~~Add the nice curvy line as a 2nd panel~~

Filtering now set to remove all 0's, false positives, and keep only data that are significant CovGE or GxE. 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Simulation_10.10.2020/CovGxE_Tradeoff.png)


# Supplementary Materials:

## Confusion Matrix 

### Full reciprocal transplant Design (n = 10000) 

**Raw Data**

| | CovGE Permutation | CovGE Bootstrap | GxE Permutation | GxE Bootstrap |
| ---| ---| ---| ---| ---|
| True Positive | 821, 8.21% | 6226, 62.26% | 4330, 43.3%  | 9479, 94.79% |
| True Negative | 71, 0.71% | 62, 0.62% | 74, 0.74% | 20, 0.2% |
| False Positive | 0, 0% | 9, 0.09% | 0, 0% | 54, 0.54% |
| False Negative | 9108, 91.08% | 3703, 37.03%| 5596, 56.96% | 447, 4.47%|
|---|---|---|---|---|
| False Negative Rate | 0.917 | 0.373| 0.564 | 0.045|
| False Positive Rate | 0 | 0.127 | 0 | 0.730 |

**Means Data**

| | CovGE Permutation | CovGE Bootstrap | GxE Permutation | GxE Bootstrap |
| ---| ---| ---| ---| ---|
| True Positive |709, 7.09%|6054, 60.54% | 7006, 70.06%| 9496, 94.96%|
| True Negative |71, 0.71%|66, 0.66%| 68, 0.68%| 21, 0.21%|
| False Positive |0|5, 0.05% |6, 0.06%| 53, 0.53%|
| False Negative |9220, 92.28%|3875, 38.75% |2920, 29.2% | 430, 4.3%|
|---|---|---|---|---|
| False Negative Rate |0.929|0.390|0.294|0.043|
| False Positive Rate |0|0.070|0.081|0.716|

### Paired Common Garden Design (n = 10000) 

**Raw Data**

| | CovGE Permutation | CovGE Bootstrap | GxE Permutation | GxE Bootstrap |
| ---| ---| ---| ---| ---|
| True Positive |621,6.21% |5850, 58.5% |4460,44.6% |9436, 94.36% |
| True Negative |87, 0.87% |79, 0.79% | 89, 0.89%|24, 0.24% |
| False Positive | 0,0% | 8, 0.08%|0,0% |65, 0.65% |
| False Negative |9292, 92.29%| 4063, 40.63%| 5451, 54.51%|475, 4.75% |
|---|---|---|---|---|
| False Negative Rate |0.937|0.41 |0.55|0.048|
| False Positive Rate |0 |0.092|0 |0.73 |  

**Means Data**

| | CovGE Permutation | CovGE Bootstrap | GxE Permutation | GxE Bootstrap |
| ---| ---| ---| ---| ---|
| True Positive |526, 5.26% |5728, 57.28% |7070, 70.7% |9469, 94.69% |
| True Negative | 87, 0.87%| 82, 0.82%| 81, 0.81%|32, 0.32% |
| False Positive |0, 0% |5, 0.05% |8, 0.08% |57, 0.57% |
| False Negative | 9387, 93.87%| 4185, 41.85%| 2841, 28.41%| 442, 4.42%|
|---|---|---|---|---|
| False Negative Rate |0.947 | 0.422| 0.287| 0.045|
| False Positive Rate |0|0.057|0.09 | 0.640|

## Confusion Plots for Means Data
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Simulation_10.10.2020/Confusion2Panel_means.png)

## Parameter Coverage
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Simulation_10.10.2020/HexPlotPanel.png)

## CovGE vs. GxE
~~The setup here is confusing, because the actual number of total samples is very different for the top and bottom graphs. Might be easier to chat about.~~

![image]()

## Heatmap for Std. Dev = 0.5
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Simulation_10.10.2020/Heatmap_Panel_LowStdDev.png)

## Omega^2 vs. GxE with Estimated Marginal Means
~~Remove the 1:1 line. I don't think there is a reason to expect these should be related 1:1~~ 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Simulation_10.10.2020/OmegaVsEMM_10.20.2020.png)

## Confusion Panels - Raw Data
~~Add a legend.~~ Also where are the false positives for bootstrap? They're in the confusion matrix. I think there's a bug in your code where the FP and FN are based on the same thing, and not the criteria in panel of the plot.

**Full Reciprocal Transplant**
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Simulation_10.10.2020/16panel_confusion_Scen1_raw.png)

**Paired Common Garden**
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Simulation_10.10.2020/16panel_confusion_Scen2_raw.png)

## Confusion Panels - Mean Data
~~Add a legend~~
Something weird is going on here. Why aren't the CI centered on the estimates for the CovGE?

**Full Reciprocal Transplant**
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Simulation_10.10.2020/16panel_confusion_Scen1_means.png)

**Paired Common Garden**
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Simulation_10.10.2020/16panel_confusion_Scen2_means.png)


## Covariance: Means Vs. Raw Checks
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Simulation_10.10.2020/Cov_MeansVsRaw_panel.png)

## GxE: Means Vs. Raw Checks
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Simulation_10.10.2020/GxE_MeansVsRaw_panel.png)







