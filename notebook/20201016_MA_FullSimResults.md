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

KEL: change to "Non-significant GxE" and "Significant GxE"

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Simulation_10.10.2020/PhenotypePanel.png)


## Figure 4
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Simulation_10.10.2020/Confusion2Panel.png)


## Figure 5

I can't read the numbers on the dark blue background

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Simulation_10.10.2020/HeatMaps.png)


## Figure 6
Both of these have the same COV-Ge and delta-GxE in the text
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Simulation_10.10.2020/Real_DataPanel.png)


## Figure 7
What is the filtering on this plot? Did you exclude cases in which one or both equaled zero?
Add the nice curvy line as a 2nd panel. 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Simulation_10.10.2020/GxE_Cov_Tradeoff.png)


# Supplementary Materials:

## Confusion Matrix 
These percents and rates are the same for the means data. 

KEL note: I think there is a mistake if this is the same for the means data.  The CI for the means data is not the same, so the percentages should be different. Also, the CovGE Bootstrap has a higher false positive rate that you reported before. Where are these false positives coming from? Is it from the 2 environment scenarios? I think we need to dig into this some more, because if the false positive rate is >5% that is higher than we would expect based on random!

Also what are the counts these percentages are based on? I just want to confirm we have good coverage of the cov = 0 and gxe = 0 cases.

| | CovGE Permutation | CovGE Bootstrap | GxE Permutation | GxE Bootstrap |
| ---| ---| ---| ---| ---|
| True Positive | 8.21% | 62.26% | 43.3 % | 94.79% |
| True Negative | 0.71% | 0.62% | 0.74% | 0.20% |
| False Positive | 0% | 0.09% | 0% | 0.54% |
| False Negative | 91.1% | 37.03% | 55.96% | 4.47%|
|---|---|---|---|---|
| False Negative Rate | 91.73% | 37.3% | 56.38% | 4.5% |
| False Positive Rate | 0 % | 12.68% | 0% | 72% |

## Parameter Coverage
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Simulation_10.10.2020/HexPlotPanel.png)

## CovGE vs. GxE
The setup here is confusing, because the actual number of total samples is very different for the top and bottom graphs. Might be easier to chat about.
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Simulation_10.10.2020/Cov_GxE.png)

## Heatmap for Std. Dev = 0.5
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Simulation_10.10.2020/HeatMap_LowStdDev.png)

## Omega^2 vs. GxE with Estimated Marginal Means
Remove the 1:1 line. I don't think there is a reason to expect these should be related 1:1
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Simulation_10.10.2020/OmegaVsGxE.png)

## Confusion Panels - Raw Data
Add a legend. Also where are the false positives for bootstrap? They're in the confusion matrix. I think there's a bug in your code where the FP and FN are based on the same thing, and not the criteria in panel of the plot.
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Simulation_10.10.2020/ConfusionPanels_16Plots.png)

## Confusion Panels - Mean Data
Add a legend
Something weird is going on here. Why aren't the CI centered on the estimates for the CovGE?
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Simulation_10.10.2020/ConfusionPanels_16Plots_means.png)

## Covariance: Means Vs. Raw Checks
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Simulation_10.10.2020/Cov_MeansVsRaw_panel.png)

## GxE: Means Vs. Raw Checks
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Simulation_10.10.2020/GxE_MeansVsRaw_panel.png)







