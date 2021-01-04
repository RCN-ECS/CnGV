

# Results from most recent simulation runs 

I ran a (hopefully) final round of simulations last week. The plots here are arranged as they would appear in the paper. As a reminder, figure 1 and 2 are heuristic plots and not generated from simulated data. 
Supplementary materials are located below the main 7 plots.

## Figure 3

* KEL note: looks good. Add CI for permutation and ANOVA p-values for GxE. I worry the font will be too small for a 6-inch page. The significant GxE and non-sig GxE appear to be reversed?

**Molly: pick a few to show and leave it - maybe 2x2 and 4x4**

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/PhenotypePlotPanel.png)

## Figure 4

* KEL note: see notes in other post.

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/12.21.FPR_panel.png)

## Figure 5

* KEL note: I need more information to interpret this graph - effecct size? power? permutation or bootstrap?

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/12.28.HeatmapPanel.png)

## Figure 6

* KEL note: this is beautiful. Font size is great. Suggest adding CI to cov G,E. I thought we had a sig GxE for your data before? Did you have that in your study? 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/12.21.RealData.panel.png)

## Figure 7

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/12.21.Cov_GxE_tradeoff_panel.png)

## Supplemental Figures 

### Supplemental Figure 1

* KEL note: need context on data filtering to interpret.
 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/12.21.Hex.ParameterCoverage.png)

### Supplemental Figure 2


* KEL note: how did you evaluate power for CovGE near 1 for common garden with 16 genotypes? There doesn't seem to be any data simulated there. Likewise, how did you evaluate power for GxE near 1? Doesn't seem to be any sims there.
 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/12.21.FalseNeg.FullPanel_effectSizes.png)

### Supplemental Figure 3

* KEL note: I need more info to interpret the legend of this graph. Would probably be easier to just talk through this one. Also why is there only one line?

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/12.21.OmegaVsGxE_panel.png)

### Supplemental Figure 4

* KEL note: Is it possible the lenght of the CI is the same, but there is some kind of bias? (e.g. asymmetric CI b/c bounded by 0 and 1?). Should we plot lower and upper CI instead?
* are you going to do the same thing for GxE?

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/12.21.MeansVsRaw_panel.png)

### KEL Other Thoughts

I'm trying to think of other things we've plotted in the past that have helped us (i) find bugs: 
- e.g. P-value comparisons between ANOVA and GxE
- population vs. sample estimates

and (ii) go with other supplemental parts of the paper:
- the variance partitioning plots and examples (https://github.com/RCN-ECS/CnGV/blob/master/notebook/20201019_VarPart_KEL.Rmd)
- the PL statistic comparison
- 
