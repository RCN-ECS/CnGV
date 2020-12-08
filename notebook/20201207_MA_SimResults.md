#Sanity Checks of Results

Over the course of last week I ran 10 replicates of the new parameters with a more targeted approach to improve sampling of CovGE = 0; GxE = 0; and also to gain more representation of the moderate estimates for Common Garden design. 

For reciprocal transplant, I now have 1241 outcomes in which true covariance is 0, and 1242 outcomes in which true GxE = 0. 
For common garden I have 2018 outcomes in which true covariance is 0, and 1834 outcomes that true GxE = 0.

In terms of sampling moderate values of each: 
Full Reciprocal transplant has 4844 results with CovGE between 0.2 and 0.6 or -0.2 and -0.6, and GxE has 4476 with GxE_emm between 0.3 and 0.6. 
Common Garden has 7305 results with CovGE between 0.2 and 0.6 or -0.2 and -0.6, and 13215 results with GxE_emm between 0.3 and 0.6.

###Parameter Coverage 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.1.20/Figures/12.7.ParameterCoverage_hex.png)

###Power Heatmap
With low total samples excluded

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.1.20/Figures/12.7.Heatmap.png)

####Full Reciprocal Transplant -- Overall
|---| CovGE - Permutation | CovGE - Bootstrap | GxE - Anova | GxE - Permutation | GxE - Bootstrap |
|---| --- | --- | --- | --- | --- |
False Negative | 4632  82.55 | 2517 - 44.86 |2066 - 36.82|2183 -38.91 |202 - 3.60 |
False Positive   | 1  0.02 | 60    1.07  |29 - 0.52|3 - 0.05 | 558-0.94 |
True Negative | 616  10.98  | 557  9.93  | 587 -10.46| 613 -10.92 | 58 -1.03  |
True Positive | 362   6.45 |2477 44.15 | 2929 -52.20 |2812-50.12 |4793- 85.42 |
False Negative Rate | 0.93 | 0.5|  0.41|0.44| 0.04|
False Positive Rate |0.00 |0.1|0.05|0.00| 0.91|


####Common Garden Design
|---| CovGE - Permutation | CovGE - Bootstrap | GxE - Anova | GxE - Permutation | GxE - Bootstrap |
|---| --- | --- | --- | --- | --- |
False Negative | | | | | |
False Positive  | | | | | |
True Negative |   | | | | |
True Positive |  | | | | |
False Negative Rate | | | | | |
False Positive Rate | | | | | |

#### False Positive Rate 

I filtered the data to the reflect the same window as the power analysis (CovGE between 0.2 and 0.6 or -0.2 and -0.6; and GxE between 0.3 and 0.6). When divided this way, 

Labeled with the number of observations within each category

**Full Reciprocal Transplant**
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.1.20/Figures/12.7.FRT_FalsePos_Trimmed.png)

**Paired Common Garden**
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.1.20/Figures/12.7.CG_FalsePos_trimmed.png)
