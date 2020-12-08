# Sanity Checks of Results 

I ran 10 replicates of the new parameters with a more targeted approach to improve sampling of CovGE = 0; GxE = 0; and also to gain more representation of the moderate estimates for Common Garden design. 

For reciprocal transplant, I now have 1241 outcomes in which true covariance is 0, and 1242 outcomes in which true GxE = 0. 
For common garden I have 2018 outcomes in which true covariance is 0, and 1834 outcomes that true GxE = 0.

In terms of sampling moderate values of each: 
Full Reciprocal transplant has 4844 results with CovGE between 0.2 and 0.6 or -0.2 and -0.6, and GxE has 4476 with GxE_emm between 0.3 and 0.6. 
Common Garden has 7305 results with CovGE between 0.2 and 0.6 or -0.2 and -0.6, and 13215 results with GxE_emm between 0.3 and 0.6.

### Parameter Coverage 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.1.20/Figures/12.7.ParameterCoverage_hex.png)

### Power Heatmap
Low total samples excluded

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.1.20/Figures/12.7.Heatmap.png)


#### False Positive Rate 

I filtered the data to the reflect the same window as the power analysis (CovGE between 0.2 and 0.6 or -0.2 and -0.6; and GxE between 0.3 and 0.6). When divided this way, I miss nearly 100% of the results in which trueCovGE and trueGxE = 0. This wildly inflates the false positive rate because there are few, if any, true negatives. It doesn't make sense to filter these data to the same window as the heatmap.

**Full Reciprocal Transplant**
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.1.20/Figures/12.7.FRT_FalsePos_untrimmed.png)

**Paired Common Garden**
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.1.20/Figures/12.7.CG_FalsePositive_untrimmed.png)

#### Full Reciprocal Transplant -- Table Format
|---| CovGE - Permutation | CovGE - Bootstrap | GxE - Anova | GxE - Permutation | GxE - Bootstrap |
|---| --- | --- | --- | --- | --- |
False Negative | 4632  82.55 | 2517 - 44.86 |2066 - 36.82|2183 -38.91 |202 - 3.60 |
False Positive   | 1  0.02 | 60    1.07  |29 - 0.52|3 - 0.05 | 558-0.94 |
True Negative | 616  10.98  | 557  9.93  | 587 -10.46| 613 -10.92 | 58 -1.03  |
True Positive | 362   6.45 |2477 44.15 | 2929 -52.20 |2812-50.12 |4793- 85.42 |
False Negative Rate | 0.93 | 0.5|  0.41|0.44| 0.04|
False Positive Rate |0.00 |0.1|0.05|0.00| 0.91|


#### Common Garden Design -- Table Format
|---| CovGE - Permutation | CovGE - Bootstrap | GxE - Anova | GxE - Permutation | GxE - Bootstrap |
|---| --- | --- | --- | --- | --- |
False Negative |8082 -77.61 | 5301 - 50.90 | 1315 - 29.33 | 2670 -   25.64 | 0 - 0.00 |
False Positive  |  12 - 0.12 | 72 - 0.69 | 18 - 0.40|   0 - 0.00 | 814 - 7.82 |
True Negative | 891 - 8.56  | 831 - 7.98|  469 - 10.46|  814 - 7.82 | 0 - 0.00 |
True Positive |1429 - 13.72   | 4210 - 40.43 | 2682- 59.81|  6930 - 66.55 | 9600 - 92.18|
False Negative Rate |  0.85| 0.56| 0.33|0.28 | 0|
False Positive Rate | 0.01| 0.08|  0.04| 0.00| 1|




