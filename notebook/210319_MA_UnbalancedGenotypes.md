# Exploring Unbalanced Common Garden Designs

It was identified that several studies in the meta-analysis have unbalanced designs, particularly in common garden designs, where there can be different numbers of genotypes across the 2 environments. 

To see whether estimated marginal means are robust to unbalanced designs, I ran a few quick tests across different magnitudes of CovGE and different total numbers of genotypes to see how removing genotypes affects CovGE estimates. I tested both population and sample estimates

Answer: Yes, it affects the Population AND sample estimate. The more missing genotypes relative to total, the more variation in the sample estimatE, which makes sense. 

Does removing a random genotype to balance out designs fix the problem? No. I played around with omitting genotypes from the other native environment to see if restoring "pairs" would bring the estimate back the original. It does not. (Also makes sense - removing data will affect CovGE!) 

In the plots below, I am showing how removing 1 (top), 2 (middle), or 3 (bottom) genotypes affects the population or sample estimate. The variation gets larger as the number of genotypes removed gets larger relative to the total number.

### Population Estimate

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_3.10.21/3.21.Pop_MissingGenotypes.png)


### Sample Estimate
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_3.10.21/3.19.MissingGenotypes.png)
