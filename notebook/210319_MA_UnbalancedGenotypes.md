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

## KEL comments

I worry about just adding a discussion of imbalanced designs, because this type of imbalance actually introduces a systematic bias depending on whether genotypes are missing from the lower or higher category. It’s not just introducing random error. There may be a third option, which is to expand the calculation for these special cases, but I can’t know for sure how to guide you without more information on exactly how you made those plots. When you removed 3 genotypes was it always from the same kind of genotype (e.g. always saltwater)? Or was it random? I want to make sure we are simulating the former, so we can investigate bias. I was also wondering, what does comparing the sampling estimate relative to itself tell us, and is it confounded with the random error we expect in the sampling estimate? 

Here’s what I’m thinking: could you re-run your results with a couple of small tweaks to what you already ran? Run the results once for imbalance in the lower of the two categories of genotypes (e.g., only saltwater), then again for imbalance in the higher of the two categories genotypes ((e.g., only freshwater), which will help us understand systematic bias. (Does that make sense?) For each of those two types of imbalance, could you re-plot those results all against the truth-known population parameter in the following three ways:

(i) difference between population param balanced (e.g. 8x2 is 4x2 of one genotype category and 4x2 of another genotype category) and population param UNbalanced (e.g. e.g. 5x2 is 1x2 of one genotype category and 4x2 of another genotype category) ,  (this should show the bias, as long as all the replicates remove the same category of genotypes)

(ii) difference between population param balanced vs. sampling estimate balanced (e.g. both e.g. 8x2 is 4x2 of one genotype category and 4x2 of another genotype category)  (this should show the random sampling error in the absence of imbalance)

(iii) difference between population param balanced (e.g. 8x2 is 4x2 of one genotype category and 4x2 of another genotype category)  vs. sampling estimate UNbalanced, (e.g. e.g. 5x2 is 1x2 of one genotype category and 4x2 of another genotype category) , (this should be a combination of the systematic bias introduced by imbalance and random sampling error )

Let me know what you think and if any of this makes sense to you. Messages can be hard sometimes
