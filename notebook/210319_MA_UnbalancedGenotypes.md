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



### MA options:

ultimately the options are to: Remove all Common Garden studies that are imbalanced, or include some discussion of imbalanced designs and discuss that imbalance will affect CovGE. I lean toward the latter option, but given the uncertainty around it, I’d understand if we needed to simply remove all imbalanced studies.

## KEL comments

There may be a third option, which is to expand the calculation for these special cases, but I can’t know for sure how to guide you without more information on exactly how you made those plots. When you removed 3 genotypes was it always from the same kind of genotype (e.g. always saltwater)? **MA: YES, THE SAME GENOTYPE(s) EACH TIME** Or was it random **MA: NOT RANDOM**? I want to make sure we are simulating the former, so we can investigate bias. I was also wondering, what does comparing the sampling estimate relative to itself tell us, and is it confounded with the random error we expect in the sampling estimate? 

I worry about just adding a discussion of imbalanced designs, because this type of imbalance actually introduces a systematic bias depending on whether genotypes are missing from the lower or higher category. It’s not just introducing random error. If the bias is systematic (e.g. brings CovGE toward 0), then these imbalanced studies would be underpowered. If the bias is systematic in the other direction (e.g. larger CovGE), then these imbalanced studies would be more likely to lead to false positives. **MA: MY INITIAL READ WAS THAT REMOVALS BROUGHT ESTIMATES CLOSER TO ZERO BUT LETS SEE IF THIS BEARS OUT**

Here’s what I’m thinking we could do to better understand the bias: could you re-run your results with a couple of small tweaks to what you already ran? Run the results once for imbalance in the lower of the two categories of genotypes (e.g., only saltwater), then again for imbalance in the higher of the two categories genotypes ((e.g., only freshwater), which will help us understand systematic bias. (Does that make sense?) For each of those two types of imbalance, could you re-plot those results all against the truth-known population parameter in the following three ways:

(i) comparison of population param balanced (e.g. 8x2 is 4x2 of one genotype category and 4x2 of another genotype category) and population param UNbalanced (e.g. e.g. 5x2 is 1x2 of one genotype category and 4x2 of another genotype category) ,  (this should show the bias, as long as all the replicates remove the same category of genotypes)

## Population Balanced vs. Unbalanced: 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_3.10.21/3.23.Population_GenotypeRemoval.png)
 
 **CAPTION** In the above plot, I have plotted the population CovGE (black) with the Population CovGEs after genotypes have been removed. I remove a genotype from both "sides" (i.e., freshwater vs. saltwater groups). I call them either "A-side" or "B-side". The effects of genotype removal on CovGE differ according to where they are removed. The more genotypes that are removed, the greater the effect on CovGE.  
 
(ii) comparison of population param balanced vs. sampling estimate balanced (e.g. both e.g. 8x2 is 4x2 of one genotype category and 4x2 of another genotype category)  (this should show the random sampling error in the absence of imbalance)


(iii) comparison of population param balanced (e.g. 8x2 is 4x2 of one genotype category and 4x2 of another genotype category)  vs. sampling estimate UNbalanced, (e.g. e.g. 5x2 is 1x2 of one genotype category and 4x2 of another genotype category) , (this should be a combination of the systematic bias introduced by imbalance and random sampling error )

## Population vs. Sample - Balanced and Unbalanced: 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_3.10.21/3.23.PopV.Sample_Unbalanced.png)
 **CAPTION** The first panel shows population CovGE (black) with 100 sample estimates (grey points). The consecutive plots show how sample estimates are affected after removing 1, 2, or 3 genotypes from either group (A-side group is Dark Blue, B-side group is Lighter blue). The variation in CovGE sample estimates seems to increase, but I'm not super sure 1-2 genotype removals has that big of an effect on variation in CovGE sample estimates.
 
* Note that if we need to be able to see the direction of the bias for the two things being compared in the plots.* Let's take the case that there is systmatic bias, in which the CovGE will be closer to 0 in the imbalaneced cases. If you simulate postivie CovGE, the difference between the balanced and unbalanced cases will be a decrease. If you simulate negative CovGE, the difference between the balanced and unbalanced cases will be an incnrease. The bias is systematic (CovGE goes towards 0), but if you only plot the histgram of the difference, the bias would not be evident.

Let me know what you think and if any of this makes sense to you. 
