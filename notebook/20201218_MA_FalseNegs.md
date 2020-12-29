
I have gone back and forth on whether to use bootstrapping or permutation to assign significance to CovGE. So here are some plots to guide decisions. 

## KEL notes
I like all of these - great job!

- Overall thoughts
  - How were the effect sizes binned? was it >0-0.25, >.25-0.5, etc? Some of the GxE had > 1, just making sure we didn't miss them
  - the same general hard thing of how to represent the total sample sizes - e.g. the yellow bars under 2 populations on the reciprocal transplant is not the same total sample size as the yellow bars under 2 populations on the bottom.
  - It might help to have the GxE permutation next to ANOVa so we can show their performance. ANOVA generally has higher power and lower false positive rates. We should discuss how to deal with this for the metanalysis.
  - For the paper I'm thinking (happy to chat more):
    - for CovGE we should use bootstrap and permutation to display and discuss results 
    - for GxE we should use permutation P-value and ANOVA P-value to display and discuss results

- False positive graphs
  - I like how these are displayed because it helps to see the different approaches. but it makes it hard to compare to the true positive graphs.
  - should we remove the common garden 8 pops x 2 sample size? I remember we discussed a lower boundary for total sample size, but don't remember what it was.

- True positive graphs
  - I'm wondering if there is a better way to represent these with line plots instead of the bars. Typically power (y axis) is a function of effect size (x-axis) with different designs or methods represented by different lines. If we go by design, that could work well since we don't have balance with population number and sample number.
  - A time-saving approach I often use is to draw out some different ways of representing the data with colored pens and discussing which one works the best with different people, which can save time with writing code for a million different graphs.
  
 Other thoughts
   - It would be helpful to see the false positives and power plots together, one on top of the other.  Could the design for the false positives could be made to be parallel to the design for the true positives? It took me a while to get my brain around both types of graphs. Up for discussion or maybe lab opinion. For example, if someone wanted to plan a design, they should be able to see what power and FP rates they would have. I'm not sure if this is feasible as it will depend alot on how we visualize power. Just something to keep in mind.


## False Positives 

Per Katie's instruction - False positives are subject to residual error and sampling design. 

Below I have plotted false positive rates for each metric (Covariance + bootstrap, Covariance + permutation, etc.) using the larger residual error (standard deviation = 1) according to sampling design. These are as expected. 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/12.18.FRT_FalsePos.png)

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/12.15.CG_FalsePositive.png)

# False Negatives -- Effect Size
Per Katie's instruction - False positives are subject to residual error, sampling design, and effect size. With that in mind, I've plotted Power (1-beta) as bar plots to visualize how power changes according to effect size (divided into 4 effect sizes), design, number of populations, and sample size. 

**Full Reciprocal Transplant**
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/12.18.FRT.GxE.Power.png)

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/12.18.FRT.Cov.Power.png)

**Common Garden Design**
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/12.18.CG.Power.GxE.png)

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/12.18.CG.Cov_Power.png)

# False Negatives -- Total Sample Size
As suggested by Katie, I've also plotted Power (1-beta) according total sample size. Same plots as above, just showing 128 total samples

**Full Reciprocal Transplant**
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/12.18.FRT_GxE_128samples.png)

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/12.18.FRT_Cov_128Samples.png)

**Common Garden Design**
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/12.18.CG.GxE.128Samples.png)

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_12.15.20/12.18.CG.Cov.128Samples.png)


**Conclusions**

I think using permutation for CovGE and GxE are both good for the paper


