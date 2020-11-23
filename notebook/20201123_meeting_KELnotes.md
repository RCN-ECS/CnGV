Katie's notes from 10am meeting on 11/23/2020 with Molly

- Molly updated code to make sure we had population and sample covariance and GxE

- Reran a subset of sims

- For confusion matrix tables - the sanity check is that FP + TN should add up to the same number for COVGE (or GxE/anova), 
and FN + TP should add up to the same number

- We double checked sample CovGE and population CovGE, since we updated that code
  - looks like our sample estimate tends to overestimate CovGE, dig into why some of them are way off?

- Katie will send Molly gist on R code to put at top of scripts to automatically update packages
  - https://gist.github.com/DrK-Lo/a945a29d6606b899022d0f03109b9483 
  - I put this at the top of all my scripts

- How confident are you that you will be able to get the full parameter space for the common garden?
  - Clustering term - Katie doesn't understand how this works
  - Molly explained e_pop is a random normal that affects the intercepts of the genotypes
  - Katie suggests focusing on creating scenarios with large true_cov between -1 and 1, so we can compare to reciprocal transplant
  - We need 1000 simulations within a given range (e.g 0.2- 0.6) between reciprocal transplant and common garden
  - Maybe set three levels of sd for e_pop? sd = 0 (should give -1,1), a medium number that gives covGE={-0.7,0.7}

- For publication, we need to "dig in" to the power/FN and FP rates. The "overall" tables are useful for us, but not for publication
  - Lotterhos and Whitlock 2014 - https://onlinelibrary.wiley.com/doi/epdf/10.1111/mec.12725 
      - For example, Fig 2. This paper has been cited alot, so I must have done something right in describing the results?
      - Happy to discuss this paper before tackling the stuff below, if it helps
  - False positive rates will be a function of sampling design and random noise
  - True positive rates (or false negative rates) will be a function of sampling design, random noise, and effect size
  - Power = 1- FNR (so don't need to present both)
  - We need to make sure we have enough sims within a cell (e.g., high resolution) so that differences in the rates in different situations are not caused by random error
      - for example, if there are only 20 sims within a cell, then our resolution is only 0.05
      - 1000 is a good aim
  - How does power to detect a specific effect increase with total sample size?
  - How does power increase with effect size, as a function of a single total sample size?
  - Is permutation really that bad in a well designed experiment with strong covGE in the system?
  - Happy to help sketch out some plots and talk about them before getting into the weeds
   
