# 20201109 Meeting notes

* Make the manuscript as complete as possible through Results section, as if we were submitting for peer review

* Check that all approaches for raw data vs. means data are as parallel as possible 
  * Standardization by group means vs. all phenotypic means - what we had talked about back in August vs. in manuscript - Molly confirmed that everything is standardized by group means
  * Hypothesis test as delta-GxE=0 (raw data?) vs. delta-GxE = G+E (means data?)
    * At one point I remember digging into the weird null distribution for GxE being created by the permutation: https://github.com/RCN-ECS/CnGV/blob/master/notebook/20200616_KEL_ExploreMeansPermutation.R 
    * Review this post together with Molly
  
* Confusion matrix presentation - things to discuss and decide
  * While the "overall" confusion matrices are useful for us, I am hesitant to present them in the main text because the true positive / 
  false negative rates depend on the study design and magnitude of the CovGE
    * For example, Fig 4 could be limited to all studies with 128 total sample size or something like that
  * Determine what designs false positives and false negatives are coming from
  * Compare ANOVA confusion matrix to our GxE approach
  * Confusion matrix rates are not consistent across Recip. Trans vs. Paired common garden

* Make a full list of all other things that we need to check and discuss - or not passing sanity checks
  * Write out our equations for the population parameter calculation vs. the sample estimates - especially for Cov-GE.
  * Explain the clustering approach for paired common garden sims more clearly
  * Make an outline for the discussion with Geoff
  * Make sure to have 1000 sims of CovGE = 0 for FPRs?
  * Make sure to have 1000 sims of GxE = 0 for FPRs?
  * Paired common garden why actual cov GE not covering large values?
  

* List of questions for https://github.com/RCN-ECS/CnGV/blob/master/notebook/20201016_MA_FullSimResults.md
  * https://github.com/RCN-ECS/CnGV/blob/master/notebook/20201016_MA_FullSimResults.md#confusion-matrix Error rates for full reciprocal transplant and common garden not matching
  * https://github.com/RCN-ECS/CnGV/blob/master/notebook/20201016_MA_FullSimResults.md#parameter-coverage Paired common garden does not appear to be covering large values of CovGE. What are the actual sample sizes for each cell in this figure used to calculate the power from: https://github.com/RCN-ECS/CnGV/blob/master/notebook/20201016_MA_FullSimResults.md#figure-5
    * Also in this figure, within a panel, the power for a design should increase with sample size. However, it's not always consistently the case, which means we have a mistake or some kind of sampling error due to a small number of simulatiosn in that cell.
  * I don't understand this graph https://github.com/RCN-ECS/CnGV/blob/master/results/Sim_10.30.2020/11.1.Cov_RawVsMeans_Env2.png
  * https://github.com/RCN-ECS/CnGV/blob/master/notebook/20201016_MA_FullSimResults.md#gxe-permutation-p-values-vs-anova-p-values I'm confused about the upsidedown "L"-shaped pattern here, and why the patterns in the P-values are not consistent with the raw data vs. group means approaches.
 
  
* Make a new notebook post so that we have the results to compare to last post


# 20201110 Meeting notes

* we talked about manuscript writing is good idea
* manuscript could be good without GxE, Molly definitely wants it in there

GxE issue
* permutation approach is different for means data vs. raw data
* the GxE hypothesis test that we have for the means data is specifically for the case when means and SE are provided
* it doesn't appear that we have a function to create the null distribution for raw data
* start with one dataset - look at the null distributions for GxE hypothesis test that you're using
 * choose a dataset that has a significant ANOVA P-value, but the the P-values are reall different for means vs. raw data - make sure it has a decent total sample size - what does this show? discuss check with Katie
 * use that dataset to develop and test code
 * then run a test set to see if solves the problem
 * if necessary, rerun sims
 
Confusion matrices
* Come up with a way to present it for specific designs that we recommend
* Do we have enough sims to estimate false positives? Typical sim studies recommend 1000
* Molly thinks this should be a representation of random sampling - we discussed how to structure it so that 
* Start with Figure 5, remake with false positive and false negative rates. Specify number of sims in each cell. Then we can talk about how we want to "hone in" if we low samples for the false positive rates.
* Are the 2x2 driving the high FPR?
* Discussed some of the differences in the confusion matrices between the reciprocal transplant vs. common garden data, is it different in the way the parameter space is sampled or something else about the design?
* Discussed benefit of comparing ANOVA confusion matrix to our GxE approach, in terms of convincing reviewers
* Let's figure out how we're going to present it, then figure out what sims we would need to make that presentation, then just run those sims.


population parameter calculation vs. the sample estimates







