## Looking back

We've done a lot, and it can be hard to make decisions about how to focus the manusccript given how much we have. 
Therefore I propose that we create a list of questions that we think the manuscript should answer (with sanity checks, 
which we may or may not include in the manuscript). Then we should go through the plots and test if they answer these questions.

Here is a start to that list (up for discussion) (not necessarily in order they will be in the manuscript):

### 1. Why do we need a measure of CovGE? (I think this is well addressed in the manuscript: we don't have one)

### 2. Why do we need a meausure of delta GxE, when we already have an ANOVA framework which is widely used?
  - Answer: equivalent effect sizes of GxE may have different w^2 from an ANOVA depending on the residual variance. Therefore we need a measure that is comparable among datasets with different residual error and sampling designs.
     - which of our plots illustrates this?
 - Also: CovGE can exist in the presence of GxE, but based on first principles is also inversely related, so we need a measure of both that is 
 comparable across different datasets.

### IMPT: Are the sampling estimates accurate and unbiased estimations of the population parameters?
- which of our plots show this?

### 3. How should people test for the significance of CovGE and GxE estimates? How should people describe the uncertainy of their estimates?
  
CovGE
  - Answers:
    - permutation can be used to test for significance, and bootstrap for CI
      - which of our plots illustrates this?
    - all plots which report P-values should also report CI in the estimate in the plot, a table, the main text, or the figure legend
      - for designs that are underpowered for permutation, we commonly find the CI do not include zero. Studies that find this pattern are likely to have a true effect size, and should pursue larger designs for significance.
  
  GxE
  - Answers: permutation is OK, but the bootstrap estimate is unreliable. So we don't have a good method for GxE uncertainty at this time.
    - which of our plots illustrates this?
  - how does the FPR our metric compare to ANOVA?
     - which of our plots illustrates this?
  - how does the power our metric compare to ANOVA? 
    - which of our plots illustrates this?

### 4. What are the costs and benefits of variance partitioning vs. effect size?
  - datasets with equivalent CovGE may be driven by Vg explaining more variance than Ve, or by Ve explaining more variance than Vg. 
  Understanding whether CovGE relates to questions regarding the evolution of fixed genetic effects (Vg) vs. plastic effects (Ve) and may also affect predictions for how populations with equivalent levels of CovGE, but different contributions of Vg vs Ve, respond to changing environments
  - When CovGE is maximized, the contributions of Vg and Ve to phenotypic variance are equal
      - Which of our figures illustrates the above concepts?
      - once the rest of the manusccript and supp is together, KEL will help with supp figures for this
  - Our measure of CovGE is robust to differences in residual variation among datasets
      - Which of our figures illustrates this?

### 5. What kinds of designs do we recommend for different scenarios?
  - Is FP related to total sample size? For the same total sample size, are certain designs worse than others?
    - which of our plots illustrates these?
  - How does power increase with total sample size? For the same total sample size, are certain designs worse than others?
      - which of our plots illustrates these?
  - If people are working in a system with a really strong pattern of CovGE, how small a study can they get away with?
    - which of our plots illustrates how power increases with effect size for boostrap and permutation?
  - How do the sample designs/sizes of the Albecker and Trussel studies relate to our power analysis? I think we've got this nailed for Trussel, I'm just curious about whether Molly's design is predicted to have power under our analysis.
      - **MOLLY NOTE: I have 7 genotypes and 2 experimental environments with an average of 27 samples per treatment (variable due to survival). Thus my total sample size is 378, plenty of power.**

### 6. For the purposes of metanalysis, how do the means data and raw data compare?
  - same questions go for FP and power for the means data. But the simple thing may be to calculate these for each scenario in the raw data and just compare
  in a scatterplot to the means data. For things that fall off the 1:1 line we may need to dig in.
  - in order to get this paper into review, I suggest that we discuss moving this to the metanalysis paper. If we come across issues, we can resolve them in the revisions of the raw data paper.
  - does it matter that the CI for the GxE means matches the raw data, if they are not accurate in the first place?
  
### 7. How does CovGE compare to other stats that have been used for CovGE (PL)?
   - in order to get this paper into review, I suggest that we discuss moving this to the metanalysis paper. Especially since the PL was a metanalysis and our results will almost certainly overturn their conclusions.

### Additional things to consider:

- include a legend for a color scheme
- use different color schemes for different things
  - Example: In your FP bar plots, the color scheme is used for sampling design. But in the power heatmaps below it, the same color scheme is used to indicate power between 0 and 1. 
- use the same color scheme for the same thing across multiple plots
- before you share plots with me again, please add a figure legend. Eg for heatmaps, indicate the effect size that is being visualized and the filtering criteria (is it all sims or only the ones with greater variation?). Alternatively, just writing the manuscript might save time.

- Let's come up with a concise code for the sampling designs that includes the total sample size. The reason is that total sample size, rather 
than nuances in sampling design, seem to drive everything (although I may be wrong about this). Example:
"RT
2Gx2Ex8rep
32"
"CG
8Gx2Ex4rep
64"
For barplots, should we order by total sample size?

- Keep in mind that KEL started a LaTeX supp mat in the google drive folder, based on conversations about the supp mat format that we had last summer.


