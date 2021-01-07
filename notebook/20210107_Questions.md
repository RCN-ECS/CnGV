## Looking back

We've done a lot, and it can be hard to make decisions about how to focus the manusccript given how much we have. 
Therefore I propose that we create a list of questions that we think the manuscript should answer (with sanity checks, 
which we may or may not include in the manuscript). Then we should go through the plots and test if they answer these questions.

Here is a start to that list (up for discussion) (not necessarily in order):

1. Why do we need a measure of CovGE? (I think this is well addressed in the manuscript, we don't have one)

2. Why do we need a meausure of delta GxE, when we already have an ANOVA framework which is widely used?
  - Answer: equivalent effect sizes of GxE may have different w^2 from an ANOVA depending on the residual variance. Therefore we need a measure that is comparable among datasets with different residual error and sampling designs.
     - which of our plots illustrates this?


3. How should people test for the significance of CovGE and GxE estimates? How should people describe the uncertainy of their estimates?
  CovGE
  - Answers:
    - permutation can be used to test for significance, and bootstrap for CI
      - all plots which report P-values should also report CI in the estimate in the plot, a table, the main text, or the figure legend
      - for designs that are underpowered for permutation, we commonly find the CI do not include zero. Studies that find this pattern are likely to have a true effect size, and should pursue larger designs for significance.
  
  GxE
  - Answers: permutation is Ok
  - how does the FPR our metric compare to ANOVA?
     - which of our plots illustrates this?
  - how does the power our metric compare to ANOVA? 
    - which of our plots illustrates this?

4. What are the costs and benefits of variance partitioning vs. 

4. What kinds of designs would we recommend for different scenarios?
