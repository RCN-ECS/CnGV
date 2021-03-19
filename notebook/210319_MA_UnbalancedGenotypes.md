# Exploring Unbalanced Common Garden Designs

It was identified that several studies in the meta-analysis have unbalanced designs, particularly in common garden designs, where there can be different numbers of genotypes across the 2 environments. 

To see whether estimated marginal means are robust to unbalanced designs, I ran a few quick tests across different magnitudes of CovGE and different total numbers of genotypes to see how removing genotypes affects CovGE estimates. 

Answer: It affects the estimate. The more missing genotypes relative to total, the more variation, as expected. I played around with omitting genotypes from the other native environment to see if equalizing the number of genotypes returned the estimate to the  full desiggn
