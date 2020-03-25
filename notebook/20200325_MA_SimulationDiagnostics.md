# Simulation Diagnostic Explorations: 

The results of my simulations are giving me some strange values. GxE estimates near zero are significant (and they shouldn't be), while other GxE values are high, but not significant (and they should be). Additionally, there is some weirdness with the covariance estimates. 

I chose the following 3 scenarios to look further at the null distributions as a way of figuring out these weird p-value situations. 

Index | delta_env | delta_gen | sample_size | n_pop | std_dev | interaction | reason 
---|---|---|---|---|---|---|---
1 | 0.0| -1|5|5|0.1|0| The GxE is near zero yet p <0.01. Also, the true_cov is -.21 with p = 0.27 but is weird because index = 2 is higher true_cov but lower pvalue
2 | 0.5| 0|5|5|0.1|0| The cov estimate is higher than index=1 but the pvalue is lower
3 | 1| -1|5|5|0.1|1.5|This one has a high GxE and is significant - this one is behaving the way I would expect for GxE. 
