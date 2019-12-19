#Expand Meta-analysis Data Collection

A present frustration of the meta-analysis is that its dying the death of a thousand qualifications. 

Right now we can only deal with data that are categorical (excludes all continuous/common garden designs), have matching numbers of genotypes and environments, and unless I come up with something to fix the means/SE code, we can only extract studies that provide raw data. 

To reduce some of these restrictions, I came up with a way to use data where the number of genotypes does not match the number of environments. 

Several studies have multiple genotypes across a different environmental gradient, but only expose the genotypes to 2 environments. These can't be simply analyzed because covariance requires an Emean for every Gmean.

For example, 

| Native_env_cat|	Exp_env_cat|
|---|---|
|High_warm	|Warm|
|High_warm	|Cool|
|High_cool |Warm|
|High_cool	|Cool|
|Low_warm	|Warm|
|Low_warm	|Cool|
|Low_cool	|Warm|
|Low_cool	|Cool|

