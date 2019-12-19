# Expand Meta-analysis Data Collection

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

In this case, there are multiple genotypes from low and high locations that are each exposed to warm and cool conditions. 
This would be excluded from the current approach because the number of genotypes doesn't match the number of experimental environments.

Because these 4 genotypes are plotted against the same 2 categorical environments, I decided that if I divide the dataset into the 2 groups of genotypes (in this case, genotypes from High and Low), then I can then separately estimate gmeans and emeans for each genotype, then put the datasets back together to create one large Gmean/Emean dataset to estimate covariance.

This requires that when data is extracted from papers with environments like this, we need to divide up the Native_env_cat. 

I have adjusted the code to work if we add a column called `Native_env_cat2` *JUST TO THE DATASETS THAT NEED IT* 

Enter MAIN environment (one that matches to the `exp_env_cat`) in `Native_env_cat`. 

Then put the second qualifier in the `Native_env_cat2` column. If there is a third genotype pair, add another column called “`Native_env_cat3`”. 

So the example above should look like: 
|Native_env_cat	|Native_env_cat	|Nat_env_factor	|Exp_env_cat|
|---|---|---|---|
|Warm|	High|	N_1|	Warm|
|Warm|	High|	N_1|	Cool|
|Cool|  High| N_2|  Warm|
|Cool|	High|	N_2|	Cool|
|Warm|	Low|	N_3|	Warm|
|Warm|	Low|	N_3|	Cool|
|Cool|	Low|	N_4|	Warm|
|Cool|	Low|	N_4|	Cool|

