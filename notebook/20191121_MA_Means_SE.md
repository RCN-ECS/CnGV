#Test of Meta-Analysis Means and Standard Error

Several studies present summary statistics rather than raw data (9/20 extracted studies as of 11/15) that would be good to include in the analyses and manuscript. But to estimate covariance with confidence intervals, we need a greater resolution of the sampling distribution. 

The sampled distribution of means is normally distributed with the standard error equivalent to the standard deviation that mean to the sampled distribution, so if we assume the central limit theorem, we may be able to simulate all combinations of probable slopes, estimate means and confidence intervals around those slopes, and use that data to estimate covariance.

I used simulated data to test whether this would reliably mimic original data. I took the following approach:

1. Simulated data - both categorical and continuous
```#start

```
2. Estimate means and standard error of simulated data.
3. Bootstrap slopes and confidence intervals
4. Calculate covariance of bootstrapped slopes 
5. Compare covariances. If it works, the covariances should match. 

