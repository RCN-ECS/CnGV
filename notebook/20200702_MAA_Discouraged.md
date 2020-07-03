
Well, I've encountered another problem in the sims, and I'm starting to wonder how its possible to have been working on a project for over a year and still be so... wrong. Am I just going in circles? Each edit creates a new error? Am I sisyphus, doomed to my coding purgatory for all time? 

The issue now is with the confidence intervals and pvalues. As you said, confidence intervals that do not cross zero should be significant. Except we have a whole bunch that aren't significant that should be: 

### Covariance

Shown are estimates with 95% CIs from a single replicate. 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.2.CovConfint_discrep.png)

### GxE

This one has its own set of problems: in the top left, the damn GxE estimates don't overlap with the confidence intervals?! WTF. 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.2.GxEConfint_discreps.png)


### Estimates with CI's and True estimate

I already looked at whether confidence intervals overlap with true values. I found that no, they don't, but at least with covariance it is close enough it doesn't seem to be alarming. Here are the plots again, this time with "x" as the estimate (colors signify std. deviation), and red circle as the true value. 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.2.ConfintOverlap_true.png)

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/7.2.GxEoverlapper_true.png)

From that plot of GxE, it is clear that its just the first sets of parameters that generate such wonky CIs/estimates. Here is a subset of the parameter sets that generate those weird values:
One clear pattern is that these occur in the n_pop = 10 scenario. But why would that drive this pattern while other n_pop values don't?

```{r}
 delta_env delta_gen sample_size n_pop std_dev interaction true_GxE_emm GxE_emm GxE_emm_lwrCI GxE_emm_uprCI
0.01	-1.00	 5	10	0.5	0	0	0.05	0.06	0.08
0.50	-1.00	 5	10	0.5	0	0	0.05	0.06	0.07
1.00	-1.00	 5	10	0.5	0	0	0.04	0.04	0.06
0.01	0.01	 5	10	0.5	0	0	0.33	0.39	0.49
0.50	0.01	 5	10	0.5	0	0	0.10	0.12	0.16
1.00	0.01	 5	10	0.5	0	0	0.05	0.06	0.08
0.01	1.00	 5	10	0.5	0	0	0.05	0.06	0.08
0.50	1.00	 5	10	0.5	0	0	0.05	0.06	0.07
1.00	1.00	 5	10	0.5	0	0	0.04	0.04	0.06
0.01	-1.00	10	10	0.5	0	0	0.04	0.04	0.06
0.50	-1.00	10	10	0.5	0	0	0.04	0.04	0.05
1.00	-1.00	10	10	0.5	0	0	0.03	0.03	0.04
0.01	0.01	10	10	0.5	0	0	0.24	0.27	0.36
0.50	0.01	10	10	0.5	0	0	0.08	0.09	0.11
1.00	0.01	10	10	0.5	0	0	0.04	0.04	0.06
0.01	1.00	10	10	0.5	0	0	0.04	0.04	0.06
0.50	1.00	10	10	0.5	0	0	0.04	0.04	0.05
1.00	1.00	10	10	0.5	0	0	0.03	0.03	0.04
0.01	-1.00	 5	10	1.5	0	0	0.14	0.17	0.21
0.50	-1.00	 5	10	1.5	0	0	0.13	0.15	0.20
1.00	-1.00	 5	10	1.5	0	0	0.11	0.13	0.16
0.01	0.01	 5	10	1.5	0	0	0.33	0.39	0.49
0.50	0.01	 5	10	1.5	0	0	0.23	0.27	0.35
1.00	0.01	 5	10	1.5	0	0	0.15	0.17	0.22
0.01	1.00	 5	10	1.5	0	0	0.15	0.17	0.22
0.50	1.00	 5	10	1.5	0	0	0.13	0.16	0.20
1.00	1.00	 5	10	1.5	0	0	0.11	0.13	0.16
0.01	-1.00	10	10	1.5	0	0	0.11	0.12	0.16
0.50	-1.00	10	10	1.5	0	0	0.10	0.11	0.15
1.00	-1.00	10	10	1.5	0	0	0.08	0.09	0.12
0.01	0.01	10	10	1.5	0	0	0.24	0.27	0.36
0.50	0.01	10	10	1.5	0	0	0.17	0.19	0.26
1.00	0.01	10	10	1.5	0	0	0.11	0.12	0.16
0.01	1.00	10	10	1.5	0	0	0.11	0.12	0.16
0.50	1.00	10	10	1.5	0	0	0.10	0.11	0.15

```

I'll continue to think on this. It seems a bit excessive to be happening due to random chance. 
