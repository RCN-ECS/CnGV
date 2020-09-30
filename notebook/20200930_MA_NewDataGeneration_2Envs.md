# New Simulation approach

We simulated 2 scenarios. 1 in which the number of environments matched the number of genotypes, and the second, in which there were 2 environments but multiple genotypes. 

In the second scenario results, my covariance estimates were bound near 0.4 and -0.4 and never reached 1 or -1 as they should.

I figured out that the way we are simulating the data is confining the bounds. Specifically, delta_gen applied to each genotype, thus spacing out genotypes that shared an enviroment and limiting the ability of genotypes from the same environment to cluster together. 
