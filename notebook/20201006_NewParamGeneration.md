# New Code to estimate parameter space

In my previous attempts to generate starting parameters that cover the full parameter space, my method worked from the ground up. I started with 'expand.grid' for basic parameters and then added more on top of that grid. As you may imagine, it was difficult to limit the size of the resulting dataframe since every added variable would double or triple the size. 
