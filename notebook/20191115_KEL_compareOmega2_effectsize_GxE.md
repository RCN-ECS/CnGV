Can you compare omega^2 with effect_size (the one calculated with this code, 
which I updated to fix the mistake https://github.com/RCN-ECS/CnGV/blob/master/notebook/20191104_KEL_interaction_effect.md) 
for the following phenotypic means for each genotype and each environment, and NO error? 

The link above describes the function:
`calc_interaction <- function(G11, G12, G21, G22)`
```
G11 <- 0.5 #Env 1, G1
G12 <- 1   #Env 2, G1
G21 <- 0   #Env 1, G2
G22 <- 0.5 #Env 2, G2
```

These are the two cases I think we should compare:
```
calc_interaction(-1, 1, 1, -1)
calc_interaction(-0.25, 0.25, 0.25, -0.25)
```
My prediction is that both should give omega^2 = 1 (if omega^2 is explaining a proportion of the variance, which I 
think it does), 
but only the effect size estimate will correctly return a larger effect size for the first case.
