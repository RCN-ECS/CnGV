# Example interaction effect
```
# Means of genotypes in each environmnent
G11 <- 0.5 #Env 1, G1
G12 <- 1   #Env 2, G1
G21 <- 0   #Env 1, G2
G22 <- 0.5 #Env 2, G2

calc_interaction <- function(G11, G12, G21, G22){
  
  plot(c(1,2), c(G11, G12), type="l", col="red", ylim=c(0,2))
  points(c(1,2), c(G21, G22), type="l", col="blue", ylim=c(0,2))
  
  
  (overall_mean <- mean(c(G11,G12,G21,G22)))
  
  # Marginal means
  (G1_mean <- overall_mean - mean(c(G11,G12)))
  (G2_mean <- overall_mean - mean(c(G22,G21)))
  (E1_mean <- overall_mean - mean(c(G11,G21)))
  (E2_mean <- overall_mean - mean(c(G12,G22)))
  
  ## Interaction effect for the i_th level of factor G (G2)
  ## and jth level of factor E (E2)
  print(c("Interaction effect for G22:",  
          overall_mean - G2_mean - E2_mean - G22))
  
  print(c("Interaction effect for G11:", 
          overall_mean - G1_mean - E1_mean - G11))
  
  print(c("Interaction effect for G12:",
          overall_mean - G1_mean - E2_mean - G12))
  
  print(c("Interaction effect for G21:", 
          overall_mean - G2_mean - E1_mean - G21))
}

calc_interaction(0.5, 1, 0, 0.5)
calc_interaction(0.5, 1, 0, 0.6)
calc_interaction(0.5, 1, 1, 0.5)

# whether it is positive or negative will depend on the constrast,
# suggest using absolute value for the paper
```
