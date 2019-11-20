# Predicting CoGV and CnGV

Cogradient variation can be defined as positive covariance, while countergradient variation can be defined as negative covariance.

Specifying which scenario requires knowledge of each genotype's native environment. In our simulations, G1 is always native to E1, while G2 is always native to G2. 
Predicting whether CoGV or CnGV depends on both the slope and intercepts of each genotype. 

For example, each of the following are examples of *CoGV*. Native environment is denoted by the X. If slope is positive, then the intercept of G2 should be greater than G1. If the slope is negative, the intercept of G1 should be greater than G2. 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/plotA.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/PlotC.png)

Conversely, in these examples of *CnGV*, if the slope is positive, then the intercept of G1 should be greater than G2. If the slope is negative, the intercept of G2 should be greater than G1. 

![image](https://github.com/RCN-ECS/CnGV/blob/master/results/PlotB.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/plotD.png)

In pure GxE scenarios, we should expect the slopes of each genotype to be exactly opposite but have the same intercept: 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/PlotE.png)
