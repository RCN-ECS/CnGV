#Another Day, Another Problem

The permutation code is not giving me accurate p-values for the GxE means values. 

I have tried everything I can think of (i've been debugging for 3 days straight) and I cannot find the problem. 

FOR EXAMPLE. Here is a 5 populations scenario that SHOULD have significant GxE across the board: 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/notebook_figs/609_meanProbplot.png)

however, when I run my code, here are my GxE results: 



GxE_emm_loop              0.65 # New way of estimating using for-loop
true_GxE_emm_loop         0.66
true_GxE_emm_loop_pvalue  0.00
GxE_emm_loop_lwrCI        0.64
GxE_emm_loop_uprCI        0.66
GxE_emm_loop_pvalue       0.00

true_omega_GxE            0.71 # Omega^2
true_omega_GxE_pvalue     0.00
GxE_omega_estimate        0.69
GxE_omega_lwrCI           0.68
GxE_omega_uprCI           0.71
GxE_omega_pvalue          0.00

true_eta_GxE              0.71 # Eta^2
true_eta_GxE_pvalue       0.00
GxE_eta                   0.70
GxE_eta_lwrCI             0.69
GxE_eta_uprCI             0.71
GxE_eta_pvalue            0.00

True_resid_variation      0.71 # Residual Variation
Resid_variation           0.70
GxE_residVar_lwrCI        0.69
GxE_residVar_uprCI        0.72

*true_GxE_means            0.65
true_GxE_means_pvalue     0.50
GxE_means                 0.64
GxE_means_lwrCI           0.63
GxE_means_uprCI           0.65
GxE_means_pvalue          0.50*
